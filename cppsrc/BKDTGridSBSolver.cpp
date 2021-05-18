//
//  BKDTGridSBSolver<KDTType, dist_type>.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/23/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//


/*
 Build Version 1
 Inaccurate and waste of resource.
 Precalculate the side length in HavDist. Then apply that to each cell.
 */

#include <thread>
//#include <omp.h>
#include "BKDTGridSBSolver.hpp"

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
BKDTGridSBSolver<KDTType, dist_type>::BKDTGridSBSolver(dist_type alpc) : GridSBSolver<dist_type>(alpc) {}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
typename std::vector<typename KDTree<dist_type, 3, const SBLoc<dist_type>*, Point<dist_type, 3>::DistType::EUC>::node_type>::iterator
BKDTGridSBSolver<KDTType, dist_type>::cacheAllPossibleLocsOneCell(size_t r0, size_t c0, dist_type diff,
                                                                  typename std::vector<typename KDTree<dist_type, 3, const SBLoc<dist_type>*, Point<dist_type, 3>::DistType::EUC>::node_type>::iterator begin) {
    dist_type rowDistCellCtrGridCtr = ((r0+ 0.5 - this->rowSize/2 )*this->sideLen),
           colDistCellCtrGridCtr = ((c0 + 0.5 -this->colSize/2)*this->sideLen);
    dist_type cellCtrLat = this->midLat + SBLoc<dist_type>::deltaLatOnSameLngFromHavDist(rowDistCellCtrGridCtr),
           cellCtrLng = SBLoc<dist_type>::lngFromHavDist(colDistCellCtrGridCtr,
                                              this->midLng, cellCtrLat);
    return sbKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt({cellCtrLat,
        cellCtrLng}), diff, begin);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void BKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
    gridTreeCache = std::vector<KDTree<dist_type, 3, const SBLoc<dist_type>*, Point<dist_type, 3>::DistType::EUC>>
                    (this->rowSize*this->colSize);
    gridSingleCache = std::vector<const SBLoc<dist_type>*>(this->rowSize*this->colSize);
    std::vector<typename KDTree<dist_type, 3, const SBLoc<dist_type>*, Point<dist_type, 3>::DistType::EUC>::node_type> ptLocPairs(this->numLocs);
    size_t totalTreeSize = 0, singleLocs = 0;
    dist_type diff = xyzDistFromSideLen();
//#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
//default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
//reduction(+:totalTreeSize, singleLocs) collapse(2)
    for (size_t r = 0; r < this->rowSize; ++r) {
        for (size_t c = 0; c < this->colSize; ++c) {
            size_t idx = r*this->colSize+c;
            auto locsEnd = cacheAllPossibleLocsOneCell(r, c, diff,
                                                       ptLocPairs.begin());
            size_t locsSize = locsEnd - ptLocPairs.begin();
            if (locsSize > 1) {
                gridTreeCache[idx] = KDTree<dist_type, 3, const SBLoc<dist_type>*, Point<dist_type, 3>::DistType::EUC>
                    (ptLocPairs.begin(), locsEnd);
            } else {
                gridSingleCache[idx] = ptLocPairs[0].value;
                singleLocs++;
            }
            totalTreeSize += locsSize;
        }
    }
    size_t multiLocs = this->rowSize*this->colSize - singleLocs;
    std::cout << "ave tree size: " << totalTreeSize/(this->rowSize*this->colSize)
              << "\nSingle loc cells: " << singleLocs
              << "\nMulti-loc cells:" << multiLocs << std::endl;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
dist_type BKDTGridSBSolver<KDTType, dist_type>::xyzDistFromSideLen() {
    dist_type lat2 = SBLoc<dist_type>::deltaLatOnSameLngFromHavDist(this->sideLen*sqrt(2));
    return Point<dist_type, 3>::template
    dist<Point<dist_type, 3>::DistType::EUC>(SBLoc<dist_type>::geoPtToCart3DPt({0.0, 0.0}),
                                             SBLoc<dist_type>::geoPtToCart3DPt({lat2, 0.0}));
}
 
template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void BKDTGridSBSolver<KDTType, dist_type>::build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    std::vector<typename KDTree<dist_type, 3, const SBLoc<dist_type>*, Point<dist_type, 3>::DistType::EUC>::node_type> kdtData;
    kdtData.reserve(locData->size());
    std::transform(locData->begin(), locData->end(), std::back_inserter(kdtData),
                   [&](const SBLoc<dist_type>& loc) -> typename KDTree<dist_type, 3, const SBLoc<dist_type>*, Point<dist_type, 3>::DistType::EUC>::node_type { return
                       {SBLoc<dist_type>::geoPtToCart3DPt(loc.geoPt), &loc};});
    sbKdt = KDTree<dist_type, 3, const SBLoc<dist_type>*, Point<dist_type, 3>::DistType::EUC>(kdtData.begin(), kdtData.end());
    this->numLocs = sbKdt.size();
    this->findKeyLngLat(locData);
    this->rowSize = sqrt(locData->size()) / this->AVE_LOC_PER_CELL;
    this->sideLen = 
    SBLoc<dist_type>::havDist({this->minLat, 0.0}, {this->maxLat, 0.0}) / this->rowSize;
    dist_type lowestLatCircleRadius = SBLoc<dist_type>::EARTH_RADIUS * cos(this->minLat < 0 ? 0 : this->minLat);
    dist_type longestColDistSpan = 2 * std::numbers::pi_v<dist_type> * lowestLatCircleRadius *
                                (std::fabs(this->maxLng - this->minLng)/(2*std::numbers::pi_v<dist_type>));
    this->colSize = longestColDistSpan/this->sideLen + 1;
    
    fillGridCache();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* BKDTGridSBSolver<KDTType, dist_type>::findNearest(const Point<dist_type, 2>& geoSearchPt) const {
    auto idxPr = this->getIdx(geoSearchPt[1], geoSearchPt[0]);
    size_t idx = idxPr.first*this->colSize+idxPr.second;
    auto singleLoc = gridSingleCache[idx];
    return singleLoc != nullptr ? singleLoc :
           gridTreeCache[idx].kNNValue(SBLoc<dist_type>::geoPtToCart3DPt(geoSearchPt), 1);
}


template class BKDTGridSBSolver<KDTree, double>;
template class BKDTGridSBSolver<KDTree, float>;

/*
template class BKDTGridSBSolver<KDTreeCusMem, double>;
template class BKDTGridSBSolver<KDTreeCusMem, float>;

template class BKDTGridSBSolver<KDTreeExpandLongest, double>;
template class BKDTGridSBSolver<KDTreeExpandLongest, float>;

template class BKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class BKDTGridSBSolver<KDTreeExpandLongestVec, float>;
*/



