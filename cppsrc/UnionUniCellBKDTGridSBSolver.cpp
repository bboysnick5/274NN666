//
//  UnionUniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/20/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniCellBKDTGridSBSolver.hpp"



template <template <size_t, class, typename Point<3>::DistType> class KDTType>
UnionUniCellBKDTGridSBSolver<KDTType>::
UnionUniCellBKDTGridSBSolver(double alpc, size_t maxCacheCellVecSize)
: UnionUniLatLngBKDTGridSBSolver<KDTType>(alpc, maxCacheCellVecSize) {}


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UnionUniCellBKDTGridSBSolver<KDTType>::fillGridCache() {
    this->gridCache.reserve(this->locKdt.size()*1.2/this->AVE_LOC_PER_CELL);
    thisRowStartIdx.reserve(this->rowSize);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs;
    ptLocPairs.reserve(this->MAX_CACHE_CELL_VEC_SIZE);
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    double thisLat = -0.5 * M_PI, thisCtrLat = thisLat + 0.5 * this->latInc;
    for (size_t r = 0, idx = 0; r < this->rowSize;
         ++r, thisCtrLat += this->latInc, thisLat += this->latInc) {
        size_t thisColSize = static_cast<size_t>(2.0*M_PI * SBLoc::EARTH_RADIUS *
                             cos(thisLat > 0.0 ? thisLat-this->latInc : thisLat)/
                             this->sideLen) + 2;
        double thisColSizeInverse = 1.0/thisColSize;
        double rawThisLngInc = 2.0*M_PI*thisColSizeInverse;
        double thisLngInc = rawThisLngInc + rawThisLngInc*thisColSizeInverse/0xFFFF;
        thisRowStartIdx.emplace_back(idx, 1.0/thisLngInc);
        double thisDiff = SBLoc::xyzDistFromLngLat(r*this->latInc - 0.5*M_PI,
                          (r+1)*this->latInc - 0.5*M_PI, thisLngInc),
               thisCtrLng = 0.5 * thisLngInc - M_PI;
        for (size_t thisEndIdx = idx + thisColSize; idx < thisEndIdx;
             ++idx, thisCtrLng += thisLngInc) {
            UnionUniLatLngBKDTGridSBSolver<KDTType>::fillCacheCell
            (thisCtrLng, thisCtrLat, thisDiff, ptLocPairs);
        }
    }
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
const SBLoc* UnionUniCellBKDTGridSBSolver<KDTType>::
findNearest(double lat, double lng) const {
    const auto &[startIdx, thisLngIncInverse] = thisRowStartIdx[(lat+0.5*M_PI)*this->latIncInverse];
    return UnionUniLatLngBKDTGridSBSolver<KDTType>::returnNNLocFromCacheVariant(lat,
        lng, this->gridCache[startIdx + static_cast<size_t>((lng+M_PI)*thisLngIncInverse)]);
}




template class UnionUniCellBKDTGridSBSolver<KDTree>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec>;

