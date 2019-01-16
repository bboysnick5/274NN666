//
//  UniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniLatLngBKDTGridSBSolver.hpp"
//#include <omp.h>

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
UniLatLngBKDTGridSBSolver<KDTType>::
UniLatLngBKDTGridSBSolver(double alpc, size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType>(), AVE_LOC_PER_CELL(alpc),
  MAX_CACHE_CELL_VEC_SIZE(maxCacheCellVecSize) {}


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UniLatLngBKDTGridSBSolver<KDTType>::printSolverInfo() const {
    std::cout << "Total tree nodes: " << totalNodeSize
    << "\nAve tree size: " << totalNodeSize/gridCache.size()
    << "\nAve tree height: "
    << static_cast<size_t>(log2(totalNodeSize/gridCache.size() + 1)) + 1
    << "\nRatio of tree nodes over num locs: "
    << totalNodeSize/this->locKdt.size()
    << "\nTotal num of loc cells: " << gridCache.size()
    << "\nSingle loc cells: " << singleLocs
    << "\nVector loc cells: " << vecLocs
    << "\nKd-tree loc cells: " << gridCache.size() -singleLocs -vecLocs << "\n";
}


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UniLatLngBKDTGridSBSolver<KDTType>::
fillCacheCell(double thisCtrLng, double thisCtrLat, double thisDiff,
              std::vector<std::pair<Point<3>, const SBLoc*>>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc::latLngToCart3DPt(thisCtrLng, thisCtrLat),
                                   thisDiff, std::back_inserter(ptLocPairs));
    size_t locsSize = ptLocPairs.size();
    if (locsSize == 1) {
        this->gridCache.emplace_back(ptLocPairs[0].second);
        this->singleLocs++;
    } else if (locsSize < MAX_CACHE_CELL_VEC_SIZE) {
        ptLocPairs.shrink_to_fit();
        this->gridCache.emplace_back(std::move(ptLocPairs));
        ptLocPairs.reserve(std::sqrt(this->locKdt.size()));
        this->vecLocs++;
    } else {
        this->gridCache.emplace_back(std::in_place_type<KDT<KDTType>>,
                                     ptLocPairs.begin(), ptLocPairs.end());
    }
    ptLocPairs.clear();
    this->totalNodeSize += locsSize;
}


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UniLatLngBKDTGridSBSolver<KDTType>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2*M_PI/colSize + 2*M_PI/(colSize*colSize*0xFFFF);
    gridCache.reserve(this->locKdt.size()*1.2/AVE_LOC_PER_CELL);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs(this->locKdt.size());
    
    double thisCtrLat = 0.5 * (latInc - M_PI);
    for (size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc) {
        double thisCtrLng = 0.5 * lngInc - M_PI;
        double thisDiff = SBLoc::xyzDistFromLngLat(r*latInc- 0.5*M_PI,
                          (r+1)*latInc-0.5*M_PI, lngInc);
        for (size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            fillCacheCell(thisCtrLng, thisCtrLat, thisDiff, ptLocPairs);
        }
    }
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UniLatLngBKDTGridSBSolver<KDTType>::calcSideLenFromAlpc() {
    double surfaceArea = 4*M_PI*SBLoc::EARTH_RADIUS*SBLoc::EARTH_RADIUS;
    double numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UniLatLngBKDTGridSBSolver<KDTType>::
build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    BKDTSBSolver<KDTType>::generateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc::latFromHavDist(sideLen, 0));
    rowSize = std::ceil(M_PI/(latInc - latInc*latInc/(M_PI*0xFFFF)));
    fillGridCache();
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
const SBLoc* UniLatLngBKDTGridSBSolver<KDTType>::
findNearest(double lng, double lat) const {
    const auto &v = gridCache[static_cast<size_t>((lat+0.5*M_PI)/latInc)*colSize
                              + static_cast<size_t>((lng+M_PI)/lngInc)];
    switch (v.index()) {
        case 0: {
            const auto p = SBLoc::latLngToCart3DPt(lng, lat);
            return std::min_element(std::get<0>(v).begin(), std::get<0>(v).end(),
                                    [&](const auto& p1, const auto& p2){
                                        return Point<3>::template dist<Point<3>::DistType::EUCSQ>(p1.first, p) < Point<3>::template dist<Point<3>::DistType::EUCSQ>(p2.first, p);
                                    })->second;
        }
        case 1:
            return std::get<1>(v);
        default:
            return std::get<2>(v).kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
    }
}





template class UniLatLngBKDTGridSBSolver<KDTree>;
template class UniLatLngBKDTGridSBSolver<KDTreeCusMem>;


