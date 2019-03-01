//
//  UnionUniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/20/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniCellBKDTGridSBSolver.hpp"



template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniCellBKDTGridSBSolver<KDTType, dist_type>::
UnionUniCellBKDTGridSBSolver(dist_type alpc, size_t maxCacheCellVecSize)
: UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>(alpc, maxCacheCellVecSize) {}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniCellBKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
    thisRowStartIdx.reserve(this->rowSize);
    this->gridCache.reserve(this->locKdt.size()*1.2/this->AVE_LOC_PER_CELL);
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(this->MAX_CACHE_CELL_VEC_SIZE);
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    dist_type thisLat = -0.5 * M_PI, thisCtrLat = thisLat + 0.5 * this->latInc;
    for (size_t r = 0, idx = 0; r < this->rowSize;
         ++r, thisCtrLat += this->latInc, thisLat += this->latInc) {
        size_t thisColSize = static_cast<size_t>(2.0*M_PI * SBLoc<dist_type>::EARTH_RADIUS *
                             cos(thisLat > 0.0 ? thisLat-this->latInc : thisLat)/
                             this->sideLen) + 2;
        dist_type thisColSizeInverse = 1.0/thisColSize;
        dist_type rawThisLngInc = 2.0*M_PI*thisColSizeInverse;
        dist_type thisLngInc = rawThisLngInc + rawThisLngInc*thisColSizeInverse/0xFFFF;
        thisRowStartIdx.emplace_back(idx, 1.0/thisLngInc);
        dist_type thisDiff = SBLoc<dist_type>::xyzDistFromLngLat(r*this->latInc - 0.5*M_PI,
                             (r+1)*this->latInc - 0.5*M_PI, thisLngInc),
                  thisCtrLng = 0.5 * thisLngInc - M_PI;
        for (size_t thisEndIdx = idx + thisColSize; idx < thisEndIdx;
             ++idx, thisCtrLng += thisLngInc) {
            UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::fillCacheCell
            ({thisCtrLat, thisCtrLng}, thisDiff, thisColSize, ptLocPairs);
        }
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UnionUniCellBKDTGridSBSolver<KDTType, dist_type>::
findNearest(const Point<dist_type, 2>& geoPt) const {
    //const auto &[lat, lng] = geoPt;
    const auto &[startIdx, thisLngIncInverse] = thisRowStartIdx[(geoPt[0]+0.5*M_PI)*this->latIncInverse];
    return UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::returnNNLocFromCacheVariant(geoPt,
        this->gridCache[startIdx + static_cast<size_t>((geoPt[1]+M_PI)*thisLngIncInverse)]);
}



template class UnionUniCellBKDTGridSBSolver<KDTree, double>;
template class UnionUniCellBKDTGridSBSolver<KDTree, float>;

template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, double>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, float>;

template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, double>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, float>;

template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, float>;



