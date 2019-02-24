//
//  UniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniCellBKDTGridSBSolver.hpp"

//#include <omp.h>


template <template <class value_type, size_t, class, typename Point<value_type, 3>::DistType> class KDTType>
UniCellBKDTGridSBSolver<KDTType>::
UniCellBKDTGridSBSolver(double alpc, size_t maxCacheCellVecSize)
: UniLatLngBKDTGridSBSolver<KDTType>(alpc, maxCacheCellVecSize) {}


template <template <class value_type, size_t, class, typename Point<value_type, 3>::DistType> class KDTType>
void UniCellBKDTGridSBSolver<KDTType>::fillGridCache() {
    this->gridCache.reserve(this->locKdt.size()*1.2/this->AVE_LOC_PER_CELL);
    thisRowStartIdx.reserve(this->rowSize);
    std::vector<std::pair<Point<double, 3>, const SBLoc*>> ptLocPairs;
    ptLocPairs.reserve(this->MAX_CACHE_CELL_VEC_SIZE);
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    double thisLat = -0.5 * M_PI, thisCtrLat = thisLat + 0.5 * this->latInc;
    for (size_t r = 0, idx = 0; r < this->rowSize;
         ++r, thisCtrLat += this->latInc, thisLat += this->latInc) {
        size_t thisColSize = static_cast<size_t>(2*M_PI * SBLoc::EARTH_RADIUS *
                             cos(thisLat > 0 ? thisLat-this->latInc : thisLat)/
                             this->sideLen) + 2;
        double thisLngInc = 2*M_PI/thisColSize + 2*M_PI/(thisColSize*thisColSize*0xFFFF);
        thisRowStartIdx.emplace_back(idx, 1.0/thisLngInc);
        double thisDiff = SBLoc::xyzDistFromLngLat(r*this->latInc - 0.5*M_PI,
                          (r+1)*this->latInc - 0.5*M_PI, thisLngInc),
               thisCtrLng = 0.5 * thisLngInc - M_PI;
        for (size_t thisEndIdx = idx + thisColSize; idx < thisEndIdx;
             ++idx, thisCtrLng += thisLngInc) {
            UniLatLngBKDTGridSBSolver<KDTType>::fillCacheCell
            (thisCtrLng, thisCtrLat, thisDiff, ptLocPairs);
        }
    }
}

template <template <class value_type, size_t, class, typename Point<value_type, 3>::DistType> class KDTType>
const SBLoc* UniCellBKDTGridSBSolver<KDTType>::
findNearest(double lat, double lng) const {
    const auto &[startIdx, thisLngIncInverse] = thisRowStartIdx[(lat+0.5*M_PI)*this->latIncInverse];
    return UniLatLngBKDTGridSBSolver<KDTType>::returnNNLocFromCacheVariant(lat,
    lng, this->gridCache[startIdx + static_cast<size_t>((lng+M_PI)*thisLngIncInverse)]);
}




template class UniCellBKDTGridSBSolver<KDTree>;
template class UniCellBKDTGridSBSolver<KDTreeCusMem>;
template class UniCellBKDTGridSBSolver<KDTreeExpandLongest>;
template class UniCellBKDTGridSBSolver<KDTreeExpandLongestVec>;




/* TO DO
 
 
 Instead of using vector of pair of Point loc* cache,
 use a more compact pointer to pair pool cache, will likely reduce build time
 but not search time.
 */
