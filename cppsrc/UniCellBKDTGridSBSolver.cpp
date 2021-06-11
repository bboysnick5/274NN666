//
//  UniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniCellBKDTGridSBSolver.hpp"

//#include <omp.h>


template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type>
UniCellBKDTGridSBSolver<KDTType, dist_type>::
UniCellBKDTGridSBSolver(dist_type alpc, std::size_t maxCacheCellVecSize)
: UniLatLngBKDTGridSBSolver<KDTType, dist_type>(alpc, maxCacheCellVecSize) {}


template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type>
void UniCellBKDTGridSBSolver<KDTType, dist_type>::FillGridCache() {
    thisRowStartIdx.reserve(this->rowSize);
    this->grid_cache_.reserve(this->locKdt.size()*1.2/this->AVE_LOC_PER_CELL);
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(this->MAX_CACHE_CELL_VEC_SIZE);
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    dist_type thisLat = -0.5 * def::kMathPi<dist_type>, thisCtrLat = thisLat + 0.5 * this->latInc;
    for (std::size_t r = 0, idx = 0; r < this->rowSize;
         ++r, thisCtrLat += this->latInc, thisLat += this->latInc) {
        std::size_t thisColSize = static_cast<std::size_t>(2*def::kMathPi<dist_type> * SBLoc<dist_type>::EARTH_RADIUS *
                             std::cos(thisLat > 0 ? thisLat-this->latInc : thisLat)/
                             this->sideLen) + 2;
        dist_type thisLngInc = 2*def::kMathPi<dist_type>/thisColSize + 2*def::kMathPi<dist_type>/(thisColSize*thisColSize*65536);
        thisRowStartIdx.emplace_back(idx, 1.0/thisLngInc);
        dist_type lat1 = r*this->latInc - 0.5*def::kMathPi<dist_type>;
        dist_type diagonalDistSq3DEUC = SBLoc<dist_type>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + this->latInc, thisLngInc),
               thisCtrLng = 0.5 * thisLngInc - def::kMathPi<dist_type>;
        for (std::size_t thisEndIdx = idx + thisColSize; idx < thisEndIdx;
             ++idx, thisCtrLng += thisLngInc) {
            UniLatLngBKDTGridSBSolver<KDTType, dist_type>::FillCacheCell
            (thisCtrLng, thisCtrLat, diagonalDistSq3DEUC, ptLocPairs);
        }
    }
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UniCellBKDTGridSBSolver<KDTType, dist_type>::
FindNearestLoc(const PointND<dist_type, 2>& geoSearchPt) const {
    const auto &[startIdx, thisLngIncInverse] = thisRowStartIdx[(geoSearchPt[0]+0.5*def::kMathPi<dist_type>)*this->latIncInverse];
    return UniLatLngBKDTGridSBSolver<KDTType, dist_type>::ReturnNNLocFromCacheVariant(geoSearchPt,
        this->grid_cache_[startIdx + static_cast<std::size_t>((geoSearchPt[1]+def::kMathPi<dist_type>)*thisLngIncInverse)]);
}



template class UniCellBKDTGridSBSolver<KDTree, double>;
template class UniCellBKDTGridSBSolver<KDTree, float>;

template class UniCellBKDTGridSBSolver<KDTreeCusMem, double>;
template class UniCellBKDTGridSBSolver<KDTreeCusMem, float>;

template class UniCellBKDTGridSBSolver<KDTreeExpandLongest, double>;
template class UniCellBKDTGridSBSolver<KDTreeExpandLongest, float>;

template class UniCellBKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class UniCellBKDTGridSBSolver<KDTreeExpandLongestVec, float>;




/* TO DO
 
 
 Instead of using vector of pair of PointND loc* cache,
 use a more compact pointer to pair pool cache, will likely reduce build time
 but not search time.
 */
