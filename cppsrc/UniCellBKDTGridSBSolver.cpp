//
//  UniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniCellBKDTGridSBSolver.hpp"

//#include <omp.h>


template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
UniCellBKDTGridSBSolver<KDTType, FPType, Policy>::
UniCellBKDTGridSBSolver(FPType alpc, std::size_t maxCacheCellVecSize)
: UniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>(alpc, maxCacheCellVecSize) {}


template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
void UniCellBKDTGridSBSolver<KDTType, FPType, Policy>::FillGridCache() {
    thisRowStartIdx.reserve(this->row_size_);
    this->grid_cache_.reserve(this->loc_kdt_.size()*1.2/this->AVE_LOC_PER_CELL);
    std::vector<typename KDT<KDTType, FPType>::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(pt_loc_vec) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    FPType thisLat = -0.5 * def::kMathPi<FPType>, thisCtrLat = thisLat + 0.5 * this->lat_inc_;
    for (std::size_t r = 0, idx = 0; r < this->row_size_;
         ++r, thisCtrLat += this->lat_inc_, thisLat += this->lat_inc_) {
        std::size_t thisColSize = static_cast<std::size_t>(2*def::kMathPi<FPType> * SBLoc<FPType>::EARTH_RADIUS *
                             std::cos(thisLat > 0 ? thisLat-this->lat_inc_ : thisLat)/
                             this->side_len_) + 2;
        FPType thisLngInc = 2*def::kMathPi<FPType>/thisColSize + 2*def::kMathPi<FPType>/(thisColSize*thisColSize*65536);
        thisRowStartIdx.emplace_back(idx, 1.0/thisLngInc);
        FPType lat1 = r*this->lat_inc_ - 0.5*def::kMathPi<FPType>;
        FPType diagonalDistSq3DEUC = SBLoc<FPType>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + this->lat_inc_, thisLngInc),
               thisCtrLng = 0.5 * thisLngInc - def::kMathPi<FPType>;
        for (std::size_t thisEndIdx = idx + thisColSize; idx < thisEndIdx;
             ++idx, thisCtrLng += thisLngInc) {
            UniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::FillCacheCell
            (thisCtrLng, thisCtrLat, diagonalDistSq3DEUC, pt_loc_vec);
        }
    }
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
const SBLoc<FPType>* UniCellBKDTGridSBSolver<KDTType, FPType, Policy>::
FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    const auto &[startIdx, thisLngIncInverse] = thisRowStartIdx[(geo_search_pt[0]+0.5*def::kMathPi<FPType>)*this->lat_inc_inverse_];
    return UniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::ReturnNNLocFromCacheVariant(geo_search_pt,
        this->grid_cache_[startIdx + static_cast<std::size_t>((geo_search_pt[1]+def::kMathPi<FPType>)*thisLngIncInverse)]);
}



template class UniCellBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kSingle>;
template class UniCellBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kSingle>;

template class UniCellBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSingle>;
template class UniCellBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSingle>;

template class UniCellBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSingle>;
template class UniCellBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSingle>;

template class UniCellBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSingle>;
template class UniCellBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSingle>;




/* TO DO
 
 
 Instead of using vector of pair of PointND loc* cache,
 use a more compact pointer to pair pool cache, will likely reduce build time
 but not search time.
 */
