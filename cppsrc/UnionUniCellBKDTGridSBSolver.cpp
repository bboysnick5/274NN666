//
//  UnionUniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/20/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniCellBKDTGridSBSolver.hpp"
#include <thread>
#include <map>
#include <algorithm>
#include <omp.h>


template <template <typename FPType, std::uint_fast8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
UnionUniCellBKDTGridSBSolver<KDTType, FPType, policy>::
UnionUniCellBKDTGridSBSolver(FPType alpc, std::size_t maxCacheCellVecSize)
: UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>(alpc, maxCacheCellVecSize) {}


template <template <typename FPType, std::uint_fast8_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
void UnionUniCellBKDTGridSBSolver<KDTType, FPType, policy>::
LoopBody(def::Policy_Tag<def::ThreadingPolicy::kSingle>) {
    std::vector<typename KDT<KDTType, FPType>::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
    this->grid_cache_.reserve(totalCacheCells);
    FPType thisCtrLat = 0.5*this->lat_inc_ - 0.5*def::kMathPi<FPType>;
    FPType lat1 = - 0.5*def::kMathPi<FPType>;
    for (std::size_t r = 0; r < this->row_size_; ++r, thisCtrLat += this->lat_inc_, lat1 += this->lat_inc_) {
        auto &[thisColSize, cosThisLngInc] = col_size_CosLngIncEachRowVec[r];
        FPType thisLngInc = 1.0/thisRowStartIdxThisLngIncInverseVec[r].second;
        FPType diagonalDistSq3DEUC = UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::EUC3DDistSqFromLatCosDeltaLng(lat1, lat1 + this->lat_inc_, cosThisLngInc);
        FPType thisCtrLng = 0.5 * thisLngInc - def::kMathPi<FPType>;
        for (std::size_t c = 0; c < thisColSize; ++c, thisCtrLng += thisLngInc) {
            this->FillCacheCell({thisCtrLat, thisCtrLng}, diagonalDistSq3DEUC, thisColSize, pt_loc_vec);
        }
    }
}

template <template <typename FPType, std::uint_fast8_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
void UnionUniCellBKDTGridSBSolver<KDTType, FPType, policy>::
LoopBody(def::Policy_Tag<def::ThreadingPolicy::kMultiOmp>) {
    this->grid_cache_.resize(totalCacheCells, 0);
    FPType initCtrLat = 0.5*this->lat_inc_ - 0.5*def::kMathPi<FPType>;
    FPType initLat1 = - 0.5*def::kMathPi<FPType>;
    std::vector<typename KDT<KDTType, FPType>::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
firstprivate(pt_loc_vec) default(shared) schedule(dynamic, 1) 
    for (std::size_t idx = 0; idx < totalCacheCells; ++idx) {
        const auto it = std::upper_bound(thisRowStartIdxThisLngIncInverseVec.cbegin(), thisRowStartIdxThisLngIncInverseVec.cend(), idx, [](const std::size_t &idx, const std::pair<std::size_t, FPType> &p){return idx < p.first;}) - 1;
        std::size_t r = it - thisRowStartIdxThisLngIncInverseVec.begin();
        std::size_t c = idx - it->first;
        FPType thisLngInc = 1.0/it->second;
        auto &[thisColSize, cosThisLngInc] = col_size_CosLngIncEachRowVec[r];
        FPType lat1 = r*this->lat_inc_ + initLat1;
        FPType thisCtrLat = initCtrLat + r*this->lat_inc_;
        FPType diagonalDistSq3DEUC = UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::EUC3DDistSqFromLatCosDeltaLng(lat1, lat1 + this->lat_inc_, cosThisLngInc);
        FPType initThisCtrLng = 0.5 * thisLngInc - def::kMathPi<FPType>;
        this->FillCacheCell(idx, {thisCtrLat, initThisCtrLng + c*thisLngInc}, diagonalDistSq3DEUC, thisColSize, pt_loc_vec);
    }
}

template <template <typename FPType, std::uint_fast8_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
void UnionUniCellBKDTGridSBSolver<KDTType, FPType, policy>::FillGridCache() {
    col_size_CosLngIncEachRowVec.reserve(this->row_size_);
    thisRowStartIdxThisLngIncInverseVec.reserve(this->row_size_);
    FPType earthPerimeterOverSideLen = 2.0*def::kMathPi<FPType> * SBLoc<FPType>::EARTH_RADIUS /this->side_len_;
    FPType thisLat = -0.5*def::kMathPi<FPType>;
    totalCacheCells = 0;
    for (std::size_t r = 0; r < this->row_size_ - 1; ++r, thisLat += this->lat_inc_) {
        std::size_t thisColSize = static_cast<std::size_t>(earthPerimeterOverSideLen * std::cos(thisLat)) + 1;
        FPType thisLngInc = 2.0*def::kMathPi<FPType>/thisColSize + std::numeric_limits<FPType>::epsilon();
        col_size_CosLngIncEachRowVec.emplace_back(thisColSize, std::cos(thisLngInc));
        thisRowStartIdxThisLngIncInverseVec.emplace_back(totalCacheCells, 1.0/thisLngInc);
        totalCacheCells += thisColSize;
    }
    col_size_CosLngIncEachRowVec.emplace_back(col_size_CosLngIncEachRowVec.back());
    thisRowStartIdxThisLngIncInverseVec.emplace_back(totalCacheCells, thisRowStartIdxThisLngIncInverseVec.back().second);
    totalCacheCells += col_size_CosLngIncEachRowVec.back().first;
    UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::LoopBodyThreadingPolicyDispatch();
    col_size_CosLngIncEachRowVec.clear();
    col_size_CosLngIncEachRowVec.shrink_to_fit();
}

template <template <typename FPType, std::uint_fast8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
const SBLoc<FPType>* UnionUniCellBKDTGridSBSolver<KDTType, FPType, policy>::
FindNearestLoc(PointND<FPType, 2> geo_search_pt) const {
  //  return UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::FindNearestLoc(geoPt);
    const auto &[startIdx, thisLngIncInverse] = thisRowStartIdxThisLngIncInverseVec[(geo_search_pt[0]+0.5*def::kMathPi<FPType>)*this->lat_inc_inverse_];
    return UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::ReturnNNLocFromCacheVariant(geo_search_pt,
        this->grid_cache_[startIdx + static_cast<std::size_t>((geo_search_pt[1]+def::kMathPi<FPType>)*thisLngIncInverse)]);
}


template class UnionUniCellBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniCellBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kMultiHand>;

template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiHand>;

template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiHand>;

template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiHand>;
