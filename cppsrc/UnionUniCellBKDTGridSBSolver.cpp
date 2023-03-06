//
//  UnionUniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/20/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniCellBKDTGridSBSolver.hpp"

#include <omp.h>

#include <algorithm>
#include <map>
#include <thread>

template <SolverConfig Config>
UnionUniCellBKDTGridSBSolver<Config>::UnionUniCellBKDTGridSBSolver(
    FPType alpc, std::size_t maxCacheCellVecSize)
    : UnionUniLatLngBKDTGridSBSolver<Config>(alpc, maxCacheCellVecSize) {}

template <SolverConfig Config>
void UnionUniCellBKDTGridSBSolver<Config>::LoopBody(
    def::ThreadingPolicyTag<def::ThreadingPolicy::kSingle>) {
    std::vector<typename KDTType::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
    this->grid_cache_.reserve(totalCacheCells);
    FPType thisCtrLat = 0.5 * this->lat_inc_ - 0.5 * def::kMathPi<FPType>;
    FPType lat1 = -0.5 * def::kMathPi<FPType>;
    for (std::size_t r = 0; r < this->row_size_;
         ++r, thisCtrLat += this->lat_inc_, lat1 += this->lat_inc_) {
        auto &[thisColSize, cosThisLngInc] = col_size_CosLngIncEachRowVec[r];
        FPType thisLngInc = 1.0 / thisRowStartIdxThisLngIncInverseVec[r].second;
        FPType diagonalDistSq3DCART = UnionUniLatLngBKDTGridSBSolver<
            Config>::CART3DDistSqFromLatCosDeltaLng(lat1, lat1 + this->lat_inc_,
                                                    cosThisLngInc);
        FPType thisCtrLng = 0.5 * thisLngInc - def::kMathPi<FPType>;
        for (std::size_t c = 0; c < thisColSize;
             ++c, thisCtrLng += thisLngInc) {
            this->FillCacheCell({thisCtrLat, thisCtrLng}, diagonalDistSq3DCART,
                                thisColSize, pt_loc_vec);
        }
    }
}

template <SolverConfig Config>
void UnionUniCellBKDTGridSBSolver<Config>::LoopBody(
    def::ThreadingPolicyTag<def::ThreadingPolicy::kMultiOmp>) {
    this->grid_cache_.resize(totalCacheCells, 0);
    FPType initCtrLat = 0.5 * this->lat_inc_ - 0.5 * def::kMathPi<FPType>;
    FPType initLat1 = -0.5 * def::kMathPi<FPType>;
    std::vector<typename KDTType::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
    firstprivate(pt_loc_vec) default(shared) schedule(dynamic, 1)
    for (std::size_t idx = 0; idx < totalCacheCells; ++idx) {
        const auto it =
            std::upper_bound(
                thisRowStartIdxThisLngIncInverseVec.cbegin(),
                thisRowStartIdxThisLngIncInverseVec.cend(), idx,
                [](std::size_t idx, const std::pair<std::size_t, FPType> &p) {
                    return idx < p.first;
                }) -
            1;
        std::size_t r = it - thisRowStartIdxThisLngIncInverseVec.begin();
        std::size_t c = idx - it->first;
        FPType thisLngInc = 1.0 / it->second;
        auto &[thisColSize, cosThisLngInc] = col_size_CosLngIncEachRowVec[r];
        FPType lat1 = r * this->lat_inc_ + initLat1;
        FPType thisCtrLat = initCtrLat + r * this->lat_inc_;
        FPType diagonalDistSq3DCART = UnionUniLatLngBKDTGridSBSolver<
            Config>::CART3DDistSqFromLatCosDeltaLng(lat1, lat1 + this->lat_inc_,
                                                    cosThisLngInc);
        FPType initThisCtrLng = 0.5 * thisLngInc - def::kMathPi<FPType>;
        this->FillCacheCell(idx, {thisCtrLat, initThisCtrLng + c * thisLngInc},
                            diagonalDistSq3DCART, thisColSize, pt_loc_vec);
    }
}

template <SolverConfig Config>
void UnionUniCellBKDTGridSBSolver<Config>::FillGridCache() {
    col_size_CosLngIncEachRowVec.reserve(this->row_size_);
    thisRowStartIdxThisLngIncInverseVec.reserve(this->row_size_);
    FPType earthPerimeterOverSideLen = 2.0 * def::kMathPi<FPType> *
                                       SBLoc<FPType>::EARTH_RADIUS /
                                       this->side_len_;
    FPType thisLat = -0.5 * def::kMathPi<FPType>;
    totalCacheCells = 0;
    for (std::size_t r = 0; r < this->row_size_ - 1;
         ++r, thisLat += this->lat_inc_) {
        std::size_t thisColSize =
            static_cast<std::size_t>(earthPerimeterOverSideLen *
                                     std::cos(thisLat)) +
            1;
        FPType thisLngInc = 2.0 * def::kMathPi<FPType> / thisColSize +
                            std::numeric_limits<FPType>::epsilon();
        col_size_CosLngIncEachRowVec.emplace_back(thisColSize,
                                                  std::cos(thisLngInc));
        thisRowStartIdxThisLngIncInverseVec.emplace_back(totalCacheCells,
                                                         1.0 / thisLngInc);
        totalCacheCells += thisColSize;
    }
    col_size_CosLngIncEachRowVec.emplace_back(
        col_size_CosLngIncEachRowVec.back());
    thisRowStartIdxThisLngIncInverseVec.emplace_back(
        totalCacheCells, thisRowStartIdxThisLngIncInverseVec.back().second);
    totalCacheCells += col_size_CosLngIncEachRowVec.back().first;
    UnionUniLatLngBKDTGridSBSolver<Config>::LoopBodyThreadingPolicyDispatch();
    col_size_CosLngIncEachRowVec.clear();
    col_size_CosLngIncEachRowVec.shrink_to_fit();
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *UnionUniCellBKDTGridSBSolver<Config>::FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    //  return UnionUniLatLngBKDTGridSBSolver<KDTType, FPType,
    //  ParPolicy>::FindNearestLoc(geo_pt);
    const auto &[startIdx, thisLngIncInverse] =
        thisRowStartIdxThisLngIncInverseVec[(geo_search_pt[0] +
                                             0.5 * def::kMathPi<FPType>)*this
                                                ->lat_inc_inverse_];
    return UnionUniLatLngBKDTGridSBSolver<Config>::ReturnNNLocFromCacheVariant(
        geo_search_pt,
        this->grid_cache_[startIdx +
                          static_cast<std::size_t>(
                              (geo_search_pt[1] +
                               def::kMathPi<FPType>)*thisLngIncInverse)]);
}

template class UnionUniCellBKDTGridSBSolver<kConfigStSoaNoSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigMtOmpSoaNoSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStAosNoSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigMtOmpAosNoSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStAosoaNoSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStSoaSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigMtOmpSoaSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStAosSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigMtOmpAosSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStAosoaSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<kConfigMtOmpAosoaSimdKDTree<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeCusMem<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongest<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongestVec<float>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<float>>;

template class UnionUniCellBKDTGridSBSolver<kConfigStSoaNoSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStAosNoSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStAosoaNoSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStSoaSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<kConfigMtOmpSoaSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStAosSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<kConfigMtOmpAosSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<kConfigStAosoaSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTree<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeCusMem<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongest<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongestVec<double>>;
template class UnionUniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<double>>;