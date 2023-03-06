//
//  UniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniCellBKDTGridSBSolver.hpp"

// #include <omp.h>

template <SolverConfig Config>
UniCellBKDTGridSBSolver<Config>::UniCellBKDTGridSBSolver(
    FPType alpc, std::size_t maxCacheCellVecSize)
    : UniLatLngBKDTGridSBSolver<Config>(alpc, maxCacheCellVecSize) {}

template <SolverConfig Config>
void UniCellBKDTGridSBSolver<Config>::FillGridCache() {
    thisRowStartIdx.reserve(this->row_size_);
    this->grid_cache_.reserve(this->loc_kdt_.size() * 1.2 /
                              this->AVE_LOC_PER_CELL);
    std::vector<typename KDTType::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(pt_loc_vec) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    FPType thisLat = -0.5 * def::kMathPi<FPType>,
           thisCtrLat = thisLat + 0.5 * this->lat_inc_;
    for (std::size_t r = 0, idx = 0; r < this->row_size_;
         ++r, thisCtrLat += this->lat_inc_, thisLat += this->lat_inc_) {
        std::size_t thisColSize =
            static_cast<std::size_t>(
                2 * def::kMathPi<FPType> * SBLoc<FPType>::EARTH_RADIUS *
                std::cos(thisLat > 0 ? thisLat - this->lat_inc_ : thisLat) /
                this->side_len_) +
            2;
        FPType thisLngInc =
            2 * def::kMathPi<FPType> / thisColSize +
            2 * def::kMathPi<FPType> / (thisColSize * thisColSize * 65536);
        thisRowStartIdx.emplace_back(idx, 1.0 / thisLngInc);
        FPType lat1 = r * this->lat_inc_ - 0.5 * def::kMathPi<FPType>;
        FPType diagonalDistSq3DCART =
                   SBLoc<FPType>::CART3DDistSqFromLatDeltaLng(
                       lat1, lat1 + this->lat_inc_, thisLngInc),
               thisCtrLng = 0.5 * thisLngInc - def::kMathPi<FPType>;
        for (std::size_t thisEndIdx = idx + thisColSize; idx < thisEndIdx;
             ++idx, thisCtrLng += thisLngInc) {
            UniLatLngBKDTGridSBSolver<Config>::FillCacheCell(
                thisCtrLng, thisCtrLat, diagonalDistSq3DCART, pt_loc_vec);
        }
    }
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *UniCellBKDTGridSBSolver<Config>::FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    const auto &[startIdx, thisLngIncInverse] =
        thisRowStartIdx[(geo_search_pt[0] +
                         0.5 * def::kMathPi<FPType>)*this->lat_inc_inverse_];
    return UniLatLngBKDTGridSBSolver<Config>::ReturnNNLocFromCacheVariant(
        geo_search_pt,
        this->grid_cache_[startIdx +
                          static_cast<std::size_t>(
                              (geo_search_pt[1] +
                               def::kMathPi<FPType>)*thisLngIncInverse)]);
}

template class UniCellBKDTGridSBSolver<kConfigStSoaNoSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpSoaNoSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigStAosNoSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosNoSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigStAosoaNoSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosoaNoSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigStSoaSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpSoaSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigStAosSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigStAosoaSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosoaSimdKDTree<float>>;
template class UniCellBKDTGridSBSolver<kConfigStSoaNoSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<kConfigStAosNoSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<kConfigStAosoaNoSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<kConfigStSoaSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpSoaSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<kConfigStAosSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<kConfigStAosoaSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeCusMem<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongest<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongestVec<float>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<float>>;

template class UniCellBKDTGridSBSolver<kConfigStSoaNoSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpSoaNoSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigStAosNoSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosNoSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigStAosoaNoSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosoaNoSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigStSoaSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpSoaSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigStAosSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigStAosoaSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosoaSimdKDTree<double>>;
template class UniCellBKDTGridSBSolver<kConfigStSoaNoSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<kConfigStAosNoSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<kConfigStSoaSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpSoaSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<kConfigStAosSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<kConfigMtOmpAosSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<kConfigStAosoaSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeCusMem<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongest<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongestVec<double>>;
template class UniCellBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<double>>;