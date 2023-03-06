//
//  UniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniLatLngBKDTGridSBSolver.hpp"

#include <vector>

#include "Algorithm.hpp"
#include "Utility.hpp"
// #include <omp.h>

template <SolverConfig Config>
UniLatLngBKDTGridSBSolver<Config>::UniLatLngBKDTGridSBSolver(
    FPType alpc, std::size_t maxCacheCellVecSize)
    : AVE_LOC_PER_CELL(alpc), kMaxCacheCellVecSize_(maxCacheCellVecSize) {}

template <SolverConfig Config>
void UniLatLngBKDTGridSBSolver<Config>::PrintSolverInfo() const {
    std::cout << "Total cache locs: " << totalNodeSize
              << "\nAve tree size: " << totalNodeSize / grid_cache_.size()
              << "\nAve tree height: "
              << static_cast<std::size_t>(
                     log2(totalNodeSize / grid_cache_.size() + 1)) +
                     1
              << "\nRatio of cache locs over actual num locs: "
              << totalNodeSize / totalLocSize
              << "\nTotal num of loc cells: " << grid_cache_.size()
              << "\nSingle loc cells: " << singleLocs
              << "\nVector loc cells: " << vecLocs << "\nKd-tree loc cells: "
              << grid_cache_.size() - singleLocs - vecLocs << "\n";
}

template <SolverConfig Config>
void UniLatLngBKDTGridSBSolver<Config>::FillCacheCell(
    FPType thisCtrLng, FPType thisCtrLat, FPType diagonalDistSq3DCART,
    std::vector<typename KDTType::node_type> &pt_loc_vec) {
    this->loc_kdt_.NNsWithFence(
        SBLoc<FPType>::GeoPtToCartPt({thisCtrLat, thisCtrLng}),
        diagonalDistSq3DCART, std::back_inserter(pt_loc_vec));
    std::size_t locsSize = pt_loc_vec.size();
    this->totalNodeSize += locsSize;
    if (locsSize == 1) {
        this->grid_cache_.emplace_back(pt_loc_vec[0].value);
        this->singleLocs++;
    } else if (locsSize < kMaxCacheCellVecSize_) {
        this->grid_cache_.emplace_back(pt_loc_vec);
        this->vecLocs++;
    } else {
        this->grid_cache_.emplace_back(std::in_place_type<KDTType>,
                                       pt_loc_vec.begin(), pt_loc_vec.end());
    }
    pt_loc_vec.clear();
}

template <SolverConfig Config>
void UniLatLngBKDTGridSBSolver<Config>::FillGridCache() {
    col_size_ = row_size_;
    lng_inc_ = 2.0 * def::kMathPi<FPType> / col_size_ +
               2.0 * def::kMathPi<FPType> / (col_size_ * 65536);
    grid_cache_.reserve(this->loc_kdt_.size() * 1.2 / AVE_LOC_PER_CELL);
    std::vector<typename KDTType::node_type> pt_loc_vec;
    pt_loc_vec.reserve(kMaxCacheCellVecSize_);

    FPType thisCtrLat = 0.5 * (lat_inc_ - def::kMathPi<FPType>);
    for (std::size_t r = 0; r < row_size_; ++r, thisCtrLat += lat_inc_) {
        FPType thisCtrLng = 0.5 * lng_inc_ - def::kMathPi<FPType>;
        FPType lat1 = r * this->lat_inc_ - 0.5 * def::kMathPi<FPType>;
        FPType diagonalDistSq3DCART =
            SBLoc<FPType>::CART3DDistSqFromLatDeltaLng(lat1, lat1 + lat_inc_,
                                                       lng_inc_);
        for (std::size_t c = 0; c < col_size_; ++c, thisCtrLng += lng_inc_) {
            FillCacheCell(thisCtrLng, thisCtrLat, diagonalDistSq3DCART,
                          pt_loc_vec);
        }
    }
}

template <SolverConfig Config>
void UniLatLngBKDTGridSBSolver<Config>::calcSideLenFromAlpc() {
    FPType surfaceArea = 4 * def::kMathPi<FPType> *
                         SBLoc<FPType>::EARTH_RADIUS *
                         SBLoc<FPType>::EARTH_RADIUS;
    FPType numCells = this->loc_kdt_.size() / AVE_LOC_PER_CELL;
    side_len_ = sqrt(surfaceArea / numCells);
}

template <SolverConfig Config>
void UniLatLngBKDTGridSBSolver<Config>::Build(
    std::span<const SBLoc<FPType>> loc_data_span) {
    totalLocSize = loc_data_span.size();
    BKDTSBSolver<Config>::GenerateKDT(loc_data_span);
    calcSideLenFromAlpc();
    lat_inc_ =
        std::fabs(SBLoc<FPType>::deltaLatOnSameLngFromHavDist(side_len_));
    lat_inc_inverse_ = 1.0 / lat_inc_;
    row_size_ = std::ceil(def::kMathPi<FPType> / lat_inc_);
    FillGridCache();
    this->loc_kdt_.Clear();
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *UniLatLngBKDTGridSBSolver<Config>::ReturnNNLocFromCacheVariant(
        const typename SBLoc<FPType>::GeoPtType &geo_search_pt,
        const std::variant<std::vector<typename KDTType::node_type>,
                           const SBLoc<FPType> *, KDTType> &v) const {
    switch (v.index()) {
        case 0: {
            const auto cart_search_pt =
                SBLoc<FPType>::GeoPtToCartPt(geo_search_pt);
            const auto &vec = std::get<0>(v);
            return algo::LinearNNSearch<Config.par_policy,
                                        def::DistType::kEucSq>(
                       vec.cbegin(), vec.cend(), cart_search_pt,
                       [](auto node_vec_it) { return node_vec_it->key; })
                ->value;
        }
        case 1:
            return std::get<1>(v);
        default:
            return std::get<2>(v).kNNValue(
                SBLoc<FPType>::GeoPtToCartPt(geo_search_pt), 1);
    }
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *UniLatLngBKDTGridSBSolver<Config>::FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return ReturnNNLocFromCacheVariant(
        geo_search_pt,
        grid_cache_[static_cast<std::size_t>(
                        (geo_search_pt[0] +
                         0.5 * def::kMathPi<FPType>)*lat_inc_inverse_) *
                        col_size_ +
                    static_cast<std::size_t>(
                        (geo_search_pt[1] + def::kMathPi<FPType>) / lng_inc_)]);
}

template class UniLatLngBKDTGridSBSolver<kConfigStSoaNoSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpSoaNoSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosNoSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpAosNoSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosoaNoSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpAosoaNoSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStSoaSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpSoaSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpAosSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosoaSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpAosoaSimdKDTree<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStSoaNoSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosNoSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStSoaSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosoaSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeCusMem<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongest<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongestVec<float>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<float>>;

template class UniLatLngBKDTGridSBSolver<kConfigStSoaNoSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpSoaNoSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosNoSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpAosNoSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosoaNoSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpAosoaNoSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigStSoaSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpSoaSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpAosSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosoaSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigMtOmpAosoaSimdKDTree<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigStSoaSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<kConfigStAosSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeCusMem<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongest<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongestVec<double>>;
template class UniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<double>>;
