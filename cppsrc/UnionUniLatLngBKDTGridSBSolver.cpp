//
//  UnionUniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/19/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniLatLngBKDTGridSBSolver.hpp"

#include <omp.h>

#include <memory>
#include <thread>

#include "Algorithm.hpp"
#include "Utility.hpp"

// TO DO:
// 2. Time per search, store each search time in the test loc data;

template <SolverConfig Config>
UnionUniLatLngBKDTGridSBSolver<Config>::UnionUniLatLngBKDTGridSBSolver(
    FPType alpc, std::size_t maxCacheCellVecSize)
    : kAveActualLocsPerCell_(alpc),
      kMaxCacheCellVecSize_(maxCacheCellVecSize) {}

template <SolverConfig Config>
void UnionUniLatLngBKDTGridSBSolver<Config>::PrintSolverInfo() const {
    std::size_t num_single_locs = 0, num_vec_locs = 0, num_tree_locs = 0,
                num_unique_vec_locs = 0;
    std::size_t total_num_tree_nodes = 0, total_num_vec_locs = 0;
    std::for_each(grid_cache_.cbegin(), grid_cache_.cend(),
                  [&](const BitCell &cell) mutable {
                      auto cellPtr = cell.GetRawPtr();
                      if (auto cellSize = cell.size(cellPtr); cellSize == 1) {
                          num_single_locs++;
                      } else if (cellSize < kMaxCacheCellVecSize_) {
                          num_vec_locs++;
                          total_num_vec_locs += cellSize;
                          if (cell.IsUniqueVecLoc(cellPtr))
                              num_unique_vec_locs++;
                      } else {
                          num_tree_locs++;
                          total_num_tree_nodes += cellSize;
                      }
                  });
    std::size_t totalCacheLocs =
        num_single_locs + total_num_vec_locs + total_num_tree_nodes;

    std::cout << "Total cached locs: " << totalCacheLocs << std::endl
              << "Ratio of cache locs over actual num locs: "
              << static_cast<FPType>(totalCacheLocs) / num_actual_locs_
              << std::endl
              << "Total num of loc cells: " << grid_cache_.size() << std::endl
              << "Ave cached locs per cell: "
              << static_cast<FPType>(totalCacheLocs) / grid_cache_.size()
              << std::endl
              << "Single loc cells: " << num_single_locs << std::endl
              << "Vector loc cells: " << num_vec_locs << std::endl
              << "Unique Vector loc cells: " << num_unique_vec_locs << std::endl
              << "Ave vec loc size: "
              << static_cast<FPType>(total_num_vec_locs) / num_vec_locs
              << std::endl
              << "Kd-tree loc cells: " << num_tree_locs << std::endl
              << "Ave tree size : "
              << static_cast<FPType>(total_num_tree_nodes) / num_tree_locs
              << std::endl
              << "Ave tree height: "
              << std::log2(static_cast<FPType>(total_num_tree_nodes) /
                               num_tree_locs +
                           1.0) +
                     1.0
              << std::endl;
}

template <SolverConfig Config>
void UnionUniLatLngBKDTGridSBSolver<Config>::LoopBody(
    def::ThreadingPolicyTag<def::ThreadingPolicy::kSingle>) {
    std::vector<typename KDTType::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
    grid_cache_.reserve(row_size_ * col_size_);
    FPType lat1 = -0.5 * def::kMathPi<FPType>;
    FPType this_ctr_lat = 0.5 * (lat_inc_ - def::kMathPi<FPType>);
    for (std::size_t r = 0; r < row_size_;
         ++r, this_ctr_lat += lat_inc_, lat1 += lat_inc_) {
        FPType thisCtrLng = 0.5 * lng_inc_ - def::kMathPi<FPType>;
        FPType diagonalDistSq3DCART =
            SBLoc<FPType>::CART3DDistSqFromLatDeltaLng(lat1, lat1 + lat_inc_,
                                                       lng_inc_);
        for (std::size_t c = 0; c < col_size_; ++c, thisCtrLng += lng_inc_) {
            FillCacheCell({this_ctr_lat, thisCtrLng}, diagonalDistSq3DCART,
                          col_size_, pt_loc_vec);
        }
    }
}

template <SolverConfig Config>
void UnionUniLatLngBKDTGridSBSolver<Config>::LoopBody(
    def::ThreadingPolicyTag<def::ThreadingPolicy::kMultiOmp>) {
    grid_cache_.resize(row_size_ * col_size_, 0);
    FPType initCtrLat = 0.5 * lat_inc_ - 0.5 * def::kMathPi<FPType>;
    FPType initCtrLng = 0.5 * lng_inc_ - def::kMathPi<FPType>;
    FPType initLat1 = -0.5 * def::kMathPi<FPType>;
    std::vector<typename KDTType::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
// #pragma ompdeclare reduction (merge : std::vector<BitCell> :
//  omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
//  initializer(omp_priv = omp_orig) //only one thread is used
#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
    firstprivate(pt_loc_vec) default(shared) schedule(dynamic, 1) collapse(2)
    for (std::size_t r = 0; r < row_size_; ++r) {
        for (std::size_t c = 0; c < col_size_; ++c) {
            FPType lat1 = r * lat_inc_ + initLat1;
            FPType diagonalDistSq3DCART =
                SBLoc<FPType>::CART3DDistSqFromLatDeltaLng(
                    lat1, lat1 + lat_inc_, lng_inc_);
            FPType thisCtrLat = initCtrLat + r * lat_inc_;
            std::size_t startIdxThisRow = r * col_size_;
            FillCacheCell(startIdxThisRow + c,
                          {thisCtrLat, initCtrLng + c * lng_inc_},
                          diagonalDistSq3DCART, col_size_, pt_loc_vec);
        }
    }
}

template <SolverConfig Config>
void UnionUniLatLngBKDTGridSBSolver<Config>::LoopBody(
    def::ThreadingPolicyTag<def::ThreadingPolicy::kMultiHand>) {
    std::vector<typename KDTType::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
    std::size_t totalCacheCells = row_size_ * col_size_;
    grid_cache_.resize(totalCacheCells, 0);
    FPType initCtrLat = 0.5 * lat_inc_ - 0.5 * def::kMathPi<FPType>;
    FPType initCtrLng = 0.5 * lng_inc_ - def::kMathPi<FPType>;
    FPType initLat1 = -0.5 * def::kMathPi<FPType>;

    std::size_t numThreads = std::thread::hardware_concurrency();
    std::size_t chunkSize = totalCacheCells / numThreads;

    for (std::size_t r = 0; r < row_size_; ++r) {
        FPType lat1 = r * lat_inc_ + initLat1;
        FPType diagonalDistSq3DCART =
            SBLoc<FPType>::CART3DDistSqFromLatDeltaLng(lat1, lat1 + lat_inc_,
                                                       lng_inc_);
        FPType thisCtrLat = initCtrLat + r * lat_inc_;
        std::size_t startIdxThisRow = r * col_size_;
        for (std::size_t c = 0; c < col_size_; ++c) {
            FillCacheCell(startIdxThisRow + c,
                          {thisCtrLat, initCtrLng + c * lng_inc_},
                          diagonalDistSq3DCART, col_size_, pt_loc_vec);
        }
    }
}

template <SolverConfig Config>
void UnionUniLatLngBKDTGridSBSolver<Config>::FillGridCache() {
    col_size_ = row_size_;
    lng_inc_ = 2.0 * def::kMathPi<FPType> / col_size_ +
               std::numeric_limits<FPType>::epsilon();
    lng_inc_inverse_ = 1.0 / lng_inc_;
    LoopBodyThreadingPolicyDispatch();
}

template <SolverConfig Config>
void UnionUniLatLngBKDTGridSBSolver<Config>::LoopBodyThreadingPolicyDispatch() {
    LoopBody(def::ThreadingPolicyTag<Config.par_policy.thread_policy>{});
}

template <SolverConfig Config>
void UnionUniLatLngBKDTGridSBSolver<Config>::calcSideLenFromAlpc() {
    FPType surfaceArea = 4.0 * def::kMathPi<FPType> *
                         SBLoc<FPType>::EARTH_RADIUS *
                         SBLoc<FPType>::EARTH_RADIUS;
    FPType numCells = this->loc_kdt_.size() / kAveActualLocsPerCell_;
    side_len_ = sqrt(surfaceArea / numCells);
}

template <SolverConfig Config>
void UnionUniLatLngBKDTGridSBSolver<Config>::Build(
    std::span<const SBLoc<FPType>> loc_data_span) {
    num_actual_locs_ = loc_data_span.size();
    BKDTSBSolver<Config>::GenerateKDT(loc_data_span);
    calcSideLenFromAlpc();
    lat_inc_ =
        std::fabs(SBLoc<FPType>::deltaLatOnSameLngFromHavDist(side_len_));
    lat_inc_inverse_ = 1.0 / lat_inc_;
    row_size_ = static_cast<std::size_t>(def::kMathPi<FPType> / lat_inc_) +
                1;  // equals std::ceil()
    FillGridCache();
    // auto t = std::thread(&KDTType::clear, this->loc_kdt_);
    // t.detach();
    this->loc_kdt_ = {};
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *UnionUniLatLngBKDTGridSBSolver<Config>::ReturnNNLocFromCacheVariant(
        const typename SBLoc<FPType>::GeoPtType &geo_pt,
        const BitCell &cell) const {
    std::uintptr_t ptr = cell.GetRawPtr();
    if (std::size_t raw_cell_size = cell.RawSizeBits(ptr); raw_cell_size == 1) {
        return cell.GetSingleLoc(ptr);
    } else if (raw_cell_size != 0) {
        return algo::LinearNNSearch<Config.par_policy, def::DistType::kEucSq>(
                   cell.GetLocPairs(ptr), cell.GetLocPairs(ptr) + raw_cell_size,
                   SBLoc<FPType>::GeoPtToCartPt(geo_pt),
                   [](auto loc_pair_it) { return loc_pair_it->key; })
            ->value;
    } else [[unlikely]] {
        return cell.GetCacheTree(ptr)->kNNValue(
            SBLoc<FPType>::GeoPtToCartPt(geo_pt), 1);
    }
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *UnionUniLatLngBKDTGridSBSolver<Config>::FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return ReturnNNLocFromCacheVariant(
        geo_search_pt,
        grid_cache_[static_cast<std::size_t>(
                        (geo_search_pt[0] +
                         0.5 * def::kMathPi<FPType>)*lat_inc_inverse_) *
                        col_size_ +
                    static_cast<std::size_t>(
                        (geo_search_pt[1] +
                         def::kMathPi<FPType>)*lng_inc_inverse_)]);
}

template class UnionUniLatLngBKDTGridSBSolver<kConfigStSoaNoSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigStAosNoSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigStSoaSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigMtOmpSoaSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigStAosSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigMtOmpAosSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigStAosoaSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTree<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeCusMem<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongest<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongestVec<float>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<float>>;

template class UnionUniLatLngBKDTGridSBSolver<kConfigStSoaNoSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigStAosNoSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigStSoaSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigStAosSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<kConfigStAosoaSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTree<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeCusMem<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongest<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStSoaSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpSoaSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigStAosoaSimdKDTreeExpandLongestVec<double>>;
template class UnionUniLatLngBKDTGridSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<double>>;
