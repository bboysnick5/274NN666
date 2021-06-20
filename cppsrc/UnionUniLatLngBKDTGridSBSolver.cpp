//
//  UnionUniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/19/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniLatLngBKDTGridSBSolver.hpp"
#include "Utility.hpp"
#include <omp.h>
#include <memory>
#include <thread>


// TO DO:
// 2. Time per search, store each search time in the test loc data;

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
UnionUniLatLngBKDTGridSBSolver(FPType alpc, std::size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType, FPType>(), kAveActualLocsPerCell_(alpc),
kMaxCacheCellVecSize_(maxCacheCellVecSize) {}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::PrintSolverInfo() const {
    std::size_t num_single_locs = 0, num_vec_locs = 0, num_tree_locs = 0, num_unique_vec_locs = 0;
    std::size_t total_num_tree_nodes = 0, total_num_vec_locs = 0;
    std::for_each(grid_cache_.cbegin(), grid_cache_.cend(), [&](const BitCell& cell) mutable {
        auto cellPtr = cell.GetPtr();
        if (auto cellSize = cell.size(cellPtr);
            cellSize == 1) {
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
    std::size_t totalCacheLocs = num_single_locs + total_num_vec_locs + total_num_tree_nodes;
    
    std::cout << "Total cached locs: " << totalCacheLocs << std::endl
              << "Ratio of cache locs over actual num locs: " << static_cast<FPType>(totalCacheLocs)/num_actual_locs_ << std::endl
              << "Total num of loc cells: " << grid_cache_.size() << std::endl
              << "Ave cached locs per cell: " << static_cast<FPType>(totalCacheLocs)/grid_cache_.size() << std::endl
              << "Single loc cells: " << num_single_locs << std::endl
              << "Vector loc cells: " << num_vec_locs << std::endl
              << "Unique Vector loc cells: " << num_unique_vec_locs << std::endl
              << "Ave vec loc size: " << static_cast<FPType>(total_num_vec_locs)/num_vec_locs << std::endl
              << "Kd-tree loc cells: " << num_tree_locs << std::endl
              << "Ave tree size : " << static_cast<FPType>(total_num_tree_nodes)/num_tree_locs << std::endl
              << "Ave tree height: " << log2(static_cast<FPType>(total_num_tree_nodes)/num_tree_locs + 1.0) + 1.0 << std::endl;
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
LoopBody(def::Policy_Tag<def::ThreadingPolicy::kSingle>) {
    std::vector<typename KDT<KDTType, FPType>::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
    grid_cache_.reserve(row_size_*col_size_);
    FPType lat1 = -0.5*def::kMathPi<FPType>;
    FPType this_ctr_lat = 0.5 * (lat_inc_ - def::kMathPi<FPType>);
    for (std::size_t r = 0; r < row_size_; ++r, this_ctr_lat += lat_inc_, lat1 += lat_inc_) {
        FPType thisCtrLng = 0.5 * lng_inc_ - def::kMathPi<FPType>;
        FPType diagonalDistSq3DEUC = SBLoc<FPType>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + lat_inc_, lng_inc_);
        for (std::size_t c = 0; c < col_size_; ++c, thisCtrLng += lng_inc_) {
            FillCacheCell({this_ctr_lat, thisCtrLng}, diagonalDistSq3DEUC, col_size_, pt_loc_vec);
        }
    }
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
LoopBody(def::Policy_Tag<def::ThreadingPolicy::kMultiOmp>) {
    grid_cache_.resize(row_size_*col_size_, 0);
    FPType initCtrLat = 0.5*lat_inc_ - 0.5*def::kMathPi<FPType>;
    FPType initCtrLng = 0.5*lng_inc_ - def::kMathPi<FPType>;
    FPType initLat1 = - 0.5*def::kMathPi<FPType>;
    std::vector<typename KDT<KDTType, FPType>::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
//#pragma ompdeclare reduction (merge : std::vector<BitCell> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) initializer(omp_priv = omp_orig) //only one thread is used
#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
firstprivate(pt_loc_vec) default(shared) schedule(dynamic, 1) collapse(2) 
    for (std::size_t r = 0; r < row_size_; ++r) {
        for (std::size_t c = 0; c < col_size_; ++c) {
            FPType lat1 = r*lat_inc_ + initLat1;
            FPType diagonalDistSq3DEUC = SBLoc<FPType>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + lat_inc_, lng_inc_);
            FPType thisCtrLat = initCtrLat + r*lat_inc_;
            std::size_t startIdxThisRow = r*col_size_;
            FillCacheCell(startIdxThisRow + c, {thisCtrLat, initCtrLng + c*lng_inc_},
                          diagonalDistSq3DEUC, col_size_, pt_loc_vec);
        }
    }
}


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
LoopBody(def::Policy_Tag<def::ThreadingPolicy::kMultiHand>) {
    std::vector<typename KDT<KDTType, FPType>::node_type> pt_loc_vec;
    pt_loc_vec.reserve(this->kMaxCacheCellVecSize_);
    std::size_t totalCacheCells = row_size_*col_size_;
    grid_cache_.resize(totalCacheCells, 0);
    FPType initCtrLat = 0.5*lat_inc_ - 0.5*def::kMathPi<FPType>;
    FPType initCtrLng = 0.5*lng_inc_ - def::kMathPi<FPType>;
    FPType initLat1 = - 0.5*def::kMathPi<FPType>;
    
    std::size_t numThreads = std::thread::hardware_concurrency();
    std::size_t chunkSize = totalCacheCells/numThreads;

    
    for (std::size_t r = 0; r < row_size_; ++r) {
        FPType lat1 = r*lat_inc_ + initLat1;
        FPType diagonalDistSq3DEUC = SBLoc<FPType>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + lat_inc_, lng_inc_);
        FPType thisCtrLat = initCtrLat + r*lat_inc_;
        std::size_t startIdxThisRow = r*col_size_;
        for (std::size_t c = 0; c < col_size_; ++c) {
            FillCacheCell(startIdxThisRow + c, {thisCtrLat, initCtrLng + c*lng_inc_},
                          diagonalDistSq3DEUC, col_size_, pt_loc_vec);
        }
    }
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::FillGridCache() {
    col_size_ = row_size_;
    lng_inc_ = 2.0*def::kMathPi<FPType>/col_size_ + std::numeric_limits<FPType>::epsilon();
    lng_inc_inverse_ = 1.0/lng_inc_;
    LoopBodyThreadingPolicyDispatch();
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
LoopBodyThreadingPolicyDispatch() {
    if constexpr (policy == def::ThreadingPolicy::kSingle) {
        LoopBody(def::Policy_Tag<def::ThreadingPolicy::kSingle>{});
    } else if (policy == def::ThreadingPolicy::kMultiOmp) {
        LoopBody(def::Policy_Tag<def::ThreadingPolicy::kMultiOmp>{});
    }
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::calcSideLenFromAlpc() {
    FPType surfaceArea = 4.0*def::kMathPi<FPType>*SBLoc<FPType>::EARTH_RADIUS*SBLoc<FPType>::EARTH_RADIUS;
    FPType numCells = this->loc_kdt_.size()/kAveActualLocsPerCell_;
    side_len_ = sqrt(surfaceArea/numCells);
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
Build(std::span<const SBLoc<FPType>> loc_data_span) {
    num_actual_locs_ = loc_data_span.size();
    BKDTSBSolver<KDTType, FPType>::GenerateKDT(loc_data_span);
    calcSideLenFromAlpc();
    lat_inc_ = std::fabs(SBLoc<FPType>::deltaLatOnSameLngFromHavDist(side_len_));
    lat_inc_inverse_ = 1.0/lat_inc_;
    row_size_ = static_cast<std::size_t>(def::kMathPi<FPType>/lat_inc_) + 1; // equals std::ceil()
    FillGridCache();
    //auto t = std::thread(&KDT<KDTType, FPType>::clear, this->loc_kdt_);
    //t.detach();
    this->loc_kdt_ = {};
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
const SBLoc<FPType>* UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
ReturnNNLocFromCacheVariant(const PointND<FPType, 2>& geoPt, const BitCell& cell) const {
    std::uintptr_t ptr = cell.GetPtr();
    if (std::size_t raw_cell_size = cell.RawSizeBits(ptr);
        raw_cell_size == 1) {
        return cell.GetSingleLoc(ptr);
    } else if (raw_cell_size != 0) {
        const auto* loc_pairs = cell.GetLocPairs(ptr);
        const auto pt_3d = SBLoc<FPType>::GeoPtTo3dEucPt(geoPt);
        return Utility::MinElementGivenDistFunc(loc_pairs, loc_pairs + raw_cell_size,
                                                [&](const auto& nh) {return pt_3d.template
                                                    dist<PointND<FPType, 3>::DistType::EUCSQ>(nh.key);},
                                                std::less())->value;
    } else [[unlikely]] {
        return cell.GetCacheTree(ptr)->kNNValue(SBLoc<FPType>::GeoPtTo3dEucPt(geoPt), 1);
    }
}


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
const SBLoc<FPType>* UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
FindNearestLoc(const PointND<FPType, 2>& geo_search_pt) const {
    return ReturnNNLocFromCacheVariant(geo_search_pt,
           grid_cache_[static_cast<std::size_t>((geo_search_pt[0]+0.5*def::kMathPi<FPType>)*lat_inc_inverse_)*col_size_ +
                       static_cast<std::size_t>((geo_search_pt[1]+def::kMathPi<FPType>)*lng_inc_inverse_)]);
}



template class UnionUniLatLngBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kSingle>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kSingle>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kMultiHand>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSingle>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSingle>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiHand>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSingle>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSingle>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiHand>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSingle>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSingle>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiHand>;







