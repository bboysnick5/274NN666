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

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
UnionUniLatLngBKDTGridSBSolver(dist_type alpc, std::size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType, dist_type>(), AVE_LOC_PER_CELL(alpc),
MAX_CACHE_CELL_VEC_SIZE(maxCacheCellVecSize) {}

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::PrintSolverInfo() const {
    std::size_t num_single_locs = 0, num_vec_locs = 0, num_tree_locs = 0, num_unique_vec_locs = 0;
    std::size_t total_num_tree_nodes = 0, total_num_vec_locs = 0;
    std::for_each(grid_cache_.cbegin(), grid_cache_.cend(), [&](const BitCell& cell) mutable {
        auto cellPtr = cell.getPtr();
        if (auto cellSize = cell.size(cellPtr);
            cellSize == 1) {
            num_single_locs++;
        } else if (cellSize < MAX_CACHE_CELL_VEC_SIZE) {
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
              << "Ratio of cache locs over actual num locs: " << static_cast<dist_type>(totalCacheLocs)/numGivenLocs << std::endl
              << "Total num of loc cells: " << grid_cache_.size() << std::endl
              << "Ave cached locs per cell: " << static_cast<dist_type>(totalCacheLocs)/grid_cache_.size() << std::endl
              << "Single loc cells: " << num_single_locs << std::endl
              << "Vector loc cells: " << num_vec_locs << std::endl
              << "Unique Vector loc cells: " << num_unique_vec_locs << std::endl
              << "Ave vec loc size: " << static_cast<dist_type>(total_num_vec_locs)/num_vec_locs << std::endl
              << "Kd-tree loc cells: " << num_tree_locs << std::endl
              << "Ave tree size : " << static_cast<dist_type>(total_num_tree_nodes)/num_tree_locs << std::endl
              << "Ave tree height: " << log2(static_cast<dist_type>(total_num_tree_nodes)/num_tree_locs + 1.0) + 1.0 << std::endl;
}

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
LoopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs, def::Policy_Tag<def::ThreadingPolicy::kSingle>) {
    grid_cache_.reserve(rowSize*colSize);
    dist_type lat1 = -0.5*def::kMathPi<dist_type>;
    dist_type thisCtrLat = 0.5 * (latInc - def::kMathPi<dist_type>);
    for (std::size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc, lat1 += latInc) {
        dist_type thisCtrLng = 0.5 * lngInc - def::kMathPi<dist_type>;
        dist_type diagonalDistSq3DEUC = SBLoc<dist_type>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
        for (std::size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            FillCacheCell({thisCtrLat, thisCtrLng}, diagonalDistSq3DEUC, colSize, ptLocPairs);
        }
    }
}

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
LoopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs, def::Policy_Tag<def::ThreadingPolicy::kMultiOmp>) {
    grid_cache_.resize(rowSize*colSize, 0);
    dist_type initCtrLat = 0.5*latInc - 0.5*def::kMathPi<dist_type>;
    dist_type initCtrLng = 0.5*lngInc - def::kMathPi<dist_type>;
    dist_type initLat1 = - 0.5*def::kMathPi<dist_type>;
    
//#pragma ompdeclare reduction (merge : std::vector<BitCell> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) initializer(omp_priv = omp_orig) //only one thread is used
#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
shared(grid_cache_, lngInc, latInc, initCtrLat, initCtrLng, initLat1) \
firstprivate(ptLocPairs) default(none) schedule(dynamic, 1) collapse(2) 
    for (std::size_t r = 0; r < rowSize; ++r) {
        for (std::size_t c = 0; c < colSize; ++c) {
            dist_type lat1 = r*latInc + initLat1;
            dist_type diagonalDistSq3DEUC = SBLoc<dist_type>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
            dist_type thisCtrLat = initCtrLat + r*latInc;
            std::size_t startIdxThisRow = r*colSize;
            FillCacheCell(startIdxThisRow + c, {thisCtrLat, initCtrLng + c*lngInc},
                          diagonalDistSq3DEUC, colSize, ptLocPairs);
        }
    }
}


template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
LoopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs, def::Policy_Tag<def::ThreadingPolicy::kMultiHand>) {
    std::size_t totalCacheCells = rowSize*colSize;
    grid_cache_.resize(totalCacheCells, 0);
    dist_type initCtrLat = 0.5*latInc - 0.5*def::kMathPi<dist_type>;
    dist_type initCtrLng = 0.5*lngInc - def::kMathPi<dist_type>;
    dist_type initLat1 = - 0.5*def::kMathPi<dist_type>;
    
    std::size_t numThreads = std::thread::hardware_concurrency();
    std::size_t chunkSize = totalCacheCells/numThreads;

    
    for (std::size_t r = 0; r < rowSize; ++r) {
        dist_type lat1 = r*latInc + initLat1;
        dist_type diagonalDistSq3DEUC = SBLoc<dist_type>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
        dist_type thisCtrLat = initCtrLat + r*latInc;
        std::size_t startIdxThisRow = r*colSize;
        for (std::size_t c = 0; c < colSize; ++c) {
            FillCacheCell(startIdxThisRow + c, {thisCtrLat, initCtrLng + c*lngInc},
                          diagonalDistSq3DEUC, colSize, ptLocPairs);
        }
    }
}

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::FillGridCache() {
    colSize = rowSize;
    lngInc = 2.0*def::kMathPi<dist_type>/colSize + std::numeric_limits<dist_type>::epsilon();
    lngIncInverse = 1.0/lngInc;
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    LoopBodyThreadingPolicyDispatch(ptLocPairs);
}

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
LoopBodyThreadingPolicyDispatch(std::vector<typename KDT<KDTType, dist_type>::node_type> &ptLocPairs) {
    if constexpr (policy == def::ThreadingPolicy::kSingle) {
        LoopBody(ptLocPairs, def::Policy_Tag<def::ThreadingPolicy::kSingle>{});
    } else if (policy == def::ThreadingPolicy::kMultiOmp) {
        LoopBody(ptLocPairs, def::Policy_Tag<def::ThreadingPolicy::kMultiOmp>{});
    }
}

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::calcSideLenFromAlpc() {
    dist_type surfaceArea = 4.0*def::kMathPi<dist_type>*SBLoc<dist_type>::EARTH_RADIUS*SBLoc<dist_type>::EARTH_RADIUS;
    dist_type numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    numGivenLocs = locData->size();
    BKDTSBSolver<KDTType, dist_type>::GenerateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc<dist_type>::deltaLatOnSameLngFromHavDist(sideLen));
    latIncInverse = 1.0/latInc;
    rowSize = static_cast<std::size_t>(def::kMathPi<dist_type>/latInc) + 1; // equals std::ceil()
    FillGridCache();
    //auto t = std::thread(&KDT<KDTType, dist_type>::clear, this->locKdt);
    //t.detach();
    this->locKdt = {};
}

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
ReturnNNLocFromCacheVariant(const PointND<dist_type, 2>& geoPt, const BitCell& cell) const {
    std::uintptr_t ptr = cell.getPtr();
    if (std::size_t raw_cell_size = cell.RawSizeBits(ptr);
        raw_cell_size == 1) {
        return cell.GetSingleLoc(ptr);
    } else if (raw_cell_size != 0) {
        const auto* loc_pairs = cell.GetLocPairs(ptr);
        const auto pt_3d = SBLoc<dist_type>::geoPtToCart3DPt(geoPt);
        return Utility::MinElementGivenDistFunc(loc_pairs, loc_pairs + raw_cell_size,
                                                [&](const auto& nh) {return pt_3d.template
                                                    dist<PointND<dist_type, 3>::DistType::EUCSQ>(nh.key);},
                                                std::less())->value;
    } else [[unlikely]] {
        return cell.GetCacheTree(ptr)->kNNValue(SBLoc<dist_type>::geoPtToCart3DPt(geoPt), 1);
    }
}


template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
FindNearestLoc(const PointND<dist_type, 2>& geoSearchPt) const {
    return ReturnNNLocFromCacheVariant(geoSearchPt,
           grid_cache_[static_cast<std::size_t>((geoSearchPt[0]+0.5*def::kMathPi<dist_type>)*latIncInverse)*colSize +
                       static_cast<std::size_t>((geoSearchPt[1]+def::kMathPi<dist_type>)*lngIncInverse)]);
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







