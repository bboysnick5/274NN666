//
//  BKDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "BKDTSBSolver.hpp"

#include <utility>




template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
void BKDTSBSolver<KDTType, FPType, Policy>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    GenerateKDT(loc_data_span);
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
void BKDTSBSolver<KDTType, FPType, Policy>::PrintSolverInfo() const {
    this->loc_kdt_.PrintTreeInfo();
}


template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
void BKDTSBSolver<KDTType, FPType, Policy>::GenerateKDT(std::span<const SBLoc<FPType>> loc_data_span) {
    std::vector<typename KDT<KDTType, FPType>::node_type> kdt_data_vec;
    kdt_data_vec.reserve(loc_data_span.size());
    std::transform(loc_data_span.rbegin(), loc_data_span.rend(), std::back_inserter(kdt_data_vec),
                   [](const SBLoc<FPType>& l) -> typename KDT<KDTType, FPType>::node_type {
                       return {l.LocTo3dEucPt(), &l};});
    this->loc_kdt_ = KDT<KDTType, FPType>(kdt_data_vec.begin(), kdt_data_vec.end());
}


template class BKDTSBSolver<KDTree, double, def::ThreadingPolicy::kSingle>;
template class BKDTSBSolver<KDTree, float, def::ThreadingPolicy::kSingle>;
//template class BKDTSBSolver<KDTree, double, def::ThreadingPolicy::kSimd>;
//template class BKDTSBSolver<KDTree, float, def::ThreadingPolicy::kSimd>;
template class BKDTSBSolver<KDTree, double, def::ThreadingPolicy::kMultiOmp>;
template class BKDTSBSolver<KDTree, float, def::ThreadingPolicy::kMultiOmp>;
template class BKDTSBSolver<KDTree, double, def::ThreadingPolicy::kMultiHand>;
template class BKDTSBSolver<KDTree, float, def::ThreadingPolicy::kMultiHand>;
//template class BKDTSBSolver<KDTree, double, def::ThreadingPolicy::kSimdMultiOmp>;
//template class BKDTSBSolver<KDTree, float, def::ThreadingPolicy::kSimdMultiOmp>;
//template class BKDTSBSolver<KDTree, double, def::ThreadingPolicy::kSimdMultiHand>;
//template class BKDTSBSolver<KDTree, float, def::ThreadingPolicy::kSimdMultiHand>;

template class BKDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSingle>;
template class BKDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSingle>;
//template class BKDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSimd>;
//template class BKDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSimd>;
template class BKDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiOmp>;
template class BKDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiOmp>;
template class BKDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiHand>;
template class BKDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiHand>;
//template class BKDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSimdMultiOmp>;
//template class BKDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSimdMultiOmp>;
//template class BKDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSimdMultiHand>;
//template class BKDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSimdMultiHand>;

template class BKDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSingle>;
template class BKDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSingle>;
//template class BKDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSimd>;
//template class BKDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSimd>;
template class BKDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiOmp>;
template class BKDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiOmp>;
template class BKDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiHand>;
template class BKDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiHand>;
//template class BKDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSimdMultiOmp>;
//template class BKDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSimdMultiOmp>;
//template class BKDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSimdMultiHand>;
//template class BKDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSimdMultiHand>;

template class BKDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSingle>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSingle>;
template class BKDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiOmp>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiOmp>;
template class BKDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiHand>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiHand>;
/*
template class BKDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSimd>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSimd>;
template class BKDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSimdMultiOmp>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSimdMultiOmp>;
template class BKDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSimdMultiHand>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSimdMultiHand>;
*/
