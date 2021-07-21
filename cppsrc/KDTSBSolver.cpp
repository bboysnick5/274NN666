//
//  KDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "KDTSBSolver.hpp"
#include <algorithm>


template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
void KDTSBSolver<KDTType, FPType, policy>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    GenerateKDT(loc_data_span);
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
void KDTSBSolver<KDTType, FPType, policy>::PrintSolverInfo() const {
    loc_kdt_.PrintTreeInfo();
}


template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
const SBLoc<FPType>* KDTSBSolver<KDTType, FPType, policy>::FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return loc_kdt_.kNNValue(SBLoc<FPType>::GeoPtTo3dEucPt(geo_search_pt), 1);
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
void KDTSBSolver<KDTType, FPType, policy>::GenerateKDT(std::span<const SBLoc<FPType>> loc_data_span) {
    std::for_each(loc_data_span.rbegin(), loc_data_span.rend(), [&](const SBLoc<FPType> &loc) mutable {
        this->loc_kdt_.insert(loc.LocTo3dEucPt(), &loc);});
}


template class KDTSBSolver<KDTree, double, def::ThreadingPolicy::kSingle>;
template class KDTSBSolver<KDTree, float, def::ThreadingPolicy::kSingle>;
template class KDTSBSolver<KDTree, double, def::ThreadingPolicy::kSimd>;
template class KDTSBSolver<KDTree, float, def::ThreadingPolicy::kSimd>;
template class KDTSBSolver<KDTree, double, def::ThreadingPolicy::kMultiOmp>;
template class KDTSBSolver<KDTree, float, def::ThreadingPolicy::kMultiOmp>;
template class KDTSBSolver<KDTree, double, def::ThreadingPolicy::kMultiHand>;
template class KDTSBSolver<KDTree, float, def::ThreadingPolicy::kMultiHand>;
template class KDTSBSolver<KDTree, double, def::ThreadingPolicy::kSimdMultiOmp>;
template class KDTSBSolver<KDTree, float, def::ThreadingPolicy::kSimdMultiOmp>;
template class KDTSBSolver<KDTree, double, def::ThreadingPolicy::kSimdMultiHand>;
template class KDTSBSolver<KDTree, float, def::ThreadingPolicy::kSimdMultiHand>;

template class KDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSingle>;
template class KDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSingle>;
template class KDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSimd>;
template class KDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSimd>;
template class KDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiOmp>;
template class KDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiOmp>;
template class KDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiHand>;
template class KDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiHand>;
template class KDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSimdMultiOmp>;
template class KDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSimdMultiOmp>;
template class KDTSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSimdMultiHand>;
template class KDTSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSimdMultiHand>;

template class KDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSingle>;
template class KDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSingle>;
template class KDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSimd>;
template class KDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSimd>;
template class KDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiOmp>;
template class KDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiOmp>;
template class KDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiHand>;
template class KDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiHand>;
template class KDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSimdMultiOmp>;
template class KDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSimdMultiOmp>;
template class KDTSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSimdMultiHand>;
template class KDTSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSimdMultiHand>;

template class KDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSingle>;
template class KDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSingle>;
template class KDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSimd>;
template class KDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSimd>;
template class KDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiOmp>;
template class KDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiOmp>;
template class KDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiHand>;
template class KDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiHand>;
template class KDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSimdMultiOmp>;
template class KDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSimdMultiOmp>;
template class KDTSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSimdMultiHand>;
template class KDTSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSimdMultiHand>;

