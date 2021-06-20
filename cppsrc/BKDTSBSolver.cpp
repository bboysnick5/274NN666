//
//  BKDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "BKDTSBSolver.hpp"

#include <utility>




template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void BKDTSBSolver<KDTType, FPType>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    GenerateKDT(loc_data_span);
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void BKDTSBSolver<KDTType, FPType>::PrintSolverInfo() const {
    this->loc_kdt_.PrintTreeInfo();
}


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void BKDTSBSolver<KDTType, FPType>::GenerateKDT(std::span<const SBLoc<FPType>> loc_data_span) {
    std::vector<typename KDT<KDTType, FPType>::node_type> kdt_data_vec;
    kdt_data_vec.reserve(loc_data_span.size());
    std::transform(loc_data_span.rbegin(), loc_data_span.rend(), std::back_inserter(kdt_data_vec),
                   [](const SBLoc<FPType>& l) -> typename KDT<KDTType, FPType>::node_type {
                       return {l.LocTo3dEucPt(), &l};});
    this->loc_kdt_ = KDT<KDTType, FPType>(kdt_data_vec.begin(), kdt_data_vec.end());
}


template class BKDTSBSolver<KDTree, double>;
template class BKDTSBSolver<KDTree, float>;

template class BKDTSBSolver<KDTreeCusMem, double>;
template class BKDTSBSolver<KDTreeCusMem, float>;

template class BKDTSBSolver<KDTreeExpandLongest, double>;
template class BKDTSBSolver<KDTreeExpandLongest, float>;

template class BKDTSBSolver<KDTreeExpandLongestVec, double>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float>;
