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
void BKDTSBSolver<KDTType, FPType>::Build(const std::shared_ptr<std::vector<SBLoc<FPType>>> &locData) {
    GenerateKDT(locData);
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void BKDTSBSolver<KDTType, FPType>::PrintSolverInfo() const {
    this->locKdt.printTreeInfo();
}


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void BKDTSBSolver<KDTType, FPType>::GenerateKDT(const std::shared_ptr<std::vector<SBLoc<FPType>>> &locData) {
    std::vector<typename KDT<KDTType, FPType>::node_type> kdtData;
    kdtData.reserve(locData->size());
    std::transform(locData->cbegin(), locData->cend(), std::back_inserter(kdtData),
                   [](const SBLoc<FPType>& l) -> typename KDT<KDTType, FPType>::node_type {
                       return {l.locToCart3DPt(), &l};});
    this->locKdt = KDT<KDTType, FPType>(kdtData.begin(), kdtData.end());
}


template class BKDTSBSolver<KDTree, double>;
template class BKDTSBSolver<KDTree, float>;

template class BKDTSBSolver<KDTreeCusMem, double>;
template class BKDTSBSolver<KDTreeCusMem, float>;

template class BKDTSBSolver<KDTreeExpandLongest, double>;
template class BKDTSBSolver<KDTreeExpandLongest, float>;

template class BKDTSBSolver<KDTreeExpandLongestVec, double>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float>;
