//
//  BKDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "BKDTSBSolver.hpp"

#include <utility>




template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type>
void BKDTSBSolver<KDTType, dist_type>::Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    GenerateKDT(locData);
}

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type>
void BKDTSBSolver<KDTType, dist_type>::PrintSolverInfo() const {
    this->locKdt.printTreeInfo();
}


template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type>
void BKDTSBSolver<KDTType, dist_type>::GenerateKDT(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    std::vector<typename KDT<KDTType, dist_type>::node_type> kdtData;
    kdtData.reserve(locData->size());
    std::transform(locData->cbegin(), locData->cend(), std::back_inserter(kdtData),
                   [](const SBLoc<dist_type>& l) -> typename KDT<KDTType, dist_type>::node_type {
                       return {l.locToCart3DPt(), &l};});
    this->locKdt = KDT<KDTType, dist_type>(kdtData.begin(), kdtData.end());
}


template class BKDTSBSolver<KDTree, double>;
template class BKDTSBSolver<KDTree, float>;

template class BKDTSBSolver<KDTreeCusMem, double>;
template class BKDTSBSolver<KDTreeCusMem, float>;

template class BKDTSBSolver<KDTreeExpandLongest, double>;
template class BKDTSBSolver<KDTreeExpandLongest, float>;

template class BKDTSBSolver<KDTreeExpandLongestVec, double>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float>;
