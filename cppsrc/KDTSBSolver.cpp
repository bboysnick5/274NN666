//
//  KDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "KDTSBSolver.hpp"
#include <algorithm>


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void KDTSBSolver<KDTType, FPType>::Build(const std::shared_ptr<std::vector<SBLoc<FPType>>> &locData) {
    GenerateKDT(locData);
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void KDTSBSolver<KDTType, FPType>::PrintSolverInfo() const {
    locKdt.printTreeInfo();
}


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
const SBLoc<FPType>* KDTSBSolver<KDTType, FPType>::FindNearestLoc(const PointND<FPType, 2>& geoSearchPt) const {
    return locKdt.kNNValue(SBLoc<FPType>::geoPtToCart3DPt(geoSearchPt), 1);
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void KDTSBSolver<KDTType, FPType>::GenerateKDT(const std::shared_ptr<std::vector<SBLoc<FPType>>> &locData) {
    std::for_each(locData->cbegin(), locData->cend(), [&](const SBLoc<FPType> &loc) mutable {
        this->locKdt.insert(loc.locToCart3DPt(), &loc);});
}


template class KDTSBSolver<KDTree, double>;
template class KDTSBSolver<KDTree, float>;

template class KDTSBSolver<KDTreeCusMem, double>;
template class KDTSBSolver<KDTreeCusMem, float>;

template class KDTSBSolver<KDTreeExpandLongest, double>;
template class KDTSBSolver<KDTreeExpandLongest, float>;

template class KDTSBSolver<KDTreeExpandLongestVec, double>;
template class KDTSBSolver<KDTreeExpandLongestVec, float>;

