//
//  KDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "KDTSBSolver.hpp"
#include <algorithm>


template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void KDTSBSolver<KDTType, FPType>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    GenerateKDT(loc_data_span);
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void KDTSBSolver<KDTType, FPType>::PrintSolverInfo() const {
    loc_kdt_.PrintTreeInfo();
}


template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
const SBLoc<FPType>* KDTSBSolver<KDTType, FPType>::FindNearestLoc(PointND<FPType, 2> geo_search_pt) const {
    return loc_kdt_.kNNValue(SBLoc<FPType>::GeoPtTo3dEucPt(geo_search_pt), 1);
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void KDTSBSolver<KDTType, FPType>::GenerateKDT(std::span<const SBLoc<FPType>> loc_data_span) {
    std::for_each(loc_data_span.rbegin(), loc_data_span.rend(), [&](const SBLoc<FPType> &loc) mutable {
        this->loc_kdt_.insert(loc.LocTo3dEucPt(), &loc);});
}


template class KDTSBSolver<KDTree, double>;
template class KDTSBSolver<KDTree, float>;

template class KDTSBSolver<KDTreeCusMem, double>;
template class KDTSBSolver<KDTreeCusMem, float>;

template class KDTSBSolver<KDTreeExpandLongest, double>;
template class KDTSBSolver<KDTreeExpandLongest, float>;

template class KDTSBSolver<KDTreeExpandLongestVec, double>;
template class KDTSBSolver<KDTreeExpandLongestVec, float>;

