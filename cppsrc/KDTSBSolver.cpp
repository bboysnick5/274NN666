//
//  KDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "KDTSBSolver.hpp"

#include <algorithm>

#include "KDTreeExpandLongest.hpp"
#include "SBSolver.hpp"

template <SolverConfig Config>
void KDTSBSolver<Config>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    GenerateKDT(loc_data_span);
}

template <SolverConfig Config>
void KDTSBSolver<Config>::PrintSolverInfo() const {
    loc_kdt_.PrintTreeInfo();
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *KDTSBSolver<Config>::FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return loc_kdt_.kNNValue(SBLoc<FPType>::GeoPtToCartPt(geo_search_pt), 1);
}

template <SolverConfig Config>
void KDTSBSolver<Config>::GenerateKDT(
    std::span<const SBLoc<FPType>> loc_data_span) {
    std::for_each(loc_data_span.rbegin(), loc_data_span.rend(),
                  [&](const SBLoc<FPType> &loc) mutable {
                      this->loc_kdt_.insert(loc.LocTo3dCartPt(), &loc);
                  });
}

template class KDTSBSolver<kConfigStSoaNoSimdKDTree<float>>;
template class KDTSBSolver<kConfigMtOmpSoaNoSimdKDTree<float>>;
template class KDTSBSolver<kConfigStAosNoSimdKDTree<float>>;
template class KDTSBSolver<kConfigMtOmpAosNoSimdKDTree<float>>;
template class KDTSBSolver<kConfigStAosoaNoSimdKDTree<float>>;
template class KDTSBSolver<kConfigMtOmpAosoaNoSimdKDTree<float>>;
template class KDTSBSolver<kConfigStSoaSimdKDTree<float>>;
template class KDTSBSolver<kConfigMtOmpSoaSimdKDTree<float>>;
template class KDTSBSolver<kConfigStAosSimdKDTree<float>>;
template class KDTSBSolver<kConfigMtOmpAosSimdKDTree<float>>;
template class KDTSBSolver<kConfigStAosoaSimdKDTree<float>>;
template class KDTSBSolver<kConfigMtOmpAosoaSimdKDTree<float>>;
template class KDTSBSolver<kConfigStSoaNoSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigStAosNoSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigMtOmpAosNoSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigStAosoaNoSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigMtOmpAosoaNoSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigStSoaSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigMtOmpSoaSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigStAosSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigMtOmpAosSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigStAosoaSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigMtOmpAosoaSimdKDTreeCusMem<float>>;
template class KDTSBSolver<kConfigStSoaNoSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigStAosNoSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigStAosoaNoSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigStSoaSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigStAosSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigMtOmpAosSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigStAosoaSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongest<float>>;
template class KDTSBSolver<kConfigStSoaNoSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigStAosNoSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigStAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigStSoaSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigStAosSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigMtOmpAosSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigStAosoaSimdKDTreeExpandLongestVec<float>>;
template class KDTSBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<float>>;

template class KDTSBSolver<kConfigStSoaNoSimdKDTree<double>>;
template class KDTSBSolver<kConfigMtOmpSoaNoSimdKDTree<double>>;
template class KDTSBSolver<kConfigStAosNoSimdKDTree<double>>;
template class KDTSBSolver<kConfigMtOmpAosNoSimdKDTree<double>>;
template class KDTSBSolver<kConfigStAosoaNoSimdKDTree<double>>;
template class KDTSBSolver<kConfigMtOmpAosoaNoSimdKDTree<double>>;
template class KDTSBSolver<kConfigStSoaSimdKDTree<double>>;
template class KDTSBSolver<kConfigMtOmpSoaSimdKDTree<double>>;
template class KDTSBSolver<kConfigStAosSimdKDTree<double>>;
template class KDTSBSolver<kConfigMtOmpAosSimdKDTree<double>>;
template class KDTSBSolver<kConfigStAosoaSimdKDTree<double>>;
template class KDTSBSolver<kConfigMtOmpAosoaSimdKDTree<double>>;
template class KDTSBSolver<kConfigStSoaNoSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigStAosNoSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigMtOmpAosNoSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigStAosoaNoSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigMtOmpAosoaNoSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigStSoaSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigMtOmpSoaSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigStAosSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigMtOmpAosSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigStAosoaSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigMtOmpAosoaSimdKDTreeCusMem<double>>;
template class KDTSBSolver<kConfigStSoaNoSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigStAosNoSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigStAosoaNoSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigStSoaSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigStAosSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigMtOmpAosSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigStAosoaSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongest<double>>;
template class KDTSBSolver<kConfigStSoaNoSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigStAosNoSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigStAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigStSoaSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigStAosSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigMtOmpAosSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigStAosoaSimdKDTreeExpandLongestVec<double>>;
template class KDTSBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<double>>;