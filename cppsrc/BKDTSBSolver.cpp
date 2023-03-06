/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-09 10:30:46
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/BKDTSBSolver.cpp
 */

#include "BKDTSBSolver.hpp"

#include <utility>

#include "Definition.hpp"

template <SolverConfig Config>
void BKDTSBSolver<Config>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    GenerateKDT(loc_data_span);
}

template <SolverConfig Config>
void BKDTSBSolver<Config>::PrintSolverInfo() const {
    this->loc_kdt_.PrintTreeInfo();
}

template <SolverConfig Config>
void BKDTSBSolver<Config>::GenerateKDT(
    std::span<const SBLoc<FPType>> loc_data_span) {
    std::vector<typename KDTType::node_type> kdt_data_vec;
    kdt_data_vec.reserve(loc_data_span.size());
    std::transform(loc_data_span.rbegin(), loc_data_span.rend(),
                   std::back_inserter(kdt_data_vec),
                   [](const SBLoc<FPType> &l) -> typename KDTType::node_type {
                       return {l.LocTo3dCartPt(), &l};
                   });
    this->loc_kdt_ = KDTType(kdt_data_vec.begin(), kdt_data_vec.end());
}

template class BKDTSBSolver<kConfigStSoaNoSimdKDTree<float>>;
template class BKDTSBSolver<kConfigMtOmpSoaNoSimdKDTree<float>>;
template class BKDTSBSolver<kConfigStAosNoSimdKDTree<float>>;
template class BKDTSBSolver<kConfigMtOmpAosNoSimdKDTree<float>>;
template class BKDTSBSolver<kConfigStAosoaNoSimdKDTree<float>>;
template class BKDTSBSolver<kConfigMtOmpAosoaNoSimdKDTree<float>>;
template class BKDTSBSolver<kConfigStSoaSimdKDTree<float>>;
template class BKDTSBSolver<kConfigMtOmpSoaSimdKDTree<float>>;
template class BKDTSBSolver<kConfigStAosSimdKDTree<float>>;
template class BKDTSBSolver<kConfigMtOmpAosSimdKDTree<float>>;
template class BKDTSBSolver<kConfigStAosoaSimdKDTree<float>>;
template class BKDTSBSolver<kConfigMtOmpAosoaSimdKDTree<float>>;
template class BKDTSBSolver<kConfigStSoaNoSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigStAosNoSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigMtOmpAosNoSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigStAosoaNoSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigMtOmpAosoaNoSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigStSoaSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigMtOmpSoaSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigStAosSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigMtOmpAosSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigStAosoaSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigMtOmpAosoaSimdKDTreeCusMem<float>>;
template class BKDTSBSolver<kConfigStSoaNoSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigStAosNoSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigStAosoaNoSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigStSoaSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigStAosSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigMtOmpAosSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigStAosoaSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongest<float>>;
template class BKDTSBSolver<kConfigStSoaNoSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigStAosNoSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigStAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigStSoaSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigStAosSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigMtOmpAosSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigStAosoaSimdKDTreeExpandLongestVec<float>>;
template class BKDTSBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<float>>;

template class BKDTSBSolver<kConfigStSoaNoSimdKDTree<double>>;
template class BKDTSBSolver<kConfigMtOmpSoaNoSimdKDTree<double>>;
template class BKDTSBSolver<kConfigStAosNoSimdKDTree<double>>;
template class BKDTSBSolver<kConfigMtOmpAosNoSimdKDTree<double>>;
template class BKDTSBSolver<kConfigStAosoaNoSimdKDTree<double>>;
template class BKDTSBSolver<kConfigMtOmpAosoaNoSimdKDTree<double>>;
template class BKDTSBSolver<kConfigStSoaSimdKDTree<double>>;
template class BKDTSBSolver<kConfigMtOmpSoaSimdKDTree<double>>;
template class BKDTSBSolver<kConfigStAosSimdKDTree<double>>;
template class BKDTSBSolver<kConfigMtOmpAosSimdKDTree<double>>;
template class BKDTSBSolver<kConfigStAosoaSimdKDTree<double>>;
template class BKDTSBSolver<kConfigMtOmpAosoaSimdKDTree<double>>;
template class BKDTSBSolver<kConfigStSoaNoSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigStAosNoSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigMtOmpAosNoSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigStAosoaNoSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigMtOmpAosoaNoSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigStSoaSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigMtOmpSoaSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigStAosSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigMtOmpAosSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigStAosoaSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigMtOmpAosoaSimdKDTreeCusMem<double>>;
template class BKDTSBSolver<kConfigStSoaNoSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigStAosNoSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigStAosoaNoSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigStSoaSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigStAosSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigMtOmpAosSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigStAosoaSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongest<double>>;
template class BKDTSBSolver<kConfigStSoaNoSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<
    kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<kConfigStAosNoSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<
    kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<kConfigStAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<
    kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<kConfigStSoaSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<kConfigStAosSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<kConfigMtOmpAosSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<kConfigStAosoaSimdKDTreeExpandLongestVec<double>>;
template class BKDTSBSolver<
    kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<double>>;