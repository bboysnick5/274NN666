/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-09 10:27:46
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/BFSBSolver.cpp
 */

#include "BFSBSolver.hpp"

#include "Algorithm.hpp"
#include "Definition.hpp"
#include "Utility.hpp"
// #include <oneapi/tbb.h>
#include <algorithm>

template <SolverConfig Config>
void BFSBSolver<Config>::Build(
    std::span<const SBLoc<typename decltype(Config)::FPType>> loc_data_span) {
    loc_span_ = loc_data_span;
    std::transform(loc_span_.rbegin(), loc_span_.rend(),
                   std::back_inserter(loc_geo_pt_vec_),
                   [&](const SBLoc<FPType> &loc) { return loc.geo_pt; });
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *BFSBSolver<Config>::FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return &loc_span_
        [algo::LinearNNSearch<Config.par_policy, def::DistType::kHavComp>(
             loc_geo_pt_vec_.rbegin(), loc_geo_pt_vec_.rend(), geo_search_pt) -
         loc_geo_pt_vec_.rbegin()];
    // return &*oneapi::dpl::min_element(oneapi::dpl::execution::par_unseq,
    // locData->cbegin(), locData->cend(), [&](const SBLoc<FPType> &l1, const
    // SBLoc<FPType> &l2){return l1.havDistComp(geo_search_pt) <
    // l2.havDistComp(geo_search_pt);});
}

template <SolverConfig Config>
void BFSBSolver<Config>::PrintSolverInfo() const {
    std::cout
        << "This is brute force solver using haversine distance metric.\n";
}

template class BFSBSolver<kConfigStSoaNoSimdNoTree<double>>;
template class BFSBSolver<kConfigStSoaNoSimdNoTree<float>>;
template class BFSBSolver<kConfigMtOmpSoaNoSimdNoTree<double>>;
template class BFSBSolver<kConfigMtOmpSoaNoSimdNoTree<float>>;
template class BFSBSolver<kConfigStAosNoSimdNoTree<double>>;
template class BFSBSolver<kConfigStAosNoSimdNoTree<float>>;
template class BFSBSolver<kConfigMtOmpAosNoSimdNoTree<double>>;
template class BFSBSolver<kConfigMtOmpAosNoSimdNoTree<float>>;
template class BFSBSolver<kConfigStAosoaNoSimdNoTree<double>>;
template class BFSBSolver<kConfigStAosoaNoSimdNoTree<float>>;
template class BFSBSolver<kConfigMtOmpAosoaNoSimdNoTree<double>>;
template class BFSBSolver<kConfigMtOmpAosoaNoSimdNoTree<float>>;
template class BFSBSolver<kConfigStSoaSimdNoTree<double>>;
template class BFSBSolver<kConfigStSoaSimdNoTree<float>>;
template class BFSBSolver<kConfigMtOmpSoaSimdNoTree<double>>;
template class BFSBSolver<kConfigMtOmpSoaSimdNoTree<float>>;
template class BFSBSolver<kConfigStAosSimdNoTree<double>>;
template class BFSBSolver<kConfigStAosSimdNoTree<float>>;
template class BFSBSolver<kConfigMtOmpAosSimdNoTree<double>>;
template class BFSBSolver<kConfigMtOmpAosSimdNoTree<float>>;
template class BFSBSolver<kConfigStAosoaSimdNoTree<double>>;
template class BFSBSolver<kConfigStAosoaSimdNoTree<float>>;
template class BFSBSolver<kConfigMtOmpAosoaSimdNoTree<double>>;
template class BFSBSolver<kConfigMtOmpAosoaSimdNoTree<float>>;