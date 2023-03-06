/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-08 20:44:38
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/BFCARTPtSolver.cpp
 */

#include "BFCartPtSolver.hpp"

#include <algorithm>

#include "Algorithm.hpp"
#include "SBLoc.hpp"
#include "Utility.hpp"

template <SolverConfig Config>
void BFCartPtSBSolver<Config>::Build(
    std::span<const SBLoc<FPType>> loc_data_span) {
    this->loc_span_ = loc_data_span;
    loc_cart_pt_vec_.reserve(loc_data_span.size());
    std::transform(loc_data_span.rbegin(), loc_data_span.rend(),
                   std::back_inserter(loc_cart_pt_vec_),
                   [&](const SBLoc<FPType> &loc) {
                       return SBLoc<FPType>::GeoPtToCartPt(loc.geo_pt);
                   });
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *BFCartPtSBSolver<Config>::FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return &this->loc_span_
                [algo::LinearNNSearch<Config.par_policy, def::DistType::kEucSq>(
                     loc_cart_pt_vec_.crbegin(), loc_cart_pt_vec_.crend(),
                     SBLoc<FPType>::GeoPtToCartPt(geo_search_pt)) -
                 loc_cart_pt_vec_.crbegin()];
}

template <SolverConfig Config>
void BFCartPtSBSolver<Config>::PrintSolverInfo() const {
    std::cout << "This is brute force solver using converted euclidean "
                 "distance metric.\n";
}

template class BFCartPtSBSolver<
    SolverConfig<double, def::kParPolicyStSoaNoSimd>{}>;
template class BFCartPtSBSolver<
    SolverConfig<float, def::kParPolicyStSoaNoSimd>{}>;