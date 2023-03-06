/*
 * @Author: Nick Liu
 * @Date: 2021-07-09 13:33:01
 * @LastEditTime: 2022-08-08 20:45:12
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/BFLocalStorageSBSolver.cpp
 */

#include "BFLocalStorageSBSolver.hpp"

#include <algorithm>

#include "Algorithm.hpp"
#include "Utility.hpp"

// #include <oneapi/tbb.h>

template <SolverConfig Config>
void BFLocalStorageSBSolver<Config>::Build(
    std::span<const SBLoc<FPType>> loc_data_span) {
    geo_pt_vec_.reserve(loc_data_span.size());
    loc_address_vec_.reserve(loc_data_span.size());
    std::for_each(loc_data_span.rbegin(), loc_data_span.rend(),
                  [&](const auto &loc) {
                      geo_pt_vec_.push_back(loc.geo_pt);
                      loc_address_vec_.push_back(&loc);
                  });
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *BFLocalStorageSBSolver<Config>::FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return loc_address_vec_
        [algo::LinearNNSearch<Config.par_policy, def::DistType::kHavComp>(
             geo_pt_vec_.begin(), geo_pt_vec_.end(), geo_search_pt) -
         geo_pt_vec_.begin()];
}

template <SolverConfig Config>
void BFLocalStorageSBSolver<Config>::PrintSolverInfo() const {
    std::cout << "This is brute force solver with local SBLoc storage using "
                 "haversine distance metric.\n";
}

template class BFLocalStorageSBSolver<
    SolverConfig<double, def::kParPolicyStSoaNoSimd>{}>;
template class BFLocalStorageSBSolver<
    SolverConfig<float, def::kParPolicyStSoaNoSimd>{}>;
// template class
// BFLocalStorageSBSolver<double,
// def::ThreadingPolicy::kSimdSoA>; template class BFLocalStorageSBSolver<float,
// def::ThreadingPolicy::kSimdSoA>;
template class BFLocalStorageSBSolver<
    SolverConfig<double, def::kParPolicyMtOmpSoaNoSimd>{}>;
template class BFLocalStorageSBSolver<
    SolverConfig<float, def::kParPolicyMtOmpSoaNoSimd>{}>;

// template class
// BFLocalStorageSBSolver<double,
// def::ThreadingPolicy::kMultiHand>; template class
// BFLocalStorageSBSolver<float, def::ThreadingPolicy::kMultiHand>; template
// class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSimdMultiOmp>;
// template class BFLocalStorageSBSolver<float,
// def::ThreadingPolicy::kSimdMultiOmp>; template class
// BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSimdMultiHand>;
// template class BFLocalStorageSBSolver<float,
// def::ThreadingPolicy::kSimdMultiHand>;
