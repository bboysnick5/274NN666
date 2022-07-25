//
//  BFLocalStorageSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 7/9/21.
//  Copyright © 2021 Yunlong Liu. All rights reserved.
//

#include "BFLocalStorageSBSolver.hpp"
#include "Algorithm.hpp"
#include "Utility.hpp"

//#include <oneapi/tbb.h>

template <typename FPType, def::ThreadingPolicy Policy>
void BFLocalStorageSBSolver<FPType, Policy>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    geo_pt_vec_.reserve(loc_data_span.size());
    loc_address_vec_.reserve(loc_data_span.size());
    std::for_each(loc_data_span.rbegin(), loc_data_span.rend(), [&](const auto& loc) {
        geo_pt_vec_.push_back(loc.geo_pt);
        loc_address_vec_.push_back(&loc);
    });
}

template <typename FPType, def::ThreadingPolicy Policy>
const SBLoc<FPType>* BFLocalStorageSBSolver<FPType, Policy>::
FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return loc_address_vec_[Algo<Policy, def::DistType::kHavComp>::LinearNNSearch(geo_pt_vec_.begin(),
                                                                                  geo_pt_vec_.end(), geo_search_pt)
                            - geo_pt_vec_.begin()];
}

template <typename FPType, def::ThreadingPolicy Policy>
void BFLocalStorageSBSolver<FPType, Policy>::PrintSolverInfo() const {
    std::cout << "This is brute force solver with local SBLoc storage using haversine distance metric.\n";
}


template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSingle>;
template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kSingle>;
//template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSimdSoA>;
//template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kSimdSoA>;
template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kMultiOmp>;
template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kMultiOmp>;
//template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kMultiHand>;
//template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kMultiHand>;
//template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSimdMultiOmp>;
//template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kSimdMultiOmp>;
//template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSimdMultiHand>;
//template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kSimdMultiHand>;
