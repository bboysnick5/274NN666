//
//  BFLocalStorageSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 7/9/21.
//  Copyright © 2021 Yunlong Liu. All rights reserved.
//

#include "BFLocalStorageSBSolver.hpp"
#include "Utility.hpp"
//#include <oneapi/tbb.h>

template <typename FPType, def::ThreadingPolicy policy>
void BFLocalStorageSBSolver<FPType, policy>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    geo_pt_vec_.reserve(loc_data_span.size());
    loc_address_vec_.reserve(loc_data_span.size());
    std::for_each(loc_data_span.rbegin(), loc_data_span.rend(), [&](const auto& loc) {
        geo_pt_vec_.push_back(loc.geo_pt);
        loc_address_vec_.push_back(&loc);
    });
}

template <typename FPType, def::ThreadingPolicy policy>
const SBLoc<FPType>* BFLocalStorageSBSolver<FPType, policy>::
FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return loc_address_vec_[utility::MinElementGivenDistFunc(geo_pt_vec_.begin(), geo_pt_vec_.end(),
                            [&geo_search_pt](const auto& loc_pt) {return geo_search_pt.template dist<SBLoc<FPType>::GeoPtType::DistType::kHavComp>(loc_pt);},
                                std::less<FPType>()) - geo_pt_vec_.begin()];
}

template <typename FPType, def::ThreadingPolicy policy>
void BFLocalStorageSBSolver<FPType, policy>::PrintSolverInfo() const {
    std::cout << "This is brute force solver with local SBLoc storage using haversine distance metric.\n";
}


template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSingle>;
template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kSingle>;
template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSimd>;
template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kSimd>;
template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kMultiOmp>;
template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kMultiOmp>;
template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kMultiHand>;
template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kMultiHand>;
template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSimdMultiOmp>;
template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kSimdMultiOmp>;
template class BFLocalStorageSBSolver<double, def::ThreadingPolicy::kSimdMultiHand>;
template class BFLocalStorageSBSolver<float, def::ThreadingPolicy::kSimdMultiHand>;
