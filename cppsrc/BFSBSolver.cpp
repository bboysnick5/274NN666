//
//  BFSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/11/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//


#include "BFSBSolver.hpp"
#include "Algorithm.hpp"
#include "Utility.hpp"
//#include <oneapi/tbb.h>

template <typename FPType, def::ThreadingPolicy Policy>
void BFSBSolver<FPType, Policy>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    loc_data_span_ = loc_data_span;
}

template <typename FPType, def::ThreadingPolicy Policy>
const SBLoc<FPType>* BFSBSolver<FPType, Policy>::FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return &*Algo<Policy, def::DistType::kHavComp>::LinearNNSearch(loc_data_span_.rbegin(), loc_data_span_.rend(), geo_search_pt, [](const typename SBLoc<FPType>::GeoPtType& geo_search_pt, const SBLoc<FPType>& loc) {return geo_search_pt.template dist<def::DistType::kHavComp>(loc.geo_pt);});
    //return &*oneapi::dpl::min_element(oneapi::dpl::execution::par_unseq, locData->cbegin(), locData->cend(), [&](const SBLoc<FPType> &l1, const SBLoc<FPType> &l2){return l1.havDistComp(geo_search_pt) < l2.havDistComp(geo_search_pt);});
}

template <typename FPType, def::ThreadingPolicy Policy>
void BFSBSolver<FPType, Policy>::PrintSolverInfo() const {
    std::cout << "This is brute force solver using haversine distance metric.\n";
}


template class BFSBSolver<double, def::ThreadingPolicy::kSingle>;
template class BFSBSolver<float, def::ThreadingPolicy::kSingle>;
