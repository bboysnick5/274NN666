//
//  BFEUCPtSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 1/16/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "BFEUCPtSolver.hpp"
#include "Algorithm.hpp"
#include "Utility.hpp"
#include <algorithm>

template <typename FPType, def::ThreadingPolicy Policy>
const SBLoc<FPType>* BFEUCPtSBSolver<FPType, Policy>::
FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return &*Algo<Policy, def::DistType::kEucSq>::LinearNNSearch(this->loc_data_span_.rbegin(), this->loc_data_span_.rend(),
           SBLoc<FPType>::GeoPtTo3dEucPt(geo_search_pt),
           [](const typename SBLoc<FPType>::GeoCart3DPtType& euc_search_pt, const SBLoc<FPType>& loc) {
               return euc_search_pt.template dist<def::DistType::kEucSq>(loc.LocTo3dEucPt());});
}

template <typename FPType, def::ThreadingPolicy Policy>
void BFEUCPtSBSolver<FPType, Policy>::PrintSolverInfo() const {
    std::cout << "This is brute force solver using converted euclidean distance metric.\n";
}

template class BFEUCPtSBSolver<double, def::ThreadingPolicy::kSingle>;
template class BFEUCPtSBSolver<float, def::ThreadingPolicy::kSingle>;
