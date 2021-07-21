//
//  BFEUCPtSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 1/16/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "BFEUCPtSolver.hpp"
#include "Utility.hpp"
#include <algorithm>

template <typename FPType, def::ThreadingPolicy policy>
const SBLoc<FPType>* BFEUCPtSBSolver<FPType, policy>::
FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    return &*utility::MinElementGivenDistFunc<policy>(this->loc_data_span_.rbegin(), this->loc_data_span_.rend(),
                                              [test_pt = SBLoc<FPType>::GeoPtTo3dEucPt(geo_search_pt)](const SBLoc<FPType>& l) {
                                                  return test_pt.template dist<PointND<FPType, 3>::DistType::kEucSq>(l.LocTo3dEucPt());},
                                              std::less<FPType>());
}

template <typename FPType, def::ThreadingPolicy policy>
void BFEUCPtSBSolver<FPType, policy>::PrintSolverInfo() const {
    std::cout << "This is brute force solver using converted euclidean distance metric.\n";
}

template class BFEUCPtSBSolver<double, def::ThreadingPolicy::kSingle>;
template class BFEUCPtSBSolver<float, def::ThreadingPolicy::kSingle>;
