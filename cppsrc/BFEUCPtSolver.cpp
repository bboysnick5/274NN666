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

template <typename FPType>
const SBLoc<FPType>* BFEUCPtSBSolver<FPType>::FindNearestLoc(const PointND<FPType, 2>& geoSearchPt) const {
    const auto testPt = SBLoc<FPType>::geoPtToCart3DPt(geoSearchPt);
    return &*Utility::MinElementGivenDistFunc(loc_data_span_.rbegin(), loc_data_span_.rend(),
                                         [testPt](const SBLoc<FPType>& l) {
                                            return testPt.template dist<PointND<FPType, 3>::DistType::EUCSQ>(l.locToCart3DPt());},
                                         std::less<FPType>());
}

template <typename FPType>
void BFEUCPtSBSolver<FPType>::PrintSolverInfo() const {
    std::cout << "This is brute force solver using converted euclidean distance metric.\n";
}

template class BFEUCPtSBSolver<double>;
template class BFEUCPtSBSolver<float>;
