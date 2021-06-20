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
const SBLoc<FPType>* BFEUCPtSBSolver<FPType>::FindNearestLoc(const PointND<FPType, 2>& geo_search_pt) const {
    const auto test_pt = SBLoc<FPType>::GeoPtTo3dEucPt(geo_search_pt);
    return &*Utility::MinElementGivenDistFunc(this->loc_data_span_.rbegin(), this->loc_data_span_.rend(),
                                         [test_pt](const SBLoc<FPType>& l) {
                                            return test_pt.template dist<PointND<FPType, 3>::DistType::EUCSQ>(l.LocTo3dEucPt());},
                                         std::less<FPType>());
}

template <typename FPType>
void BFEUCPtSBSolver<FPType>::PrintSolverInfo() const {
    std::cout << "This is brute force solver using converted euclidean distance metric.\n";
}

template class BFEUCPtSBSolver<double>;
template class BFEUCPtSBSolver<float>;
