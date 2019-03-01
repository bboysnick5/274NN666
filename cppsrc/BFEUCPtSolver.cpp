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

template <typename dist_type>
const SBLoc<dist_type>* BFEUCPtSBSolver<dist_type>::findNearest(const Point<dist_type, 2>& geoSearchPt) const {
    const auto testPt = SBLoc<dist_type>::geoPtToCart3DPt(geoSearchPt);
    /*
    return &*std::min_element(this->locData->cbegin(), this->locData->cend(),
                              [&testPt](const SBLoc<dist_type>& l1, const SBLoc<dist_type>& l2) {
                                  return Point<dist_type, 3>::template dist<Point<dist_type, 3>::DistType::EUCSQ>(l1.locToCart3DPt(), testPt)
                                  < Point<dist_type, 3>::template dist<Point<dist_type, 3>::DistType::EUCSQ>(l2.locToCart3DPt(), testPt);}); */
    return &*custom_min_element(this->locData->cbegin(), this->locData->cend(),
                                [&testPt](const SBLoc<dist_type>& l) {
                                    return testPt.template dist<Point<dist_type, 3>::DistType::EUCSQ>(l.locToCart3DPt());},
                                std::less<dist_type>());
}

template <typename dist_type>
void BFEUCPtSBSolver<dist_type>::printSolverInfo() const {
    std::cout << "This is brute force solver using converted euclidean distance metric.\n";
}

template class BFEUCPtSBSolver<double>;
template class BFEUCPtSBSolver<float>;
