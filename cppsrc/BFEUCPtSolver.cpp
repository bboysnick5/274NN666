//
//  BFEUCPtSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 1/16/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "BFEUCPtSolver.hpp"

const SBLoc* BFEUCPtSBSolver::findNearest(double lat, double lng) const {
    const auto testPt = SBLoc::latLngToCart3DPt(lng, lat);
    return &*std::min_element(locData->cbegin(), locData->cend(),
                              [&](const SBLoc& l1, const SBLoc& l2) {
                                  return Point<3>::dist<Point<3>::DistType::EUCSQ>(l1.locToCart3DPt(), testPt)
                                  < Point<3>::dist<Point<3>::DistType::EUCSQ>(l2.locToCart3DPt(), testPt);});
}


void BFEUCPtSBSolver::printSolverInfo() const {
    std::cout << "This is brute force solver using converted euclidean distance metric.\n";
}
