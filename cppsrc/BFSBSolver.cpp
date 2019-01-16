//
//  BFSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/11/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include <algorithm>
#include "BFSBSolver.hpp"


void BFSBSolver::build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    this->locData = locData;
}

const SBLoc* BFSBSolver::findNearest(double lng, double lat) const {
    return &*std::min_element(locData->cbegin(), locData->cend(),
                              [=](const SBLoc& l1, const SBLoc& l2) {
                                  return SBLoc::havDist(l1.lng,l1.lat, lng, lat)
                                  < SBLoc::havDist(l2.lng, l2.lat, lng, lat);});
}

void BFSBSolver::printSolverInfo() const {
    std::cout << "This is brute force solver using haversine distance metric.\n";
}
