//
//  KDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "KDTSBSolver.hpp"
#include <algorithm>

void KDTSBSolver::build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    generateKDT(locData);
    locKdt.printTreeInfo();
}

const SBLoc* KDTSBSolver::findNearest(double lng, double lat) const {
    return locKdt.kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}

void KDTSBSolver::generateKDT(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    std::for_each(locData->begin(), locData->end(), [&](const SBLoc &loc){
        locKdt.insert(loc.locToCart3DPt(), &loc);});
}
