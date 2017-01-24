//
//  KDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "KDTSBSolver.hpp"
#include <algorithm>

void KDTSBSolver::build() {
    std::for_each(sbData->begin(), sbData->end(), [&](const SBLoc &loc){
        kdt.insert(SBLoc::latLngToCart3DXYZ(loc.lng, loc.lat), &loc);});
    std::cout << "Tree height is " << kdt.height() << std::endl;
}

const SBLoc* KDTSBSolver::findNearest(double lng, double lat) const {
    return kdt.kNNValue(SBLoc::latLngToCart3DXYZ(lng, lat), 1);
}
