//
//  BKDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include <utility>
#include "BKDTSBSolver.hpp"


void BKDTSBSolver::build() {
    std::vector<std::pair<Point<3>, const SBLoc*>> kdtData;
    std::transform(sbData->begin(), sbData->end(), std::back_inserter(kdtData),
        [&](const SBLoc& loc){ return
        std::make_pair(SBLoc::latLngToCart3DXYZ(loc.lng, loc.lat), &loc);});
    kdt = KDTree<3, const SBLoc*, DistType::EUC>(kdtData.begin(), kdtData.end());
    std::cout << "Tree height is " << kdt.height() << std::endl;
}
