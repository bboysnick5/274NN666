//
//  BKDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include <utility>
#include "BKDTSBSolver.hpp"


void BKDTSBSolver::build(const std::vector<SBLoc> &sbData) {
    std::vector<std::pair<Point<3>, SBLoc>> kdtData;
    std::transform(sbData.begin(), sbData.end(), std::back_inserter(kdtData),
        [&](const SBLoc& loc){return std::make_pair(transLatLngToXYZPt(loc.lng,
        loc.lat), loc);});
    kdt = KDTree<3, SBLoc>(kdtData.begin(), kdtData.end());
    std::cout << "Tree height is " << kdt.height() << std::endl;
}

SBLoc BKDTSBSolver::findNearest(double lng, double lat) {
    return kdt.kNNValue(transLatLngToXYZPt(lng, lat), 1);
}
