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
    std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>> kdtData(sbData->size());
    std::transform(sbData->begin(), sbData->end(), kdtData.begin(),
                   [&](const SBLoc& l)->std::pair<Point<3, DistType::EUC>, const SBLoc*>{
                       return {l.locToCart3DPt(), &l};});
    kdt = KDTree<3, const SBLoc*, DistType::EUC>(kdtData.begin(), kdtData.end());
    std::cout << "Tree height is " << kdt.height() << std::endl
              << "Tree size is " << kdt.size() << std::endl;
}
