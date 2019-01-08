//
//  BKDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include <utility>
#include "BKDTSBSolver.hpp"


void BKDTSBSolver::build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    generateKDT(locData);
    locKdt.printTreeInfo();
}

void BKDTSBSolver::generateKDT(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    std::vector<std::pair<Point<3>, const SBLoc*>> kdtData;
    kdtData.reserve(locData->size());
    std::transform(locData->begin(), locData->end(),std::back_inserter(kdtData),
                   [&](const SBLoc& l)->std::pair<Point<3>,
                   const SBLoc*>{return {l.locToCart3DPt(), &l};});
    locKdt = KDTree<3, const SBLoc*, Point<3>::DistType::EUC>(kdtData.begin(),
                                                              kdtData.end());
}
