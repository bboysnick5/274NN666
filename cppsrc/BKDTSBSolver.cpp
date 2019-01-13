//
//  BKDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include <utility>
#include "BKDTSBSolver.hpp"




template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
void BKDTSBSolver<Tree>::build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    generateKDT(locData);
}

template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
void BKDTSBSolver<Tree>::printSolverInfo() const {
    this->locKdt.printTreeInfo();
}


template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
void BKDTSBSolver<Tree>::generateKDT(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    std::vector<std::pair<Point<3>, const SBLoc*>> kdtData;
    kdtData.reserve(locData->size());
    std::transform(locData->begin(), locData->end(),std::back_inserter(kdtData),
                   [&](const SBLoc& l)->std::pair<Point<3>,
                   const SBLoc*>{return {l.locToCart3DPt(), &l};});
    this->locKdt = Tree<3, const SBLoc*, Point<3>::DistType::EUC>(kdtData.begin(),
                                                              kdtData.end());
}


template class BKDTSBSolver<KDTree>;
template class BKDTSBSolver<KDTreeCusMem>;
