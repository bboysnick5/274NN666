//
//  KDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "KDTSBSolver.hpp"
#include <algorithm>


template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
void KDTSBSolver<Tree>::build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    generateKDT(locData);
}

template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
void KDTSBSolver<Tree>::printSolverInfo() const {
    locKdt.printTreeInfo();
}


template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
const SBLoc* KDTSBSolver<Tree>::findNearest(double lng, double lat) const {
    return locKdt.kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}

template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
void KDTSBSolver<Tree>::generateKDT(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    std::for_each(locData->begin(), locData->end(), [&](const SBLoc &loc){
        locKdt.insert(loc.locToCart3DPt(), &loc);});
}


template class KDTSBSolver<KDTree>;
template class KDTSBSolver<KDTreeCusMem>;


