//
//  KDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "KDTSBSolver.hpp"
#include <algorithm>


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void KDTSBSolver<KDTType, dist_type>::build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    generateKDT(locData);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void KDTSBSolver<KDTType, dist_type>::printSolverInfo() const {
    locKdt.printTreeInfo();
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* KDTSBSolver<KDTType, dist_type>::findNearest(const Point<dist_type, 2>& geoSearchPt) const {
    return locKdt.kNNValue(SBLoc<dist_type>::geoPtToCart3DPt(geoSearchPt), 1);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void KDTSBSolver<KDTType, dist_type>::generateKDT(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    std::for_each(locData->cbegin(), locData->cend(), [&](const SBLoc<dist_type> &loc) mutable {
        this->locKdt.insert(loc.locToCart3DPt(), &loc);});
}


template class KDTSBSolver<KDTree, double>;
template class KDTSBSolver<KDTree, float>;

template class KDTSBSolver<KDTreeCusMem, double>;
template class KDTSBSolver<KDTreeCusMem, float>;

template class KDTSBSolver<KDTreeExpandLongest, double>;
template class KDTSBSolver<KDTreeExpandLongest, float>;

template class KDTSBSolver<KDTreeExpandLongestVec, double>;
template class KDTSBSolver<KDTreeExpandLongestVec, float>;

