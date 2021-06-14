//
//  BFSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/11/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//


#include "BFSBSolver.hpp"
#include "Utility.hpp"
//#include <oneapi/tbb.h>

template <typename FPType>
void BFSBSolver<FPType>::Build(const std::shared_ptr<std::vector<SBLoc<FPType>>> &locData) {
    this->locData = locData;
}

template <typename FPType>
const SBLoc<FPType>* BFSBSolver<FPType>::FindNearestLoc(const PointND<FPType, 2>& geoSearchPt) const {
    return &*Utility::MinElementGivenDistFunc_p(locData->cbegin(), locData->cend(),
                                         [&geoSearchPt](const SBLoc<FPType>& l) {return l.havDistComp(geoSearchPt);},
                                         std::less<FPType>());
    //return &*oneapi::dpl::min_element(oneapi::dpl::execution::par_unseq, locData->cbegin(), locData->cend(), [&](const SBLoc<FPType> &l1, const SBLoc<FPType> &l2){return l1.havDistComp(geoSearchPt) < l2.havDistComp(geoSearchPt);});
}

template <typename FPType>
void BFSBSolver<FPType>::PrintSolverInfo() const {
    std::cout << "This is brute force solver using haversine distance metric.\n";
}


template class BFSBSolver<double>;
template class BFSBSolver<float>;
