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

template <typename dist_type>
void BFSBSolver<dist_type>::Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    this->locData = locData;
}

template <typename dist_type>
const SBLoc<dist_type>* BFSBSolver<dist_type>::FindNearestLoc(const PointND<dist_type, 2>& geoSearchPt) const {
    return &*Utility::MinElementGivenDistFunc_p(locData->cbegin(), locData->cend(),
                                         [&geoSearchPt](const SBLoc<dist_type>& l) {return l.havDistComp(geoSearchPt);},
                                         std::less<dist_type>());
    //return &*oneapi::dpl::min_element(oneapi::dpl::execution::par_unseq, locData->cbegin(), locData->cend(), [&](const SBLoc<dist_type> &l1, const SBLoc<dist_type> &l2){return l1.havDistComp(geoSearchPt) < l2.havDistComp(geoSearchPt);});
}

template <typename dist_type>
void BFSBSolver<dist_type>::PrintSolverInfo() const {
    std::cout << "This is brute force solver using haversine distance metric.\n";
}


template class BFSBSolver<double>;
template class BFSBSolver<float>;
