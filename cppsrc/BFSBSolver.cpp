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
void BFSBSolver<FPType>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    loc_data_span_ = loc_data_span;
}

template <typename FPType>
const SBLoc<FPType>* BFSBSolver<FPType>::FindNearestLoc(PointND<FPType, 2> geo_search_pt) const {
    return &*Utility::MinElementGivenDistFunc_p(loc_data_span_.rbegin(), loc_data_span_.rend(),
                                                [geo_search_pt](const SBLoc<FPType>& l) {return l.havDistComp(geo_search_pt);},
                                                std::less<FPType>());
    //return &*oneapi::dpl::min_element(oneapi::dpl::execution::par_unseq, locData->cbegin(), locData->cend(), [&](const SBLoc<FPType> &l1, const SBLoc<FPType> &l2){return l1.havDistComp(geo_search_pt) < l2.havDistComp(geo_search_pt);});
}

template <typename FPType>
void BFSBSolver<FPType>::PrintSolverInfo() const {
    std::cout << "This is brute force solver using haversine distance metric.\n";
}


template class BFSBSolver<double>;
template class BFSBSolver<float>;
