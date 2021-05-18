//
//  BFSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/11/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "BFSBSolver.hpp"
#include "Utility.hpp"

template <typename dist_type>
void BFSBSolver<dist_type>::build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    this->locData = locData;
}

template <typename dist_type>
const SBLoc<dist_type>* BFSBSolver<dist_type>::findNearest(const Point<dist_type, 2>& geoSearchPt) const {
    //return &*std::min_element(locData->cbegin(), locData->cend(),
      //                        [=](const SBLoc<dist_type>& l1, const SBLoc<dist_type>& l2) {
        //                          return l1.havDist(geoSearchPt) < l2.havDist(geoSearchPt);});
        
    return &*custom_min_element(locData->cbegin(), locData->cend(),
                                [&geoSearchPt](const SBLoc<dist_type>& l) {return l.havDistComp(geoSearchPt);},
                                std::less<dist_type>());
}

template <typename dist_type>
void BFSBSolver<dist_type>::printSolverInfo() const {
    std::cout << "This is brute force solver using haversine distance metric.\n";
}


template class BFSBSolver<double>;
template class BFSBSolver<float>;
