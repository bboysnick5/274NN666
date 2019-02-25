//
//  BFSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/11/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include <algorithm>
#include "BFSBSolver.hpp"

template <typename dist_type>
void BFSBSolver<dist_type>::build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    this->locData = locData;
}

template <typename dist_type>
const SBLoc<dist_type>* BFSBSolver<dist_type>::findNearest(dist_type lat, dist_type lng) const {
    return &*std::min_element(locData->cbegin(), locData->cend(),
                              [=](const SBLoc<dist_type>& l1, const SBLoc<dist_type>& l2) {
                                  return SBLoc<dist_type>::havDist(l1.lng,l1.lat, lng, lat)
                                  < SBLoc<dist_type>::havDist(l2.lng, l2.lat, lng, lat);});
}

template <typename dist_type>
void BFSBSolver<dist_type>::printSolverInfo() const {
    std::cout << "This is brute force solver using haversine distance metric.\n";
}


template class BFSBSolver<double>;
template class BFSBSolver<float>;
