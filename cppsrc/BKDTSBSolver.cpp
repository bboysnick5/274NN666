//
//  BKDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include <utility>
#include "BKDTSBSolver.hpp"




template <template <class value_type, size_t, class, typename Point<value_type, 3>::DistType> class KDTType>
void BKDTSBSolver<KDTType>::build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    generateKDT(locData);
}

template <template <class value_type, size_t, class, typename Point<value_type, 3>::DistType> class KDTType>
void BKDTSBSolver<KDTType>::printSolverInfo() const {
    this->locKdt.printTreeInfo();
}


template <template <class value_type, size_t, class, typename Point<value_type, 3>::DistType> class KDTType>
void BKDTSBSolver<KDTType>::generateKDT(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    /*
    std::vector<std::pair<Point<double, 3>, const SBLoc*>> kdtData;
    kdtData.reserve(locData->size());
    std::array<double, 6> minMaxAllDim;
    std::generate(minMaxAllDim.begin(), minMaxAllDim.end(),
                  [lowHighToggle = 0u,
                  lowHigh = std::array<double, 2>{std::numeric_limits<double>::max(),
                  std::numeric_limits<double>::min()}]() mutable {
                  return lowHigh[lowHighToggle++%2];});
    std::transform(locData->cbegin(),locData->cend(), std::back_inserter(kdtData),
                   [&](const SBLoc& l) mutable ->std::pair<Point<double, 3>, const SBLoc*>{
                       const auto p = l.locToCart3DPt();
                       for (size_t i = 0, dim = 0; i < 6; i += 2, ++dim) {
                           minMaxAllDim[i] = std::min(minMaxAllDim[i], p[dim]);
                           minMaxAllDim[i+1] = std::max(minMaxAllDim[i+1], p[dim]);
                       }
                       return {p, &l};});
    std::array<double, 3> maxSpansAllDim;
    std::generate(maxSpansAllDim.begin(), maxSpansAllDim.end(),
                  [n=0, &minMaxAllDim]() mutable {
                      n += 2;
                      return minMaxAllDim[n-1] - minMaxAllDim[n-2];});
    this->locKdt = KDT<KDTType>(kdtData.begin(), kdtData.end(), maxSpansAllDim);
    */
    
    
    
    std::vector<std::pair<Point<double, 3>, const SBLoc*>> kdtData;
    kdtData.reserve(locData->size());
    std::transform(locData->cbegin(),locData->cend(), std::back_inserter(kdtData),
                   [&](const SBLoc& l) mutable ->std::pair<Point<double, 3>, const SBLoc*>{
                       return {l.locToCart3DPt(), &l};});
    this->locKdt = KDT<KDTType>(kdtData.begin(), kdtData.end()); 
}


template class BKDTSBSolver<KDTree>;
template class BKDTSBSolver<KDTreeCusMem>;
template class BKDTSBSolver<KDTreeExpandLongest>;
template class BKDTSBSolver<KDTreeExpandLongestVec>;
