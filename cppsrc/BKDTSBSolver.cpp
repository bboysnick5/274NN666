//
//  BKDTSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include <utility>
#include "BKDTSBSolver.hpp"




template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void BKDTSBSolver<KDTType, dist_type>::build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    generateKDT(locData);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void BKDTSBSolver<KDTType, dist_type>::printSolverInfo() const {
    this->locKdt.printTreeInfo();
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void BKDTSBSolver<KDTType, dist_type>::generateKDT(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    /*
    std::vector<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>> kdtData;
    kdtData.reserve(locData->size());
    std::array<dist_type, 6> minMaxAllDim;
    std::generate(minMaxAllDim.begin(), minMaxAllDim.end(),
                  [lowHighToggle = 0u,
                  lowHigh = std::array<dist_type, 2>{std::numeric_limits<dist_type>::max(),
                  std::numeric_limits<dist_type>::min()}]() mutable {
                  return lowHigh[lowHighToggle++%2];});
    std::transform(locData->cbegin(),locData->cend(), std::back_inserter(kdtData),
                   [&](const SBLoc<dist_type>& l) mutable ->std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>{
                       const auto p = l.locToCart3DPt();
                       for (size_t i = 0, dim = 0; i < 6; i += 2, ++dim) {
                           minMaxAllDim[i] = std::min(minMaxAllDim[i], p[dim]);
                           minMaxAllDim[i+1] = std::max(minMaxAllDim[i+1], p[dim]);
                       }
                       return {p, &l};});
    std::array<dist_type, 3> maxSpansAllDim;
    std::generate(maxSpansAllDim.begin(), maxSpansAllDim.end(),
                  [n=0, &minMaxAllDim]() mutable {
                      n += 2;
                      return minMaxAllDim[n-1] - minMaxAllDim[n-2];});
    this->locKdt = KDT<KDTType>(kdtData.begin(), kdtData.end(), maxSpansAllDim);
    */
    
    
    std::vector<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>> kdtData;
    kdtData.reserve(locData->size());
    std::transform(locData->cbegin(),locData->cend(), std::back_inserter(kdtData),
                   [&](const SBLoc<dist_type>& l) mutable ->std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>{
                       return {l.locToCart3DPt(), &l};});
    this->locKdt = KDT<KDTType, dist_type>(kdtData.begin(), kdtData.end());
}


template class BKDTSBSolver<KDTree, double>;
template class BKDTSBSolver<KDTree, float>;

template class BKDTSBSolver<KDTreeCusMem, double>;
template class BKDTSBSolver<KDTreeCusMem, float>;

template class BKDTSBSolver<KDTreeExpandLongest, double>;
template class BKDTSBSolver<KDTreeExpandLongest, float>;

template class BKDTSBSolver<KDTreeExpandLongestVec, double>;
template class BKDTSBSolver<KDTreeExpandLongestVec, float>;
