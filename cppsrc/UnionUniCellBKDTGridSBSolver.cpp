//
//  UnionUniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/20/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniCellBKDTGridSBSolver.hpp"
#include <thread>
#include <map>
#include <algorithm>


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
UnionUniCellBKDTGridSBSolver<KDTType, dist_type, policy>::
UnionUniCellBKDTGridSBSolver(dist_type alpc, size_t maxCacheCellVecSize)
: UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>(alpc, maxCacheCellVecSize) {}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
void UnionUniCellBKDTGridSBSolver<KDTType, dist_type, policy>::
loopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs,
              typename UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::template Policy_Tag<def::ThreadingPolicy::kSingle>) {
    this->gridCache.reserve(totalCacheCells);
    dist_type thisCtrLat = 0.5*this->latInc - 0.5*def::kMathPi<dist_type>;
    dist_type lat1 = - 0.5*def::kMathPi<dist_type>;
    for (size_t r = 0; r < this->rowSize; ++r, thisCtrLat += this->latInc, lat1 += this->latInc) {
        auto &[thisColSize, cosThisLngInc] = colSizeCosLngIncEachRowVec[r];
        dist_type thisLngInc = 1.0/thisRowStartIdxThisLngIncInverseVec[r].second;
        dist_type diagonalDistSq3DEUC = UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::EUC3DDistSqFromLatCosDeltaLng(lat1, lat1 + this->latInc, cosThisLngInc);
        dist_type thisCtrLng = 0.5 * thisLngInc - def::kMathPi<dist_type>;
        for (std::size_t c = 0; c < thisColSize; ++c, thisCtrLng += thisLngInc) {
            this->fillCacheCell({thisCtrLat, thisCtrLng}, diagonalDistSq3DEUC, thisColSize, ptLocPairs);
        }
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
void UnionUniCellBKDTGridSBSolver<KDTType, dist_type, policy>::ompLoopCollapsePrep() {
    
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
void UnionUniCellBKDTGridSBSolver<KDTType, dist_type, policy>::
loopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs,
         typename UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
         template Policy_Tag<def::ThreadingPolicy::kMultiOmp>) {
    this->gridCache.resize(totalCacheCells, 0);
    dist_type initCtrLat = 0.5*this->latInc - 0.5*def::kMathPi<dist_type>;
    dist_type initLat1 = - 0.5*def::kMathPi<dist_type>;
#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
shared(this->gridCache, thisRowStartIdxThisLngIncInverseVec, colSizeCosLngIncEachRowVec, \
       this->latInc, initCtrLat, initLat1) \
firstprivate(ptLocPairs) default(none) schedule(dynamic, 1) \
//ordered
    for (std::size_t idx = 0; idx < totalCacheCells; ++idx) {
        const auto it = std::upper_bound(thisRowStartIdxThisLngIncInverseVec.cbegin(), thisRowStartIdxThisLngIncInverseVec.cend(), idx, [](const std::size_t &idx, const std::pair<std::size_t, dist_type> &p){return idx < p.first;}) - 1;
        std::size_t r = it - thisRowStartIdxThisLngIncInverseVec.begin();
        std::size_t c = idx - it->first;
        dist_type thisLngInc = 1.0/it->second;
        auto &[thisColSize, cosThisLngInc] = colSizeCosLngIncEachRowVec[r];
        dist_type lat1 = r*this->latInc + initLat1;
        dist_type thisCtrLat = initCtrLat + r*this->latInc;
        dist_type diagonalDistSq3DEUC = UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::EUC3DDistSqFromLatCosDeltaLng(lat1, lat1 + this->latInc, cosThisLngInc);
        dist_type initThisCtrLng = 0.5 * thisLngInc - def::kMathPi<dist_type>;
        this->fillCacheCell(idx, {thisCtrLat, initThisCtrLng + c*thisLngInc}, diagonalDistSq3DEUC, thisColSize, ptLocPairs);
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
void UnionUniCellBKDTGridSBSolver<KDTType, dist_type, policy>::fillGridCache() {
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(this->MAX_CACHE_CELL_VEC_SIZE);
    colSizeCosLngIncEachRowVec.reserve(this->rowSize);
    thisRowStartIdxThisLngIncInverseVec.reserve(this->rowSize);
    dist_type earthPerimeterOverSideLen = 2.0*def::kMathPi<dist_type> * SBLoc<dist_type>::EARTH_RADIUS /this->sideLen;
    dist_type thisLat = -0.5*def::kMathPi<dist_type>;
    totalCacheCells = 0;
    for (std::size_t r = 0; r < this->rowSize - 1; ++r, thisLat += this->latInc) {
        std::size_t thisColSize = static_cast<size_t>(earthPerimeterOverSideLen * cos(thisLat)) + 1;
        dist_type thisLngInc = 2.0*def::kMathPi<dist_type>/thisColSize + std::numeric_limits<dist_type>::epsilon();
        colSizeCosLngIncEachRowVec.emplace_back(thisColSize, std::cos(thisLngInc));
        thisRowStartIdxThisLngIncInverseVec.emplace_back(totalCacheCells, 1.0/thisLngInc);
        totalCacheCells += thisColSize;
    }
    colSizeCosLngIncEachRowVec.emplace_back(colSizeCosLngIncEachRowVec.back());
    thisRowStartIdxThisLngIncInverseVec.emplace_back(totalCacheCells, thisRowStartIdxThisLngIncInverseVec.back().second);
    totalCacheCells += colSizeCosLngIncEachRowVec.back().first;
    UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::loopBodyThreadingPolicyDispatch(ptLocPairs);
    colSizeCosLngIncEachRowVec.clear();
    colSizeCosLngIncEachRowVec.shrink_to_fit();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
const SBLoc<dist_type>* UnionUniCellBKDTGridSBSolver<KDTType, dist_type, policy>::
findNearest(const Point<dist_type, 2>& geoPt) const {
  //  return UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::findNearest(geoPt);
    const auto &[startIdx, thisLngIncInverse] = thisRowStartIdxThisLngIncInverseVec[(geoPt[0]+0.5*def::kMathPi<dist_type>)*this->latIncInverse];
    return UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::returnNNLocFromCacheVariant(geoPt,
        this->gridCache[startIdx + static_cast<size_t>((geoPt[1]+def::kMathPi<dist_type>)*thisLngIncInverse)]);
}


template class UnionUniCellBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTree, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniCellBKDTGridSBSolver<KDTree, float, def::ThreadingPolicy::kMultiHand>;

template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniCellBKDTGridSBSolver<KDTreeCusMem, float, def::ThreadingPolicy::kMultiHand>;

template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, float, def::ThreadingPolicy::kMultiHand>;

template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kSingle>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiOmp>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, double, def::ThreadingPolicy::kMultiHand>;
template class UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, float, def::ThreadingPolicy::kMultiHand>;
