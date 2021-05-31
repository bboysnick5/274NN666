//
//  UnionUniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/19/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniLatLngBKDTGridSBSolver.hpp"
#include "Utility.hpp"
#include <omp.h>
#include <memory>
#include <thread>


// TO DO:
// 2. Time per search, store each search time in the test loc data;

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type, Def::Threading_Policy policy>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
UnionUniLatLngBKDTGridSBSolver(dist_type alpc, size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType, dist_type>(), AVE_LOC_PER_CELL(alpc),
MAX_CACHE_CELL_VEC_SIZE(maxCacheCellVecSize) {}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type, Def::Threading_Policy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::printSolverInfo() const {
    std::size_t numSingleLocs = 0, numVecLocs = 0, numTreeLocs = 0, numUniqueVecLocs = 0;
    std::size_t totalTreeNodeSize = 0, totalVecLocSize = 0;
    std::for_each(gridCache.cbegin(), gridCache.cend(), [&](const BitCell& cell) mutable {
        auto cellPtr = cell.getPtr();
        if (auto cellSize = cell.size(cellPtr);
            cellSize == 1) {
            numSingleLocs++;
        } else if (cellSize < MAX_CACHE_CELL_VEC_SIZE) {
            numVecLocs++;
            totalVecLocSize += cellSize;
            if (cell.isUniqueVecLoc(cellPtr))
                numUniqueVecLocs++;
        } else {
            numTreeLocs++;
            totalTreeNodeSize += cellSize;
        }
    });
    std::size_t totalCacheLocs = numSingleLocs + totalVecLocSize + totalTreeNodeSize;
    
    std::cout << "Total cached locs: " << totalCacheLocs << std::endl
    << "Ratio of cache locs over actual num locs: " << static_cast<double>(totalCacheLocs)/numGivenLocs << std::endl
    << "Total num of loc cells: " << gridCache.size() << std::endl
    << "Ave cached locs per cell: " << static_cast<double>(totalCacheLocs)/gridCache.size() << std::endl
    << "Single loc cells: " << numSingleLocs << std::endl
    << "Vector loc cells: " << numVecLocs << std::endl
    << "Unique Vector loc cells: " << numUniqueVecLocs << std::endl
    << "Ave vec loc size: " << static_cast<double>(totalVecLocSize)/numVecLocs << std::endl
    << "Kd-tree loc cells: " << numTreeLocs << std::endl
    << "Ave tree size : " << static_cast<double>(totalTreeNodeSize)/numTreeLocs << std::endl
    << "Ave tree height: " << log2(static_cast<double>(totalTreeNodeSize)/numTreeLocs + 1.0) + 1.0 << std::endl;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType,
          class dist_type, Def::Threading_Policy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
loopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs, Policy_Tag<Def::Threading_Policy::SINGLE>) {
    gridCache.reserve(rowSize*colSize);
    dist_type lat1 = -0.5*Def::PI<dist_type>;
    dist_type thisCtrLat = 0.5 * (latInc - Def::PI<dist_type>);
    for (size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc, lat1 += latInc) {
        dist_type thisCtrLng = 0.5 * lngInc - Def::PI<dist_type>;
        dist_type diagonalDistSq3DEUC = SBLoc<dist_type>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
        for (size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            fillCacheCell({thisCtrLat, thisCtrLng}, diagonalDistSq3DEUC, colSize, ptLocPairs);
        }
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType,
          class dist_type, Def::Threading_Policy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
loopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs, Policy_Tag<Def::Threading_Policy::MULTI_OMP>) {
    gridCache.resize(rowSize*colSize, 0);
    dist_type initCtrLat = 0.5*latInc - 0.5*Def::PI<dist_type>;
    dist_type initCtrLng = 0.5*lngInc - Def::PI<dist_type>;
    dist_type initLat1 = - 0.5*Def::PI<dist_type>;
#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
shared(this->gridCache, lngInc, latInc, initCtrLat, initCtrLng, initLat1) \
firstprivate(ptLocPairs) default(none) schedule(static, 1) collapse(2) \
ordered
    for (size_t r = 0; r < rowSize; ++r) {
        dist_type lat1 = r*latInc + initLat1;
        dist_type diagonalDistSq3DEUC = SBLoc<dist_type>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
        dist_type thisCtrLat = initCtrLat + r*latInc;
        std::size_t startIdxThisRow = r*colSize;
        for (size_t c = 0; c < colSize; ++c) {
            fillCacheCell(startIdxThisRow + c, {thisCtrLat, initCtrLng + c*lngInc},
                          diagonalDistSq3DEUC, colSize, ptLocPairs);
        }
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType,
          class dist_type, Def::Threading_Policy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2.0*Def::PI<dist_type>/colSize + std::numeric_limits<dist_type>::epsilon();
    lngIncInverse = 1.0/lngInc;
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    loopBodyThreadingPolicyDispatch(ptLocPairs);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType,
          class dist_type, Def::Threading_Policy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
loopBodyThreadingPolicyDispatch(std::vector<typename KDT<KDTType, dist_type>::node_type> &ptLocPairs) {
    if constexpr (policy == Def::Threading_Policy::SINGLE) {
        loopBody(ptLocPairs, Policy_Tag<Def::Threading_Policy::SINGLE>{});
    } else if (policy == Def::Threading_Policy::MULTI_OMP) {
        loopBody(ptLocPairs, Policy_Tag<Def::Threading_Policy::MULTI_OMP>{});
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type, Def::Threading_Policy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::calcSideLenFromAlpc() {
    dist_type surfaceArea = 4.0*Def::PI<dist_type>*SBLoc<dist_type>::EARTH_RADIUS*SBLoc<dist_type>::EARTH_RADIUS;
    dist_type numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type, Def::Threading_Policy policy>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    numGivenLocs = locData->size();
    BKDTSBSolver<KDTType, dist_type>::generateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc<dist_type>::deltaLatOnSameLngFromHavDist(sideLen));
    latIncInverse = 1.0/latInc;
    rowSize = std::ceil(Def::PI<dist_type>/latInc);
    fillGridCache();
    //auto t = std::thread(&KDT<KDTType, dist_type>::clear, this->locKdt);
    //t.detach();
    this->locKdt = {};
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type, Def::Threading_Policy policy>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
returnNNLocFromCacheVariant(const Point<dist_type, 2>& geoPt, const BitCell& cell) const {
    if (uintptr_t ptr = cell.getPtr();
        cell.size(ptr) == 1) {
        return cell.getSingleLoc(ptr);
    } else if (cell.size(ptr) < MAX_CACHE_CELL_VEC_SIZE) {
        return Utility::custom_min_element(cell.getLocPairs(ptr), cell.getLocPairs(ptr) + cell.size(ptr),
                                           [&](const auto& nh) {return SBLoc<dist_type>::geoPtToCart3DPt(geoPt).template
                                                dist<Point<dist_type, 3>::DistType::EUCSQ>(nh.key);},
                                           std::less())->value;
    } else {
        return cell.getCacheTree(ptr)->kNNValue(SBLoc<dist_type>::geoPtToCart3DPt(geoPt), 1);
    }
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type, Def::Threading_Policy policy>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
findNearest(const Point<dist_type, 2>& geoSearchPt) const {
    return returnNNLocFromCacheVariant(geoSearchPt,
           gridCache[static_cast<size_t>((geoSearchPt[0]+0.5*Def::PI<dist_type>)*latIncInverse)*colSize +
                     static_cast<size_t>((geoSearchPt[1]+Def::PI<dist_type>)*lngIncInverse)]);
}



template class UnionUniLatLngBKDTGridSBSolver<KDTree, double, Def::Threading_Policy::SINGLE>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, float, Def::Threading_Policy::SINGLE>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, double, Def::Threading_Policy::MULTI_OMP>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, float, Def::Threading_Policy::MULTI_OMP>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, double, Def::Threading_Policy::MULTI_HAND>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, float, Def::Threading_Policy::MULTI_HAND>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double, Def::Threading_Policy::SINGLE>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float, Def::Threading_Policy::SINGLE>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double, Def::Threading_Policy::MULTI_OMP>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float, Def::Threading_Policy::MULTI_OMP>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double, Def::Threading_Policy::MULTI_HAND>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float, Def::Threading_Policy::MULTI_HAND>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double, Def::Threading_Policy::SINGLE>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float, Def::Threading_Policy::SINGLE>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double, Def::Threading_Policy::MULTI_OMP>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float, Def::Threading_Policy::MULTI_OMP>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double, Def::Threading_Policy::MULTI_HAND>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float, Def::Threading_Policy::MULTI_HAND>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double, Def::Threading_Policy::SINGLE>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float, Def::Threading_Policy::SINGLE>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double, Def::Threading_Policy::MULTI_OMP>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float, Def::Threading_Policy::MULTI_OMP>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double, Def::Threading_Policy::MULTI_HAND>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float, Def::Threading_Policy::MULTI_HAND>;







