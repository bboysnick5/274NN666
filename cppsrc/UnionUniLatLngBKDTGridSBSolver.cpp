//
//  UnionUniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/19/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniLatLngBKDTGridSBSolver.hpp"
#include <memory>

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::UnionCell::
UnionCell(size_t maxCacheVecSize, const std::vector<std::pair<Point<dist_type, 3>,
          const SBLoc<dist_type>*>>& bufVec) : _size(bufVec.size()) {
    if (_size == 1) {
        cacheLoc = bufVec[0].second;
    } else if (_size < maxCacheVecSize) {
        //cacheLocs = static_cast<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>*>(::operator new[](_size, std::nothrow));
        //std::uninitialized_move(bufVec.cbegin(), bufVec.cend(), cacheLocs);
        cacheLocs = new std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>[_size];
        std::copy(bufVec.cbegin(), bufVec.cend(), cacheLocs);
    } else {
        _size = 0;
        cacheTree = new KDT<KDTType, dist_type>(bufVec.begin(), bufVec.end());
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::UnionCell::~UnionCell() {
    if (_size > 1) {
        delete cacheLocs;
    } else if (_size == 0) {
        delete cacheTree;
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
size_t UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::UnionCell::size() const {
    return _size != 0 ? _size : cacheTree->size();
}



template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>

UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::
BitCell(size_t maxCacheVecSize,
        const std::vector<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>> &bufVec) {
    size_t size = bufVec.size();
    const uintptr_t MASK = ~(1ULL << 48);
    if (size == 1) {
        ptr = (reinterpret_cast<std::uintptr_t>(bufVec[0].second) & MASK) | (1ull << 48);
    } else if (size < maxCacheVecSize) {
        std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*> *cacheLocs = static_cast<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>*>(::operator new(size*sizeof(std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>), std::nothrow));
        std::uninitialized_move(bufVec.cbegin(), bufVec.cend(), cacheLocs);
        //std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*> *cacheLocs = new std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>[size];
        //std::copy(bufVec.cbegin(), bufVec.cend(), cacheLocs);
        ptr = (reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK) | (size << 48);
    } else {
        ptr = reinterpret_cast<std::uintptr_t>(new KDT<KDTType, dist_type>(bufVec.begin(), bufVec.end())) & MASK;
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::~BitCell() {
    auto size = ptr >> 48;
    if (size > 1) {
        ::operator delete(const_cast<void*>(reinterpret_cast<const void*>(getLocPairs())));
    } else if (size == 0) {
        delete getCacheTree();
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
size_t UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::size() const {
    auto size = ptr >> 48;
    return size != 0 ? size : getCacheTree()->size();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
size_t UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::rawSize() const {
    return ptr >> 48;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::getSingleLoc() const {
    return reinterpret_cast<const SBLoc<dist_type>*>(static_cast<intptr_t>(ptr << 16) >> 16);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::getLocPairs() const {
    return reinterpret_cast<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>*>(static_cast<intptr_t>(ptr << 16) >> 16);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const KDT<KDTType, dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::getCacheTree() const {
    return reinterpret_cast<const KDT<KDTType, dist_type>*>(static_cast<intptr_t>(ptr << 16) >> 16);
}



template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
UnionUniLatLngBKDTGridSBSolver(dist_type alpc, size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType, dist_type>(), AVE_LOC_PER_CELL(alpc),
MAX_CACHE_CELL_VEC_SIZE(maxCacheCellVecSize) {}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::printSolverInfo() const {
    std::cout << "Total cache locs: " << totalNodeSize
    << "\nAve tree size: " << totalNodeSize/gridCache.size()
    << "\nAve tree height: "
    << static_cast<size_t>(log2(totalNodeSize/gridCache.size() + 1)) + 1
    << "\nRatio of cache locs over actual num locs: "
    << totalNodeSize/totalLocSize
    << "\nTotal num of loc cells: " << gridCache.size()
    << "\nSingle loc cells: " << singleLocs
    << "\nVector loc cells: " << vecLocs
    << "\nKd-tree loc cells: " << gridCache.size() -singleLocs -vecLocs << "\n";
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
fillCacheCell(dist_type thisCtrLng, dist_type thisCtrLat, dist_type thisDiff,
              std::vector<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::latLngToCart3DPt(thisCtrLat, thisCtrLng),
                                   thisDiff, std::back_inserter(ptLocPairs));
    size_t vecSize = ptLocPairs.size();
    this->totalNodeSize += vecSize;
    if (vecSize == 1) {
        ++singleLocs;
    } else if (vecSize < MAX_CACHE_CELL_VEC_SIZE) {
        ++vecLocs;
    }
    this->gridCache.emplace_back(MAX_CACHE_CELL_VEC_SIZE, ptLocPairs);
    ptLocPairs.clear();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2.0*M_PI/colSize + 2.0*M_PI/(colSize*colSize*0xFFFF);
    lngIncInverse = 1.0/lngInc;
    gridCache.reserve(this->locKdt.size()*1.2/AVE_LOC_PER_CELL);
    std::vector<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    
    dist_type thisCtrLat = 0.5 * (latInc - M_PI);
    for (size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc) {
        dist_type thisCtrLng = 0.5 * lngInc - M_PI;
        dist_type thisDiff = SBLoc<dist_type>::xyzDistFromLngLat(r*latInc- 0.5*M_PI,
                                                   (r+1)*latInc-0.5*M_PI, lngInc);
        for (size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            fillCacheCell(thisCtrLng, thisCtrLat, thisDiff, ptLocPairs);
        }
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::calcSideLenFromAlpc() {
    dist_type surfaceArea = 4*M_PI*SBLoc<dist_type>::EARTH_RADIUS*SBLoc<dist_type>::EARTH_RADIUS;
    dist_type numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    totalLocSize = locData->size();
    BKDTSBSolver<KDTType, dist_type>::generateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc<dist_type>::latFromHavDist(sideLen, 0.0));
    latIncInverse = 1.0/latInc;
    rowSize = std::ceil(M_PI/(latInc - latInc*latInc/(M_PI*0xFFFF)));
    fillGridCache();
    this->locKdt.clear();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
returnNNLocFromCacheVariant(dist_type lat, dist_type lng, const BitCell& cell) const {
    if (cell.size() == 1) {
        return cell.getSingleLoc();
    } else if (cell.size() < MAX_CACHE_CELL_VEC_SIZE) {
        const auto p = SBLoc<dist_type>::latLngToCart3DPt(lat, lng);
        return std::min_element(cell.getLocPairs(), cell.getLocPairs() + cell.size(),
                                [&](const auto& p1, const auto& p2){ return Point<dist_type, 3>::template
                                    dist<Point<dist_type, 3>::DistType::EUCSQ>(p1.first, p) <
                                    Point<dist_type, 3>::template dist<Point<dist_type, 3>::DistType::EUCSQ>(p2.first, p);
                                })->second;
    } else {
        return cell.getCacheTree()->kNNValue(SBLoc<dist_type>::latLngToCart3DPt(lat, lng), 1);
    }
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
findNearest(dist_type lat, dist_type lng) const {
    return returnNNLocFromCacheVariant(lat, lng, gridCache[static_cast<size_t>(
        (lat+0.5*M_PI)*latIncInverse)*colSize+ static_cast<size_t>((lng+M_PI)*lngIncInverse)]);
}



template class UnionUniLatLngBKDTGridSBSolver<KDTree, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, float>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float>;






