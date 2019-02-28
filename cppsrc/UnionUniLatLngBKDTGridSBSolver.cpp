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
UnionCell(const std::vector<typename KDT<KDTType, dist_type>::node_type>& bufVec, size_t maxCacheVecSize) : _size(bufVec.size()) {
    if (_size == 1) {
        cacheLoc = bufVec[0].value;
    } else if (_size < maxCacheVecSize) {
        //cacheLocs = static_cast<typename KDT<KDTType, dist_type>::node_type*>(::operator new[](_size, std::nothrow));
        //std::uninitialized_move(bufVec.cbegin(), bufVec.cend(), cacheLocs);
        cacheLocs = new typename KDT<KDTType, dist_type>::node_type[_size];
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
BitCell(const std::vector<typename KDT<KDTType, dist_type>::node_type> &bufVec,
        size_t maxCacheVecSize, const BitCell& prevCell, const BitCell& upCell) {
    size_t size = bufVec.size();
    const uintptr_t MASK = ~(1ULL << 48);
    if (size == 1) {
        ptr = (reinterpret_cast<std::uintptr_t>(bufVec[0].value) & MASK) | (1ull << 48);
    } else if (size < maxCacheVecSize) {
        for (const auto &cell : {&prevCell, &upCell}) {
            if (size == cell->size() &&
                std::equal(bufVec.cbegin(), bufVec.cend(), cell->getLocPairs(),
                           [](const auto &nh1, const auto &nh2){return nh1.value == nh2.value;})) {
                ptr = cell->ptr & (std::numeric_limits<uintptr_t>::max() - 1);
                return;
            }
        }
        auto *cacheLocs = new typename KDT<KDTType, dist_type>::node_type[size];
        std::copy(bufVec.cbegin(), bufVec.cend(), cacheLocs);
        ptr = (reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK) | (size << 48) | 1ull;
    } else {
        ptr = reinterpret_cast<std::uintptr_t>(new KDT<KDTType, dist_type>(bufVec.begin(), bufVec.end())) & MASK;
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::
BitCell(const std::vector<typename KDT<KDTType, dist_type>::node_type> &bufVec,
        size_t maxCacheVecSize, const BitCell& prevCell) {
    size_t size = bufVec.size();
    const uintptr_t MASK = ~(1ULL << 48);
    if (size == 1) {
        ptr = (reinterpret_cast<std::uintptr_t>(bufVec[0].value) & MASK) | (1ull << 48);
    } else if (size < maxCacheVecSize) {
        if (size == prevCell.size() &&
            std::equal(bufVec.cbegin(), bufVec.cend(), prevCell.getLocPairs(),
                       [](const auto &nh1, const auto &nh2){return nh1.value == nh2.value;})) {
            ptr = prevCell.ptr & (std::numeric_limits<uintptr_t>::max() - 1);
            return;
        }
        auto *cacheLocs = new typename KDT<KDTType, dist_type>::node_type[size];
        std::copy(bufVec.cbegin(), bufVec.cend(), cacheLocs);
        ptr = (reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK) | (size << 48) | 1ull;
    } else {
        ptr = reinterpret_cast<std::uintptr_t>(new KDT<KDTType, dist_type>(bufVec.begin(), bufVec.end())) & MASK;
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::
BitCell(const std::vector<typename KDT<KDTType, dist_type>::node_type> &bufVec,
        size_t maxCacheVecSize) {
    size_t size = bufVec.size();
    const uintptr_t MASK = ~(1ULL << 48);
    if (size == 1) {
        ptr = (reinterpret_cast<std::uintptr_t>(bufVec[0].value) & MASK) | (1ull << 48);
    } else if (size < maxCacheVecSize) {
        //typename KDT<KDTType, dist_type>::node_type *cacheLocs = static_cast<typename KDT<KDTType, dist_type>::node_type*>(::operator new(size*sizeof(typename KDT<KDTType, dist_type>::node_type), std::nothrow));
        //std::uninitialized_move(bufVec.cbegin(), bufVec.cend(), cacheLocs);
        auto *cacheLocs = new typename KDT<KDTType, dist_type>::node_type[size];
        std::copy(bufVec.cbegin(), bufVec.cend(), cacheLocs);
        ptr = (reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK) | (size << 48) | 1ull;
    } else {
        ptr = reinterpret_cast<std::uintptr_t>(new KDT<KDTType, dist_type>(bufVec.begin(), bufVec.end())) & MASK;
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::~BitCell() {
    auto size = ptr >> 48;
    if (size > 1 && (ptr & 1ull) == 1ull) {
        delete[] getLocPairs();
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
const typename KDT<KDTType, dist_type>::node_type* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::getLocPairs() const {
    return reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(ptr << 16) >> 16) & std::numeric_limits<uintptr_t>::max() - 1);
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
    // << "\nUnique Vector loc cells: " <<
    << "\nKd-tree loc cells: " << gridCache.size() -singleLocs -vecLocs << "\n";
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
fillCacheCell(const Point<dist_type, 2>& thisCtrGeoPt, dist_type thisDiff, size_t thisColSize,
              std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt(thisCtrGeoPt),
                                   thisDiff, std::back_inserter(ptLocPairs));
    size_t vecSize = ptLocPairs.size();
    this->totalNodeSize += vecSize;
    if (vecSize == 1) {
        ++singleLocs;
    } else if (vecSize < MAX_CACHE_CELL_VEC_SIZE) {
        ++vecLocs;
    }
    // TODO: fill init optimization
    if (this->gridCache.size() >= thisColSize) {
        this->gridCache.emplace_back(ptLocPairs, MAX_CACHE_CELL_VEC_SIZE,
                                     this->gridCache.back(), this->gridCache[this->gridCache.size() - thisColSize]);
    } else if (this->gridCache.size() > 0) {
        this->gridCache.emplace_back(ptLocPairs, MAX_CACHE_CELL_VEC_SIZE, this->gridCache.back());
    } else {
        this->gridCache.emplace_back(ptLocPairs, MAX_CACHE_CELL_VEC_SIZE);
    }
    ptLocPairs.clear();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2.0*M_PI/colSize + 2.0*M_PI/(colSize*colSize*0xFFFF);
    lngIncInverse = 1.0/lngInc;
    gridCache.reserve(this->locKdt.size()*1.2/AVE_LOC_PER_CELL);
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    
    dist_type thisCtrLat = 0.5 * (latInc - M_PI);
    for (size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc) {
        dist_type thisCtrLng = 0.5 * lngInc - M_PI;
        dist_type thisDiff = SBLoc<dist_type>::xyzDistFromLngLat(r*latInc- 0.5*M_PI,
                                                   (r+1)*latInc-0.5*M_PI, lngInc);
        for (size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            fillCacheCell({thisCtrLat, thisCtrLng}, thisDiff, colSize, ptLocPairs);
        }
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::calcSideLenFromAlpc() {
    dist_type surfaceArea = 4.0*M_PI*SBLoc<dist_type>::EARTH_RADIUS*SBLoc<dist_type>::EARTH_RADIUS;
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
returnNNLocFromCacheVariant(const Point<dist_type, 2>& geoPt, const BitCell& cell) const {
    if (cell.size() == 1) {
        return cell.getSingleLoc();
    } else if (cell.size() < MAX_CACHE_CELL_VEC_SIZE) {
        const auto p = SBLoc<dist_type>::geoPtToCart3DPt(geoPt);
        return std::min_element(cell.getLocPairs(), cell.getLocPairs() + cell.size(),
                                [&](const auto& nh1, const auto& nh2){ return Point<dist_type, 3>::template
                                    dist<Point<dist_type, 3>::DistType::EUCSQ>(nh1.key, p) <
                                    Point<dist_type, 3>::template dist<Point<dist_type, 3>::DistType::EUCSQ>(nh2.key, p);
                                })->value;
    } else {
        return cell.getCacheTree()->kNNValue(SBLoc<dist_type>::geoPtToCart3DPt(geoPt), 1);
    }
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
findNearest(const Point<dist_type, 2>& geoSearchPt) const {
    return returnNNLocFromCacheVariant(geoSearchPt, gridCache[static_cast<size_t>(
        (geoSearchPt[0]+0.5*M_PI)*latIncInverse)*colSize+ static_cast<size_t>((geoSearchPt[1]+M_PI)*lngIncInverse)]);
}



template class UnionUniLatLngBKDTGridSBSolver<KDTree, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, float>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float>;






