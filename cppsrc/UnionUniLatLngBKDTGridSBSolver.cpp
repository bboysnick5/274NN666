//
//  UnionUniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/19/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniLatLngBKDTGridSBSolver.hpp"
#include "Utility.hpp"
//#include <omp.h>
#include <memory>
#include <thread>

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

/*
 template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
 template <typename... PrevBitCells>
 UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::
 BitCell(const std::vector<typename KDT<KDTType, dist_type>::node_type> &bufVec,
         size_t maxCacheVecSize, const PrevBitCells& ...prevCells) {
     size_t size = bufVec.size();
     const uintptr_t MASK = ~(1ULL << 48);
     if (size == 1) {
         ptr = (reinterpret_cast<std::uintptr_t>(bufVec[0].value) & MASK) | (1ull << 48);
     } else if (size < maxCacheVecSize) {
         auto checkPrev = [&bufVec, size, this](const BitCell& cell) {
             if (size == cell.size() &&
                 std::equal(bufVec.cbegin(), bufVec.cend(), cell.getLocPairs(),
                            [](const auto &nh1, const auto &nh2){return nh1.value == nh2.value;})) {
                 ptr = cell.ptr & (std::numeric_limits<uintptr_t>::max() - 1);
                 return;
             }
         };
         (checkPrev(prevCells), ...);
         auto *cacheLocs = new typename KDT<KDTType, dist_type>::node_type[size];
         std::copy(bufVec.cbegin(), bufVec.cend(), cacheLocs);
         ptr = (reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK) | (size << 48) | 1ull;
     } else {
         ptr = reinterpret_cast<std::uintptr_t>(new KDT<KDTType, dist_type>(bufVec.begin(), bufVec.end())) & MASK;
     }
 }

 */


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
    << static_cast<size_t>(log2(totalNodeSize/gridCache.size() + 1.0)) + 1
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
fillCacheCell(const Point<dist_type, 2>& thisCtrGeoPt, dist_type diagonalDist3DEUC, size_t thisColSize,
              std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt(thisCtrGeoPt), diagonalDist3DEUC, std::back_inserter(ptLocPairs));
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
    lngInc = 2.0*std::numbers::pi_v<dist_type>/colSize + 2.0*std::numbers::pi_v<dist_type>/(colSize*65536);
    lngIncInverse = 1.0/lngInc;
    gridCache.reserve(this->locKdt.size()*1.2/AVE_LOC_PER_CELL);
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    
    dist_type thisCtrLat = 0.5 * (latInc - std::numbers::pi_v<dist_type>);
/*#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
default(none) schedule(auto) shared(this->gridCache, lngInc, latInc) firstprivate(thisCtrLat, ptLocPairs) \
reduction(+:totalTreeSize, singleLocs) collapse(2) */
    for (size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc) {
        dist_type thisCtrLng = 0.5 * lngInc - std::numbers::pi_v<dist_type>;
        dist_type lat1 = r*latInc- 0.5*std::numbers::pi_v<dist_type>;
        dist_type diagonalDist3DEUC = SBLoc<dist_type>::EUC3DDistFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
        for (size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            fillCacheCell({thisCtrLat, thisCtrLng}, diagonalDist3DEUC, colSize, ptLocPairs);
        }
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::calcSideLenFromAlpc() {
    dist_type surfaceArea = 4.0*std::numbers::pi_v<dist_type>*SBLoc<dist_type>::EARTH_RADIUS*SBLoc<dist_type>::EARTH_RADIUS;
    dist_type numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    totalLocSize = locData->size();
    BKDTSBSolver<KDTType, dist_type>::generateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc<dist_type>::deltaLatOnSameLngFromHavDist(sideLen));
    latIncInverse = 1.0/latInc;
    rowSize = std::ceil(std::numbers::pi_v<dist_type>/latInc);
    fillGridCache();
    //auto t = std::thread(&KDT<KDTType, dist_type>::clear, this->locKdt);
    //t.detach();
    this->locKdt.clear();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
returnNNLocFromCacheVariant(const Point<dist_type, 2>& geoPt, const BitCell& cell) const {
    if (cell.size() == 1) {
        return cell.getSingleLoc();
    } else if (cell.size() < MAX_CACHE_CELL_VEC_SIZE) {
        const auto p = SBLoc<dist_type>::geoPtToCart3DPt(geoPt);
        return custom_min_element(cell.getLocPairs(), cell.getLocPairs() + cell.size(),
                                  [&p](const auto& nh){ return p.template
                                      dist<Point<dist_type, 3>::DistType::EUCSQ>(nh.key);},
                                  std::less())->value;
    } else {
        return cell.getCacheTree()->kNNValue(SBLoc<dist_type>::geoPtToCart3DPt(geoPt), 1);
    }
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
findNearest(const Point<dist_type, 2>& geoSearchPt) const {
    return returnNNLocFromCacheVariant(geoSearchPt,
           gridCache[static_cast<size_t>((geoSearchPt[0]+0.5*std::numbers::pi_v<dist_type>)*latIncInverse)*colSize +
                     static_cast<size_t>((geoSearchPt[1]+std::numbers::pi_v<dist_type>)*lngIncInverse)]);
}



template class UnionUniLatLngBKDTGridSBSolver<KDTree, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTree, float>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float>;

template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float>;






