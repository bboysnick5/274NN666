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


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::BitCell(uintptr_t otherPtr) : ptr(otherPtr) {}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::BitCell(BitCell&& rhs)
: ptr(rhs.ptr) {
    rhs.ptr = 0;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::BitCell(const BitCell& rhs) : ptr(rhs.ptr) {}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
typename UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell&
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::operator=(BitCell&& rhs) & noexcept {
    if (this != &rhs) {
        destruct();
        ptr = rhs.ptr;
        rhs.ptr = NULL;
    }
    return *this;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
typename UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell&
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::operator=(const BitCell& rhs) & noexcept {
    if (this != &rhs) {
        destruct();
        ptr = rhs.ptr;
    }
    return *this;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::
BitCell(const std::vector<typename KDT<KDTType, dist_type>::node_type> &bufVec,
        size_t maxCacheVecSize, const BitCell* left, const BitCell* up) {
    size_t size = bufVec.size();
    const uintptr_t MASK = ~(1ULL << 48);
    if (size == 1) {
        ptr = (reinterpret_cast<std::uintptr_t>(bufVec[0].value) & MASK) | (1ull << 48);
    } else if (size < maxCacheVecSize) {
        auto checkPrevCell = [size, &bufVec, this](const BitCell* cell) ->bool {
            if (uintptr_t cellPtrVal = cell->ptr;
                cellPtrVal && size == (cellPtrVal >> 48)
                && std::equal(bufVec.cbegin(), bufVec.cend(), reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(cellPtrVal << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT),
                              [](const auto &nh1, const auto &nh2){return nh1.value == nh2.value;})) {
                ptr = cellPtrVal & MASK_OUT_LEAST_SIG_BIT;
                return true;
            }
            return false;
        };
        if (left && checkPrevCell(left))
            return;
        if (up && checkPrevCell(up))
            return;
        auto *cacheLocs = static_cast<typename KDT<KDTType, dist_type>::node_type*>(::operator new(size*sizeof(typename KDT<KDTType, dist_type>::node_type), std::nothrow));
        std::uninitialized_move(bufVec.begin(), bufVec.end(), cacheLocs);
        ptr = (reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK) | (size << 48) | 1ull;
    } else {
        ptr = reinterpret_cast<std::uintptr_t>(new KDT<KDTType, dist_type>(bufVec.begin(), bufVec.end())) & MASK;
    }
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::destruct() {
    if (ptr != 0) {
        std::size_t size = ptr >> 48;
        if (size > 1 && (ptr & 1ull)) {
            ::operator delete(reinterpret_cast<void*>((static_cast<intptr_t>(ptr << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT));
        } else if (size == 0) {
            delete reinterpret_cast<KDT<KDTType, dist_type>*>(static_cast<intptr_t>(ptr << 16) >> 16);
        }
    }
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::~BitCell() {
    destruct();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
size_t UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::size() const {
    return ptr ? (ptr >> 48 ? ptr >> 48 : getCacheTree()->size()) : 0;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
size_t UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::rawSize() const {
    return ptr >> 48;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
bool UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::isUniqueVecLoc() const {
    return ptr & 1ull;
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
    std::size_t numSingleLocs = 0, numVecLocs = 0, numTreeLocs = 0, numUniqueVecLocs = 0;
    std::size_t totalTreeNodeSize = 0, totalVecLocSize = 0;
    std::for_each(gridCache.cbegin(), gridCache.cend(), [&](const BitCell& cell) mutable {
        if (auto cellSize = cell.size();
            cellSize == 1) {
            numSingleLocs++;
        } else if (cellSize < MAX_CACHE_CELL_VEC_SIZE) {
            numVecLocs++;
            totalVecLocSize += cellSize;
            if (cell.isUniqueVecLoc())
                numUniqueVecLocs++;
        } else {
            numTreeLocs++;
            totalTreeNodeSize += cellSize;
        }
    });
    std::size_t totalCacheLocs = numSingleLocs + totalVecLocSize + totalTreeNodeSize;
    
    std::cout << "Total cached locs: " << totalCacheLocs << std::endl
    << "Ratio of cache locs over actual num locs: " << static_cast<double>(totalCacheLocs)/totalLocSize << std::endl
    << "Total num of loc cells: " << gridCache.size() << std::endl
    << "Ave cell size: " << static_cast<double>(totalCacheLocs)/gridCache.size() << std::endl
    << "Single loc cells: " << numSingleLocs << std::endl
    << "Vector loc cells: " << numVecLocs << std::endl
    << "Unique Vector loc cells: " << numUniqueVecLocs << std::endl
    << "Ave vec loc size: " << static_cast<double>(totalVecLocSize)/numVecLocs << std::endl
    << "Kd-tree loc cells: " << numTreeLocs << std::endl
    << "Ave tree size : " << static_cast<double>(totalTreeNodeSize)/numTreeLocs << std::endl
    << "Ave tree height: " << log2(static_cast<double>(totalTreeNodeSize)/numTreeLocs + 1.0) + 1.0 << std::endl;
}

/*
 
 Serial Version

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

 */

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
fillCacheCell(typename std::vector<BitCell>::iterator gridCacheIt,
              const Point<dist_type, 2>& thisCtrGeoPt, dist_type diagonalDist3DEUC, size_t thisColSize,
              std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt(thisCtrGeoPt), diagonalDist3DEUC, std::back_inserter(ptLocPairs));
    *gridCacheIt = BitCell(ptLocPairs, MAX_CACHE_CELL_VEC_SIZE, &*(gridCacheIt-1), &*(gridCacheIt - thisColSize));
    ptLocPairs.clear();
}
 
template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2.0*std::numbers::pi_v<dist_type>/colSize + 2.0*std::numbers::pi_v<dist_type>/(colSize*65536);
    lngIncInverse = 1.0/lngInc;
    gridCache.resize(rowSize*colSize, NULL);
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    
    dist_type initCtrLat = 0.5*latInc - 0.5*std::numbers::pi_v<dist_type>;
    dist_type initCtrLng = 0.5*lngInc - std::numbers::pi_v<dist_type>;
    for (size_t r = 0; r < rowSize; ++r) {
        dist_type lat1 = r*latInc- 0.5*std::numbers::pi_v<dist_type>;
        dist_type diagonalDist3DEUC = SBLoc<dist_type>::EUC3DDistFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
        for (size_t c = 0; c < colSize; ++c) {
            this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt({
                initCtrLat + r*latInc, initCtrLng + c*lngInc}),
                diagonalDist3DEUC, std::back_inserter(ptLocPairs));
            std::size_t thisIdx = r*colSize + c;
            gridCache[thisIdx] = BitCell(ptLocPairs, MAX_CACHE_CELL_VEC_SIZE,
                                         &gridCache[thisIdx - 1], &gridCache[thisIdx - colSize]);
            ptLocPairs.clear();
        }
    }
}

/*
 
 Serial Verdsion
 
 
template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2.0*std::numbers::pi_v<dist_type>/colSize + 2.0*std::numbers::pi_v<dist_type>/(colSize*65536);
    lngIncInverse = 1.0/lngInc;
    gridCache.reserve(this->locKdt.size()*1.2/AVE_LOC_PER_CELL);
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    
    dist_type thisCtrLat = 0.5 * (latInc - std::numbers::pi_v<dist_type>);
    for (size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc) {
        dist_type thisCtrLng = 0.5 * lngInc - std::numbers::pi_v<dist_type>;
        dist_type lat1 = r*latInc- 0.5*std::numbers::pi_v<dist_type>;
        dist_type diagonalDist3DEUC = SBLoc<dist_type>::EUC3DDistFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
        for (size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            fillCacheCell({thisCtrLat, thisCtrLng}, diagonalDist3DEUC, colSize, ptLocPairs);
        }
    }
} */

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
        return custom_min_element(cell.getLocPairs(), cell.getLocPairs() + cell.size(),
                                  [&](const auto& nh) {return SBLoc<dist_type>::geoPtToCart3DPt(geoPt).template
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






