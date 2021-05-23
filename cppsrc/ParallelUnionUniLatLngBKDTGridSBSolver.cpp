//
//  ParallelUnionUniLatBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 5/23/21.
//  Copyright © 2021 Yunlong Liu. All rights reserved.
//

#include "ParallelUnionUniLatLngBKDTGridSBSolver.hpp"
#include "Utility.hpp"
#include <omp.h>
#include <memory>
#include <thread>


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::BitCell(uintptr_t otherPtr) : ptr(otherPtr) {}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::BitCell(BitCell&& rhs)
: ptr(rhs.ptr.load(std::memory_order_relaxed)) {
    rhs.ptr.store(0, std::memory_order_relaxed);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::BitCell(const BitCell& rhs)
: ptr(rhs.ptr.load(std::memory_order_relaxed)) {}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
typename ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell&
ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::operator=(BitCell&& rhs) & noexcept {
    if (this != &rhs) {
        if (uintptr_t ptrValue = ptr.load(std::memory_order_relaxed);
            ptrValue != 0) {
            auto size = ptrValue >> 48;
            if (size > 1 && (ptrValue & 1ull) == 1ull) {
            //delete[] getLocPairs();
                ::operator delete(reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(ptrValue << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT));
            } else if (size == 0) {
                delete reinterpret_cast<const KDT<KDTType, dist_type>*>(static_cast<intptr_t>(ptrValue << 16) >> 16);
            }
        }
        ptr.store(rhs.ptr.load(std::memory_order_relaxed), std::memory_order_relaxed);
        rhs.ptr.store(0, std::memory_order_relaxed);
    }
    return *this;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
typename ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell&
ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::operator=(const BitCell& rhs) & noexcept {
    if (this != &rhs) {
        if (uintptr_t ptrValue = ptr.load(std::memory_order_relaxed);
            ptrValue != 0) {
            auto size = ptrValue >> 48;
            if (size > 1 && (ptrValue & 1ull) == 1ull) {
                ::operator delete(reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(ptrValue << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT));
            } else if (size == 0) {
                delete reinterpret_cast<const KDT<KDTType, dist_type>*>(static_cast<intptr_t>(ptrValue << 16) >> 16);
            }
        }
        ptr.store(rhs.ptr.load(std::memory_order_relaxed), std::memory_order_relaxed);
    }
    return *this;
}



template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::
BitCell(std::vector<typename KDT<KDTType, dist_type>::node_type> &bufVec,
        size_t maxCacheVecSize, const BitCell *left, const BitCell *up) {
    if (size_t size = bufVec.size();
        size == 1) {
        ptr.store((reinterpret_cast<std::uintptr_t>(bufVec[0].value) & MASK_OUT_16TH_BIT) | (1ull << 48), std::memory_order_relaxed);
    } else if (size < maxCacheVecSize) {
        auto checkPrevCell = [size, &bufVec, this](const BitCell* cell) ->bool {
            if (uintptr_t cellPtrVal = cell->ptr.load(std::memory_order_relaxed);
                cellPtrVal && size == (cellPtrVal >> 48)
                && std::equal(bufVec.cbegin(), bufVec.cend(), reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(cellPtrVal << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT),
                              [](const auto &nh1, const auto &nh2){return nh1.value == nh2.value;})) {
                ptr.store(cellPtrVal & MASK_OUT_LEAST_SIG_BIT, std::memory_order_relaxed);
                return true;
            }
            return false;
        };
        if (left && checkPrevCell(left))
            return;
        if (up && checkPrevCell(up))
            return;
        auto *cacheLocs = static_cast<typename KDT<KDTType, dist_type>::node_type*>(
                          ::operator new(size*sizeof(typename KDT<KDTType, dist_type>::node_type), std::nothrow));
        std::uninitialized_move(bufVec.begin(), bufVec.end(), cacheLocs);
        ptr.store((reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK_OUT_16TH_BIT) | (size << 48) | 1ull,
                  std::memory_order_relaxed);
    } else {
        ptr.store(reinterpret_cast<std::uintptr_t>(new KDT<KDTType, dist_type>(bufVec.begin(), bufVec.end())) & MASK_OUT_16TH_BIT,
                  std::memory_order_relaxed);
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::~BitCell() {
    uintptr_t ptrValue = ptr.load(std::memory_order_relaxed);
    if (std::size_t size = ptrValue >> 48;
        size > 1 && (ptrValue & 1ull)) {
        //delete[] getLocPairs();
        ::operator delete(reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(ptrValue << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT));
    } else if (size == 0 && ptrValue != 0) {
        delete reinterpret_cast<const KDT<KDTType, dist_type>*>(static_cast<intptr_t>(ptrValue << 16) >> 16);
    }
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
size_t ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::size() const {
    //auto size = ptr >> 48;
    //return size != 0 ? size : getCacheTree()->size();
    uintptr_t ptrValue = ptr.load(std::memory_order_relaxed);
    return ptrValue == 0 ? 0 : ((ptrValue >> 48) != 0 ? ptrValue >> 48 : reinterpret_cast<const KDT<KDTType, dist_type>*>(static_cast<intptr_t>(ptrValue << 16) >> 16)->size());
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
bool ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::isUniqueVecLoc() const {
    return ptr & 1ull;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
size_t ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::rawSize() const {
    return ptr.load(std::memory_order_relaxed) >> 48;
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::getSingleLoc() const {
    return reinterpret_cast<const SBLoc<dist_type>*>(static_cast<intptr_t>(ptr.load(std::memory_order_relaxed) << 16) >> 16);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const typename KDT<KDTType, dist_type>::node_type*
ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::getLocPairs() const {
    return reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(ptr.load(std::memory_order_relaxed) << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const KDT<KDTType, dist_type>* ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::BitCell::getCacheTree() const {
    return reinterpret_cast<const KDT<KDTType, dist_type>*>(static_cast<intptr_t>(ptr.load(std::memory_order_relaxed) << 16) >> 16);
}



template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
ParallelUnionUniLatLngBKDTGridSBSolver(dist_type alpc, size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType, dist_type>(), AVE_LOC_PER_CELL(alpc),
MAX_CACHE_CELL_VEC_SIZE(maxCacheCellVecSize) {}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::printSolverInfo() const {
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
void ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
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
void ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
fillCacheCell(typename std::vector<BitCell>::iterator gridCacheIt,
              const Point<dist_type, 2>& thisCtrGeoPt, dist_type diagonalDist3DEUC, size_t thisColSize,
              std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt(thisCtrGeoPt), diagonalDist3DEUC, std::back_inserter(ptLocPairs));
    *gridCacheIt = BitCell(ptLocPairs, MAX_CACHE_CELL_VEC_SIZE, &*(gridCacheIt-1), &*(gridCacheIt - thisColSize));
    ptLocPairs.clear();
}
 
template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2.0*std::numbers::pi_v<dist_type>/colSize + 2.0*std::numbers::pi_v<dist_type>/(colSize*65536);
    lngIncInverse = 1.0/lngInc;
    gridCache.resize(rowSize*colSize, 0);
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    
    dist_type initCtrLat = 0.5*latInc - 0.5*std::numbers::pi_v<dist_type>;
    dist_type initCtrLng = 0.5*lngInc - std::numbers::pi_v<dist_type>;
#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
shared(this->gridCache, lngInc, latInc, initCtrLat, initCtrLng) \
firstprivate(ptLocPairs) default(none) schedule(static, 1) ordered collapse(2)
    for (size_t r = 0; r < rowSize; ++r) {
        dist_type lat1 = r*latInc- 0.5*std::numbers::pi_v<dist_type>;
        dist_type diagonalDist3DEUC = SBLoc<dist_type>::EUC3DDistFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
        for (size_t c = 0; c < colSize; ++c) {
            this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt({
                initCtrLat + r*latInc, initCtrLng + c*lngInc}),
                diagonalDist3DEUC, std::back_inserter(ptLocPairs));
            std::size_t thisIdx = r*colSize + c;
            gridCache[thisIdx] = BitCell(ptLocPairs, MAX_CACHE_CELL_VEC_SIZE,
                                         thisIdx > 0 ? &gridCache[thisIdx - 1] : nullptr,
                                         thisIdx >= colSize ? &gridCache[thisIdx - colSize] : nullptr);
            ptLocPairs.clear();
        }
    }
}

/*
 
 Serial Verdsion
 
 
template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
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
void ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::calcSideLenFromAlpc() {
    dist_type surfaceArea = 4.0*std::numbers::pi_v<dist_type>*SBLoc<dist_type>::EARTH_RADIUS*SBLoc<dist_type>::EARTH_RADIUS;
    dist_type numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    totalLocSize = locData->size();
    BKDTSBSolver<KDTType, dist_type>::generateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc<dist_type>::deltaLatOnSameLngFromHavDist(sideLen));
    latIncInverse = 1.0/latInc;
    rowSize = std::ceil(std::numbers::pi_v<dist_type>/latInc);
    fillGridCache();
    this->locKdt = {};
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
returnNNLocFromCacheVariant(const Point<dist_type, 2>& geoPt, const BitCell& cell) const {
    if (auto size = cell.size();
        size == 1) {
        return cell.getSingleLoc();
    } else if (size < MAX_CACHE_CELL_VEC_SIZE) {
        auto locPairsPt = cell.getLocPairs();
        return custom_min_element(locPairsPt, locPairsPt + size,
                                  [&](const auto& nh) {return SBLoc<dist_type>::geoPtToCart3DPt(geoPt).template
                                      dist<Point<dist_type, 3>::DistType::EUCSQ>(nh.key);},
                                  std::less())->value;
    } else {
        return cell.getCacheTree()->kNNValue(SBLoc<dist_type>::geoPtToCart3DPt(geoPt), 1);
    }
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* ParallelUnionUniLatLngBKDTGridSBSolver<KDTType, dist_type>::
findNearest(const Point<dist_type, 2>& geoSearchPt) const {
    return returnNNLocFromCacheVariant(geoSearchPt,
           gridCache[static_cast<size_t>((geoSearchPt[0]+0.5*std::numbers::pi_v<dist_type>)*latIncInverse)*colSize +
                     static_cast<size_t>((geoSearchPt[1]+std::numbers::pi_v<dist_type>)*lngIncInverse)]);
}



template class ParallelUnionUniLatLngBKDTGridSBSolver<KDTree, double>;
template class ParallelUnionUniLatLngBKDTGridSBSolver<KDTree, float>;

template class ParallelUnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, double>;
template class ParallelUnionUniLatLngBKDTGridSBSolver<KDTreeCusMem, float>;

template class ParallelUnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double>;
template class ParallelUnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float>;

template class ParallelUnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class ParallelUnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float>;



