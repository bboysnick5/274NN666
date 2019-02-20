//
//  UnionUniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 2/19/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#include "UnionUniLatLngBKDTGridSBSolver.hpp"
#include <memory>

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
UnionUniLatLngBKDTGridSBSolver<KDTType>::Cell::
Cell(size_t maxCacheVecSize, const std::vector<std::pair<Point<3>,
     const SBLoc*>>& bufVec) : _size(bufVec.size()) {
    if (_size == 1) {
        cacheLoc = bufVec[0].second;
    } else if (_size < maxCacheVecSize) {
        //cacheLocs = static_cast<std::pair<Point<3>, const SBLoc*>*>(::operator new[](_size, std::nothrow));
        //std::uninitialized_move(bufVec.cbegin(), bufVec.cend(), cacheLocs);
        cacheLocs = new std::pair<Point<3>, const SBLoc*>[_size];
        std::copy(bufVec.cbegin(), bufVec.cend(), cacheLocs);
    } else {
        _size = 0;
        cacheTree = new KDT<KDTType>(bufVec.begin(), bufVec.end());
    }
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
UnionUniLatLngBKDTGridSBSolver<KDTType>::Cell::~Cell() {
    if (_size > 1) {
        delete cacheLocs;
    } else if (_size == 0) {
        delete cacheTree;
    }
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
size_t UnionUniLatLngBKDTGridSBSolver<KDTType>::Cell::size() const {
    return _size != 0 ? _size : cacheTree->size();
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
UnionUniLatLngBKDTGridSBSolver<KDTType>::
UnionUniLatLngBKDTGridSBSolver(double alpc, size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType>(), AVE_LOC_PER_CELL(alpc),
MAX_CACHE_CELL_VEC_SIZE(maxCacheCellVecSize) {}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UnionUniLatLngBKDTGridSBSolver<KDTType>::printSolverInfo() const {
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

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UnionUniLatLngBKDTGridSBSolver<KDTType>::
fillCacheCell(double thisCtrLng, double thisCtrLat, double thisDiff,
              std::vector<std::pair<Point<3>, const SBLoc*>>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc::latLngToCart3DPt(thisCtrLng, thisCtrLat),
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

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UnionUniLatLngBKDTGridSBSolver<KDTType>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2*M_PI/colSize + 2*M_PI/(colSize*colSize*0xFFFF);
    gridCache.reserve(this->locKdt.size()*1.2/AVE_LOC_PER_CELL);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    
    double thisCtrLat = 0.5 * (latInc - M_PI);
    for (size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc) {
        double thisCtrLng = 0.5 * lngInc - M_PI;
        double thisDiff = SBLoc::xyzDistFromLngLat(r*latInc- 0.5*M_PI,
                                                   (r+1)*latInc-0.5*M_PI, lngInc);
        for (size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            fillCacheCell(thisCtrLng, thisCtrLat, thisDiff, ptLocPairs);
        }
    }
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UnionUniLatLngBKDTGridSBSolver<KDTType>::calcSideLenFromAlpc() {
    double surfaceArea = 4*M_PI*SBLoc::EARTH_RADIUS*SBLoc::EARTH_RADIUS;
    double numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UnionUniLatLngBKDTGridSBSolver<KDTType>::
build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    totalLocSize = locData->size();
    BKDTSBSolver<KDTType>::generateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc::latFromHavDist(sideLen, 0));
    rowSize = std::ceil(M_PI/(latInc - latInc*latInc/(M_PI*0xFFFF)));
    fillGridCache();
    this->locKdt.clear();
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
const SBLoc* UnionUniLatLngBKDTGridSBSolver<KDTType>::
returnNNLocFromCacheVariant(double lng, double lat, const Cell& cell) const {
    if (cell.size() == 1) {
        return cell.cacheLoc;
    } else if (cell.size() < MAX_CACHE_CELL_VEC_SIZE) {
        const auto p = SBLoc::latLngToCart3DPt(lng, lat);
        return std::min_element(cell.cacheLocs, cell.cacheLocs + cell.size(),
                                [&](const auto& p1, const auto& p2){return Point<3>::template
                                    dist<Point<3>::DistType::EUCSQ>(p1.first, p)< Point<3>::template
                                    dist<Point<3>::DistType::EUCSQ>(p2.first, p);})->second;
    } else {
        return cell.cacheTree->kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
    }
}


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
const SBLoc* UnionUniLatLngBKDTGridSBSolver<KDTType>::
findNearest(double lng, double lat) const {
    return returnNNLocFromCacheVariant(lng, lat, gridCache[static_cast<size_t>
                                                           ((lat+0.5*M_PI)/latInc)*colSize+ static_cast<size_t>((lng+M_PI)/lngInc)]);
}





template class UnionUniLatLngBKDTGridSBSolver<KDTree>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeCusMem>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest>;
template class UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec>;



