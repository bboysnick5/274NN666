//
//  UniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniLatLngBKDTGridSBSolver.hpp"
//#include <omp.h>

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
UniLatLngBKDTGridSBSolver(dist_type alpc, size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType, dist_type>(), AVE_LOC_PER_CELL(alpc),
  MAX_CACHE_CELL_VEC_SIZE(maxCacheCellVecSize) {}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::printSolverInfo() const {
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
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
fillCacheCell(dist_type thisCtrLng, dist_type thisCtrLat, dist_type thisDiff,
              std::vector<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::latLngToCart3DPt(thisCtrLat, thisCtrLng),
                                   thisDiff, std::back_inserter(ptLocPairs));
    size_t locsSize = ptLocPairs.size();
    this->totalNodeSize += locsSize;
    if (locsSize == 1) {
        this->gridCache.emplace_back(ptLocPairs[0].second);
        this->singleLocs++;
    } else if (locsSize < MAX_CACHE_CELL_VEC_SIZE) {
        this->gridCache.emplace_back(ptLocPairs);
        this->vecLocs++;
    } else {
        this->gridCache.emplace_back(std::in_place_type<KDT<KDTType, dist_type>>,
                                     ptLocPairs.begin(), ptLocPairs.end());
    }
    ptLocPairs.clear();
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2.0*M_PI/colSize + 2.0*M_PI/(colSize*colSize*0xFFFF);
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
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::calcSideLenFromAlpc() {
    dist_type surfaceArea = 4*M_PI*SBLoc<dist_type>::EARTH_RADIUS*SBLoc<dist_type>::EARTH_RADIUS;
    dist_type numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    totalLocSize = locData->size();
    BKDTSBSolver<KDTType, dist_type>::generateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc<dist_type>::latFromHavDist(sideLen, 0));
    latIncInverse = 1/latInc;
    rowSize = std::ceil(M_PI/(latInc - latInc*latInc/(M_PI*0xFFFF)));
    fillGridCache();
    this->locKdt.clear();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
returnNNLocFromCacheVariant(dist_type lat, dist_type lng,
                            const std::variant<std::vector<std::pair<Point<dist_type, 3>,
                            const SBLoc<dist_type>*>>, const SBLoc<dist_type>*, KDT<KDTType, dist_type>>& v) const {
    switch (v.index()) {
        case 0: {
            const auto p = SBLoc<dist_type>::latLngToCart3DPt(lat, lng);
            const auto &vec = std::get<0>(v);
            return std::min_element(vec.cbegin(), vec.cend(),
                [&](const auto& p1, const auto& p2){return Point<dist_type, 3>::template
                dist<Point<dist_type, 3>::DistType::EUCSQ>(p1.first, p)< Point<dist_type, 3>::template
                dist<Point<dist_type, 3>::DistType::EUCSQ>(p2.first, p);})->second;
        }
        case 1:
            return std::get<1>(v);
        default:
            return std::get<2>(v).kNNValue(SBLoc<dist_type>::latLngToCart3DPt(lat, lng), 1);
    }
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
findNearest(dist_type lat, dist_type lng) const {
    return returnNNLocFromCacheVariant(lat, lng, gridCache[static_cast<size_t>
    ((lat+0.5*M_PI)*latIncInverse)*colSize+ static_cast<size_t>((lng+M_PI)/lngInc)]);
}



template class UniLatLngBKDTGridSBSolver<KDTree, double>;
template class UniLatLngBKDTGridSBSolver<KDTree, float>;

template class UniLatLngBKDTGridSBSolver<KDTreeCusMem, double>;
template class UniLatLngBKDTGridSBSolver<KDTreeCusMem, float>;

template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double>;
template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float>;

template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float>;



/*

//
//  UniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniLatLngBKDTGridSBSolver.hpp"
//#include <omp.h>

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>

UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
UniLatLngBKDTGridSBSolver(dist_type alpc, size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType>(), AVE_LOC_PER_CELL(alpc),
MAX_CACHE_CELL_VEC_SIZE(maxCacheCellVecSize) {}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UniLatLngBKDTGridSBSolver<KDTType, dist_type>::~UniLatLngBKDTGridSBSolver() {
    std::for_each(gridCache.begin(), gridCache.end(), [](auto &v){
        if (auto tp = std::get_if<std::tuple<Point<dist_type, 3>*, const SBLoc<dist_type>**, size_t>>(&v)) {
            const auto &[ptArr, locPtrArr, size] = *tp;
            delete[] ptArr;
            delete[] locPtrArr;
        }
    });
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::printSolverInfo() const {
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
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
fillCacheCell(dist_type thisCtrLng, dist_type thisCtrLat, dist_type thisDiff,
              std::vector<std::pair<Point<dist_type, 3>, const SBLoc<dist_type>*>>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::latLngToCart3DPt(thisCtrLat, thisCtrLng),
                                   thisDiff, std::back_inserter(ptLocPairs));
    size_t locsSize = ptLocPairs.size();
    this->totalNodeSize += locsSize;
    if (locsSize == 1) {
        this->gridCache.emplace_back(ptLocPairs[0].second);
        this->singleLocs++;
    } else if (locsSize < MAX_CACHE_CELL_VEC_SIZE) {
        Point<dist_type, 3> *ptArr = new Point<dist_type, 3>[locsSize], *ptArrIt = ptArr;
        const SBLoc<dist_type>** locPtrArr = new const SBLoc<dist_type>*[locsSize], **locPtrArrIt = locPtrArr;
        for (auto ptLocPairsIt = ptLocPairs.cbegin(); ptLocPairsIt != ptLocPairs.cend();) {
            const auto &[pt, locPtr] = *ptLocPairsIt++;
            *ptArrIt++ = std::move(pt);
            *locPtrArrIt++ = locPtr;
        }
        this->gridCache.emplace_back(std::forward_as_tuple(ptArr, locPtrArr, locsSize));
        this->vecLocs++;
    } else {
        this->gridCache.emplace_back(std::in_place_type<const KDT<KDTType>>,
                                     ptLocPairs.begin(), ptLocPairs.end());
    }
    ptLocPairs.clear();
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2*M_PI/colSize + 2*M_PI/(colSize*colSize*0xFFFF);
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
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::calcSideLenFromAlpc() {
    dist_type surfaceArea = 4*M_PI*SBLoc<dist_type>::EARTH_RADIUS*SBLoc<dist_type>::EARTH_RADIUS;
    dist_type numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    totalLocSize = locData->size();
    BKDTSBSolver<KDTType>::generateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc<dist_type>::latFromHavDist(sideLen, 0));
    rowSize = std::ceil(M_PI/(latInc - latInc*latInc/(M_PI*0xFFFF)));
    fillGridCache();
    this->locKdt.clear();
}

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
returnNNLocFromCacheVariant(dist_type lng, dist_type lat,
                            const std::variant<std::tuple<Point<dist_type, 3>*,
                            const SBLoc<dist_type>**, size_t>, const SBLoc<dist_type>*, const KDT<KDTType>>& v) {
    switch (v.index()) {
        case 0: {
            const auto &[ptArr, locPtrArr, size] = std::get<0>(v);
            const auto searchPt = SBLoc<dist_type>::latLngToCart3DPt(lat, lng);
            return locPtrArr[std::min_element(ptArr, ptArr + size,
                                              [&](const auto& pt1, const auto& pt2){return Point<dist_type, 3>::template
                                                  dist<Point<dist_type, 3>::DistType::EUCSQ>(pt1, searchPt) < Point<dist_type, 3>::template
                                                  dist<Point<dist_type, 3>::DistType::EUCSQ>(pt2, searchPt);}) - ptArr];
        }
        case 1:
            return std::get<1>(v);
        default:
            return std::get<2>(v).kNNValue(SBLoc<dist_type>::latLngToCart3DPt(lat, lng), 1);
    }
}


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
findNearest(dist_type lng, dist_type lat) const {
    return returnNNLocFromCacheVariant(lng, lat, gridCache[static_cast<size_t>
                                                           ((lat+0.5*M_PI)/latInc)*colSize+ static_cast<size_t>((lng+M_PI)/lngInc)]);
}





template class UniLatLngBKDTGridSBSolver<KDTree>;
template class UniLatLngBKDTGridSBSolver<KDTreeCusMem>;
template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongest>;
template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec>;


*/
