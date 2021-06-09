//
//  UniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniLatLngBKDTGridSBSolver.hpp"
#include "Utility.hpp"
//#include <omp.h>

template <template <class DT, std::size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
UniLatLngBKDTGridSBSolver(dist_type alpc, std::size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType, dist_type>(), AVE_LOC_PER_CELL(alpc),
  MAX_CACHE_CELL_VEC_SIZE(maxCacheCellVecSize) {}


template <template <class DT, std::size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::PrintSolverInfo() const {
    std::cout << "Total cache locs: " << totalNodeSize
    << "\nAve tree size: " << totalNodeSize/grid_cache_.size()
    << "\nAve tree height: "
    << static_cast<std::size_t>(log2(totalNodeSize/grid_cache_.size() + 1)) + 1
    << "\nRatio of cache locs over actual num locs: "
    << totalNodeSize/totalLocSize
    << "\nTotal num of loc cells: " << grid_cache_.size()
    << "\nSingle loc cells: " << singleLocs
    << "\nVector loc cells: " << vecLocs
    << "\nKd-tree loc cells: " << grid_cache_.size() -singleLocs -vecLocs << "\n";
}


template <template <class DT, std::size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
FillCacheCell(dist_type thisCtrLng, dist_type thisCtrLat, dist_type diagonalDistSq3DEUC,
              std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt({thisCtrLat, thisCtrLng}),
                                   diagonalDistSq3DEUC, std::back_inserter(ptLocPairs));
    std::size_t locsSize = ptLocPairs.size();
    this->totalNodeSize += locsSize;
    if (locsSize == 1) {
        this->grid_cache_.emplace_back(ptLocPairs[0].value);
        this->singleLocs++;
    } else if (locsSize < MAX_CACHE_CELL_VEC_SIZE) {
        this->grid_cache_.emplace_back(ptLocPairs);
        this->vecLocs++;
    } else {
        this->grid_cache_.emplace_back(std::in_place_type<KDT<KDTType, dist_type>>,
                                     ptLocPairs.begin(), ptLocPairs.end());
    }
    ptLocPairs.clear();
}


template <template <class DT, std::size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::FillGridCache() {
    colSize = rowSize;
    lngInc = 2.0*def::kMathPi<dist_type>/colSize + 2.0*def::kMathPi<dist_type>/(colSize*65536);
    grid_cache_.reserve(this->locKdt.size()*1.2/AVE_LOC_PER_CELL);
    std::vector<typename KDT<KDTType, dist_type>::node_type> ptLocPairs;
    ptLocPairs.reserve(MAX_CACHE_CELL_VEC_SIZE);
    
    dist_type thisCtrLat = 0.5 * (latInc - def::kMathPi<dist_type>);
    for (std::size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc) {
        dist_type thisCtrLng = 0.5 * lngInc - def::kMathPi<dist_type>;
        dist_type lat1 = r*this->latInc - 0.5*def::kMathPi<dist_type>;
        dist_type diagonalDistSq3DEUC = SBLoc<dist_type>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + latInc, lngInc);
        for (std::size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            FillCacheCell(thisCtrLng, thisCtrLat, diagonalDistSq3DEUC, ptLocPairs);
        }
    }
}

template <template <class DT, std::size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::calcSideLenFromAlpc() {
    dist_type surfaceArea = 4*def::kMathPi<dist_type>*SBLoc<dist_type>::EARTH_RADIUS*SBLoc<dist_type>::EARTH_RADIUS;
    dist_type numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    sideLen = sqrt(surfaceArea/numCells);
}

template <template <class DT, std::size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
void UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    totalLocSize = locData->size();
    BKDTSBSolver<KDTType, dist_type>::GenerateKDT(locData);
    calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc<dist_type>::deltaLatOnSameLngFromHavDist(sideLen));
    latIncInverse = 1.0/latInc;
    rowSize = std::ceil(def::kMathPi<dist_type>/latInc);
    FillGridCache();
    this->locKdt.clear();
}

template <template <class DT, std::size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
ReturnNNLocFromCacheVariant(const Point<dist_type, 2>& geoSearchPt,
                            const std::variant<std::vector<typename KDT<KDTType, dist_type>::node_type>, const SBLoc<dist_type>*, KDT<KDTType, dist_type>>& v) const {
    switch (v.index()) {
        case 0: {
            const auto p = SBLoc<dist_type>::geoPtToCart3DPt(geoSearchPt);
            const auto &vec = std::get<0>(v);
            return Utility::MinElementGivenDistFunc(vec.cbegin(), vec.cend(),
                                               [&p](const auto& nh){return p.template
                                                    dist<Point<dist_type, 3>::DistType::EUCSQ>(nh.key);}, std::less()
                                               )->value;
        }
        case 1:
            return std::get<1>(v);
        default:
            return std::get<2>(v).kNNValue(SBLoc<dist_type>::geoPtToCart3DPt(geoSearchPt), 1);
    }
}


template <template <class DT, std::size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
const SBLoc<dist_type>* UniLatLngBKDTGridSBSolver<KDTType, dist_type>::
FindNearestLoc(const Point<dist_type, 2>& geoSearchPt) const {
    return ReturnNNLocFromCacheVariant(geoSearchPt, grid_cache_[static_cast<std::size_t>
    ((geoSearchPt[0]+0.5*def::kMathPi<dist_type>)*latIncInverse)*colSize+ static_cast<std::size_t>((geoSearchPt[1]+def::kMathPi<dist_type>)/lngInc)]);
}



template class UniLatLngBKDTGridSBSolver<KDTree, double>;
template class UniLatLngBKDTGridSBSolver<KDTree, float>;

template class UniLatLngBKDTGridSBSolver<KDTreeCusMem, double>;
template class UniLatLngBKDTGridSBSolver<KDTreeCusMem, float>;

template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double>;
template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float>;

template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float>;

