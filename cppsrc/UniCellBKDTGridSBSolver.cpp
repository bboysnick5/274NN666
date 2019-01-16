//
//  UniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniCellBKDTGridSBSolver.hpp"

//#include <omp.h>


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
UniCellBKDTGridSBSolver<KDTType>::
UniCellBKDTGridSBSolver(double alpc, size_t maxCacheCellVecSize)
: UniLatLngBKDTGridSBSolver<KDTType>(alpc, maxCacheCellVecSize) {}


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
void UniCellBKDTGridSBSolver<KDTType>::fillGridCache() {
    this->gridCache.reserve(this->locKdt.size()*1.2/this->AVE_LOC_PER_CELL);
    thisRowStartIdx.reserve(this->rowSize);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs;
    ptLocPairs.reserve(std::sqrt(this->locKdt.size()));
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    double thisLat = -0.5 * M_PI, thisCtrLat = thisLat + 0.5 * this->latInc;
    for (size_t r = 0, idx = 0; r < this->rowSize;
         ++r, thisCtrLat += this->latInc, thisLat += this->latInc) {
        size_t thisColSize = static_cast<size_t>(2*M_PI * SBLoc::EARTH_RADIUS *
                             cos(thisLat > 0 ? thisLat-this->latInc : thisLat)/
                             this->sideLen) + 2,
               thisEndIdx = idx + thisColSize;
        double thisLngInc = 2*M_PI/thisColSize +
                            2*M_PI/(thisColSize*thisColSize*0xFFFF),
               thisDiff = SBLoc::xyzDistFromLngLat(r*this->latInc - 0.5*M_PI,
                          (r+1)*this->latInc - 0.5*M_PI, thisLngInc),
               thisCtrLng = 0.5 * thisLngInc - M_PI;
        thisRowStartIdx.emplace_back(idx, thisLngInc);
        for (; idx < thisEndIdx; ++idx, thisCtrLng += thisLngInc) {
            UniLatLngBKDTGridSBSolver<KDTType>::fillCacheCell
            (thisCtrLng, thisCtrLat, thisDiff, ptLocPairs);
        }
    }
}

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
const SBLoc* UniCellBKDTGridSBSolver<KDTType>::
findNearest(double lng, double lat) const {
    const auto[startIdx, thisLngInc] = thisRowStartIdx[(lat+0.5*M_PI)/this->latInc];
    const auto &v = this->gridCache[startIdx +
                                    static_cast<size_t>((lng+M_PI)/thisLngInc)];
    switch (v.index()) {
        case 0: {
            const auto p = SBLoc::latLngToCart3DPt(lng, lat);
            const auto &vec = std::get<0>(v);
            return std::min_element
                (vec.begin(), vec.end(), [&](const auto& p1, const auto& p2) {
                    return Point<3>::template dist<Point<3>::DistType::EUCSQ>(p1.first, p)
                           < Point<3>::template dist<Point<3>::DistType::EUCSQ>(p2.first, p);
                })->second;
        }
        case 1:
            return std::get<1>(v);
        default:
            return std::get<2>(v).kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
    }
   // return std::holds_alternative<KDT<KDTType>>(v)
    //       ? std::get<0>(v).kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1)
    //       : std::get<const SBLoc*>(v);
}




template class UniCellBKDTGridSBSolver<KDTree>;
template class UniCellBKDTGridSBSolver<KDTreeCusMem>;
