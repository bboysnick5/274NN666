//
//  UniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniCellBKDTGridSBSolver.hpp"

//#include <omp.h>




template <template <size_t, typename, typename Point<3>::DistType> class Tree>
UniCellBKDTGridSBSolver<Tree>:: UniCellBKDTGridSBSolver(double alpc)
: UniLatLngBKDTGridSBSolver<Tree>(alpc) {}



template <template <size_t, typename, typename Point<3>::DistType> class Tree>
void UniCellBKDTGridSBSolver<Tree>::fillGridCache() {
    this->gridCache.reserve(this->locKdt.size()*1.2/this->AVE_LOC_PER_CELL);
    thisRowStartIdx.reserve(this->rowSize);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs(this->locKdt.size());
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    double thisLat = -0.5 * M_PI, thisCtrLat = thisLat + 0.5 * this->latInc;
    for (size_t r = 0, idx = 0; r < this->rowSize;
         ++r, thisCtrLat += this->latInc, thisLat += this->latInc) {
        size_t thisColSize = static_cast<size_t>(2*M_PI * SBLoc::EARTH_RADIUS *
                             cos(thisLat > 0 ? thisLat-this->latInc
                                 : thisLat)/this->sideLen) + 2,
               thisEndIdx = idx + thisColSize;
        double thisLngInc = 2*M_PI/thisColSize +
                            2*M_PI/(thisColSize*thisColSize*0xFFFF),
               thisDiff = SBLoc::xyzDistFromLngLat(r*this->latInc - 0.5*M_PI,
                          (r+1)*this->latInc - 0.5*M_PI, thisLngInc),
               thisCtrLng = 0.5 * thisLngInc - M_PI;
        thisRowStartIdx.emplace_back(idx, thisLngInc);
        for (; idx < thisEndIdx; ++idx, thisCtrLng += thisLngInc) {
            auto locsEnd = this->locKdt.rangeDiffKNNPairs
                           (SBLoc::latLngToCart3DPt(thisCtrLng, thisCtrLat),
                           thisDiff, ptLocPairs.begin());
            size_t locsSize = locsEnd - ptLocPairs.begin();
            switch (locsSize) {
                case 1:
                    this->gridCache.emplace_back(std::piecewise_construct,
                                                 std::forward_as_tuple(),
                                                 std::forward_as_tuple
                                                 (ptLocPairs[0].second));
                    this->singleLocs++;
                    break;
                default:
                    this->gridCache.emplace_back(std::piecewise_construct,
                                                 std::forward_as_tuple
                                                 (ptLocPairs.begin(), locsEnd),
                                                 std::forward_as_tuple(nullptr));
            }
            this->totalNodeSize += locsSize;
        }
    }
}

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
void UniCellBKDTGridSBSolver<Tree>::
build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    BKDTSBSolver<Tree>::generateKDT(locData);
    this->sideLen = UniLatLngBKDTGridSBSolver<Tree>::calcSideLenFromAlpc();
    this->latInc = std::fabs(SBLoc::latFromHavDist(this->sideLen, 0));
    this->rowSize = std::ceil(M_PI/(this->latInc -
                                    this->latInc*this->latInc/(M_PI*0xFFFF)));
    fillGridCache();
}

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
const SBLoc* UniCellBKDTGridSBSolver<Tree>::
findNearest(double lng, double lat) const {
    const auto[startIdx, thisLngInc] = thisRowStartIdx[(lat+0.5*M_PI)/this->latInc];
    const auto&[cacheTree, singleLoc] =
        this->gridCache[startIdx + static_cast<size_t>((lng+M_PI)/thisLngInc)];
    return singleLoc ? singleLoc :
           cacheTree.kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}





template class UniCellBKDTGridSBSolver<KDTree>;
template class UniCellBKDTGridSBSolver<KDTreeCusMem>;
