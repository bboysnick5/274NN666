//
//  UpgradeBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UpgradeBKDTGridSBSolver.hpp"
//#include <omp.h>

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
UpgradeBKDTGridSBSolver<Tree>::
UpgradeBKDTGridSBSolver(double alpc) : UniCellBKDTGridSBSolver<Tree>(alpc) {}

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
void UpgradeBKDTGridSBSolver<Tree>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2*M_PI/colSize + 2*M_PI/(colSize*colSize*0xFFFF);
    this->gridCache.reserve(this->locKdt.size()*1.2/this->AVE_LOC_PER_CELL);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs(this->locKdt.size());
    
    double thisCtrLat = 0.5 * (this->latInc - M_PI);
    for (size_t r = 0; r < rowSize; ++r, thisCtrLat += this->latInc) {
        double thisCtrLng = 0.5 * lngInc - M_PI;
        double thisDiff = SBLoc::xyzDistFromLngLat(r*this->latInc- 0.5*M_PI,
                          (r+1)*this->latInc-0.5*M_PI, lngInc);
        for (size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
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
void UpgradeBKDTGridSBSolver<Tree>::
build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    BKDTSBSolver<Tree>::generateKDT(locData);
    sideLen = UniCellBKDTGridSBSolver<Tree>::calcSideLenFromAlpc();
    this->latInc = std::fabs(SBLoc::latFromHavDist(sideLen, 0));
    rowSize = std::ceil(M_PI/(this->latInc - this->latInc*this->latInc/(M_PI*0xFFFF)));
    fillGridCache();
}

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
const SBLoc* UpgradeBKDTGridSBSolver<Tree>::
findNearest(double lng, double lat) const {
    const auto&[cacheTree, singleLoc] =
        this->gridCache[static_cast<size_t>((lat+0.5*M_PI)/this->latInc)*colSize
                        + static_cast<size_t>((lng+M_PI)/lngInc)];
    return singleLoc ? singleLoc :
           cacheTree.kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}





template class UpgradeBKDTGridSBSolver<KDTree>;
template class UpgradeBKDTGridSBSolver<KDTreeCusMem>;


