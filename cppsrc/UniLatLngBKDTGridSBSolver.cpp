//
//  UniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniLatLngBKDTGridSBSolver.hpp"
//#include <omp.h>

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
UniLatLngBKDTGridSBSolver<Tree>::
UniLatLngBKDTGridSBSolver(double alpc) : BKDTSBSolver<Tree>(),
                                         AVE_LOC_PER_CELL(alpc) {}


template <template <size_t, typename, typename Point<3>::DistType> class Tree>
void UniLatLngBKDTGridSBSolver<Tree>::printSolverInfo() const {
    std::cout << "Total tree nodes: " << totalNodeSize
    << "\nAve tree size: " << totalNodeSize/gridCache.size()
    << "\nAve tree height: "
    << static_cast<size_t>(log2(totalNodeSize/gridCache.size() + 1)) + 1
    << "\nRatio of tree nodes over num locs: "
    << totalNodeSize/this->locKdt.size()
    << "\nSingle loc cells: " << singleLocs << "\nMulti-loc cells: "
    << gridCache.size() - singleLocs << "\n";
}

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
void UniLatLngBKDTGridSBSolver<Tree>::fillGridCache() {
    colSize = rowSize;
    lngInc = 2*M_PI/colSize + 2*M_PI/(colSize*colSize*0xFFFF);
    gridCache.reserve(this->locKdt.size()*1.2/AVE_LOC_PER_CELL);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs(this->locKdt.size());
    
    double thisCtrLat = 0.5 * (latInc - M_PI);
    for (size_t r = 0; r < rowSize; ++r, thisCtrLat += latInc) {
        double thisCtrLng = 0.5 * lngInc - M_PI;
        double thisDiff = SBLoc::xyzDistFromLngLat(r*latInc- 0.5*M_PI,
                          (r+1)*latInc-0.5*M_PI, lngInc);
        for (size_t c = 0; c < colSize; ++c, thisCtrLng += lngInc) {
            auto locsEnd = this->locKdt.rangeDiffKNNPairs
                           (SBLoc::latLngToCart3DPt(thisCtrLng, thisCtrLat),
                            thisDiff, ptLocPairs.begin());
            size_t locsSize = locsEnd - ptLocPairs.begin();
            switch (locsSize) {
                case 1:
                    gridCache.emplace_back(std::piecewise_construct,
                                                 std::forward_as_tuple(),
                                                 std::forward_as_tuple
                                                 (ptLocPairs[0].second));
                    singleLocs++;
                    break;
                default:
                    gridCache.emplace_back(std::piecewise_construct,
                                                 std::forward_as_tuple
                                                 (ptLocPairs.begin(), locsEnd),
                                                 std::forward_as_tuple(nullptr));
            }
            totalNodeSize += locsSize;
        }
    }
}

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
double UniLatLngBKDTGridSBSolver<Tree>::calcSideLenFromAlpc() {
    double surfaceArea = 4*M_PI*SBLoc::EARTH_RADIUS*SBLoc::EARTH_RADIUS;
    double numCells = this->locKdt.size()/AVE_LOC_PER_CELL;
    return sqrt(surfaceArea/numCells);
}

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
void UniLatLngBKDTGridSBSolver<Tree>::
build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    BKDTSBSolver<Tree>::generateKDT(locData);
    sideLen = calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc::latFromHavDist(sideLen, 0));
    rowSize = std::ceil(M_PI/(latInc -
                              latInc*latInc/(M_PI*0xFFFF)));
    fillGridCache();
}

template <template <size_t, typename, typename Point<3>::DistType> class Tree>
const SBLoc* UniLatLngBKDTGridSBSolver<Tree>::
findNearest(double lng, double lat) const {
    const auto&[cacheTree, singleLoc] =
        gridCache[static_cast<size_t>((lat+0.5*M_PI)/latInc)*colSize
                        + static_cast<size_t>((lng+M_PI)/lngInc)];
    return singleLoc ? singleLoc :
           cacheTree.kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}





template class UniLatLngBKDTGridSBSolver<KDTree>;
template class UniLatLngBKDTGridSBSolver<KDTreeCusMem>;


