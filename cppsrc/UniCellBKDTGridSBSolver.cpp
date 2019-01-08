//
//  UniCellBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniCellBKDTGridSBSolver.hpp"

//#include <omp.h>


UniCellBKDTGridSBSolver::
UniCellBKDTGridSBSolver(double alpc) : BKDTSBSolver(), AVE_LOC_PER_CELL(alpc) {}

void UniCellBKDTGridSBSolver::fillGridCache(size_t rowSize, double sideLen) {
    gridCache.reserve(locKdt.size()*1.2/AVE_LOC_PER_CELL);
    thisRowStartIdx.reserve(rowSize);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs(locKdt.size());
    size_t totalTreeSize = 0, singleLocs = 0;
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency())\
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    double thisLat = -0.5*M_PI, thisCtrLat = thisLat + 0.5*latInc;
    for (size_t r = 0, idx = 0; r < rowSize;
         ++r, thisCtrLat += latInc, thisLat += latInc) {
        size_t thisColSize = static_cast<size_t>(2*M_PI * SBLoc::EARTH_RADIUS *
                             cos(thisLat>0?thisLat-latInc:thisLat)/sideLen) + 1,
               thisEndIdx = idx + thisColSize;
        double thisLngInc = 2*M_PI/thisColSize +
                            2*M_PI/(thisColSize*thisColSize*0xFF),
               thisDiff = SBLoc::xyzDistFromLngLat(r*latInc - 0.5*M_PI,
                          (r+1) * latInc - 0.5*M_PI, thisLngInc),
               thisCtrLng = 0.5 * thisLngInc - M_PI;
        thisRowStartIdx.emplace_back(idx, thisLngInc);
        for (; idx < thisEndIdx; ++idx, thisCtrLng += thisLngInc) {
            auto locsEnd = locKdt.rangeDiffKNNPairs
                           (SBLoc::latLngToCart3DPt(thisCtrLng, thisCtrLat),
                           thisDiff, ptLocPairs.begin());
            size_t locsSize = locsEnd - ptLocPairs.begin();
            if (locsSize > 1) {
                gridCache.emplace_back(KDTree<3, const SBLoc*,
                Point<3>::DistType::EUC>(ptLocPairs.begin(), locsEnd), nullptr);
            } else {
                gridCache.emplace_back(KDTree<3, const SBLoc*,
                                       Point<3>::DistType::EUC>(),
                                       ptLocPairs[0].second);
                singleLocs++;
            }
            totalTreeSize += locsSize;
        }
    }
    std::cout << "Total tree nodes: " << totalTreeSize
              << "\nAve tree size: " << totalTreeSize/gridCache.size()
              << "\nAve tree height: "
              << static_cast<size_t>(log2(totalTreeSize/gridCache.size()+1)) + 1
              << "\nRatio of tree nodes over num locs: "
              << totalTreeSize/locKdt.size()
              << "\nSingle loc cells: " << singleLocs << "\nMulti-loc cells:"
              << gridCache.size() - singleLocs << std::endl;
}

double UniCellBKDTGridSBSolver::calcSideLenFromAlpc() {
    double surfaceArea = 4*M_PI*SBLoc::EARTH_RADIUS*SBLoc::EARTH_RADIUS;
    double numCells = locKdt.size()/AVE_LOC_PER_CELL;
    return sqrt(surfaceArea/numCells);
}

void UniCellBKDTGridSBSolver::
build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    BKDTSBSolver::generateKDT(locData);
    double sideLen = calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc::latFromHavDist(sideLen, 0));
    size_t rowSize = std::ceil(M_PI/(latInc - latInc*latInc/(M_PI*0xFF)));
    fillGridCache(rowSize, sideLen);
}

const SBLoc* UniCellBKDTGridSBSolver::
findNearest(double lng, double lat) const {
    auto[startIdx, thisLngInc] = thisRowStartIdx[(lat+0.5*M_PI)/latInc];
    size_t idx = startIdx + static_cast<size_t>((lng+M_PI)/thisLngInc);
    return gridCache[idx].second ? gridCache[idx].second :
           gridCache[idx].first.kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}





