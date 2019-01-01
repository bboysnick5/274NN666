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
UniCellBKDTGridSBSolver(double alpc) : GridSBSolver(alpc) {}

void UniCellBKDTGridSBSolver::
fillGridCache(const KDTree<3, const SBLoc*, DistType::EUC>& sbKdt) {
    gridTreeCache.reserve(numLocs/AVE_LOC_PER_CELL);
    gridSingleCache.reserve(numLocs/AVE_LOC_PER_CELL);
    std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>> ptLocPairs(numLocs);
    thisRowStartIdx = std::vector<std::pair<size_t, double>>(rowSize);
    thisRowStartIdx[0].first = 0;
    size_t totalTreeSize = 0, singleLocs = 0;
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    for (size_t r = 0, c= 0, idx = 0, thisColSize = 0; r < rowSize; ++r) {
        thisColSize = std::max(static_cast<size_t>(cos(((r+1>rowSize/2?r:r+1) *
                                                        latInc-90.0)*M_PI/180.0)*
                                                   SBLoc::EARTH_RADIUS *
                                                   2*M_PI/sideLen),
                               1lu);
        if (r + 1 < rowSize)
            thisRowStartIdx[r+1].first = thisRowStartIdx[r].first + thisColSize;
        gridTreeCache.resize(gridTreeCache.size() + thisColSize);
        gridSingleCache.resize(gridSingleCache.size() + thisColSize);
        double thisLngInc = 360.0/thisColSize +
                            360.0/(thisColSize*thisColSize*0xFFFF);
        //std::numeric_limits<double>::epsilon()*0xFFFFFFF;
        thisRowStartIdx[r].second = thisLngInc;
        for (c = 0, idx = thisRowStartIdx[r].first; c < thisColSize; ++c, ++idx) {
            auto locsEnd = sbKdt.rangeDiffKNNPairs(
                           SBLoc::latLngToCart3DPt((c+0.5) * thisLngInc - 180.0,
                                                    (r+0.5) * latInc - 90.0),
                           SBLoc::xyzDistFromLngLat(r*latInc-90.0,
                                                    (r+1)*latInc-90.0,
                                                    thisLngInc),
                           ptLocPairs.begin());
            size_t locsSize = locsEnd - ptLocPairs.begin();
            if (locsSize > 1) {
                gridTreeCache[idx] = KDTree<3, const SBLoc*, DistType::EUC>
                (ptLocPairs.begin(), locsEnd);
            } else {
                gridSingleCache[idx] = ptLocPairs[0].second;
                singleLocs++;
            }
            totalTreeSize += locsSize;
        }
    }
    size_t multiLocs = gridSingleCache.size() - singleLocs;
    std::cout << "ave tree size: " << totalTreeSize/(gridTreeCache.size())
    << "\nSingle loc cells: " << singleLocs
    << "\nMulti-loc cells:" << multiLocs << std::endl;
}

double UniCellBKDTGridSBSolver::calcSideLenFromAlpc() {
    double surfaceArea = 4*M_PI*SBLoc::EARTH_RADIUS*SBLoc::EARTH_RADIUS;
    double numCells = numLocs/AVE_LOC_PER_CELL;
    return sqrt(surfaceArea/numCells);
}

void UniCellBKDTGridSBSolver::build() {
    std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>> kdtData(sbData->size());
    std::transform(sbData->begin(), sbData->end(), kdtData.begin(),
                   [&](const SBLoc& l)->std::pair<Point<3, DistType::EUC>, const SBLoc*>{return
                       {l.locToCart3DPt(), &l};});
    auto sbKdt = KDTree<3, const SBLoc*, DistType::EUC>(kdtData.begin(), kdtData.end());
    numLocs = sbKdt.size();
    sideLen = calcSideLenFromAlpc();
    latInc = std::fabs(SBLoc::latFromHavDist(sideLen, 0));
    rowSize = std::ceil(180.0/(latInc -
                               std::numeric_limits<double>::epsilon()*0xFFFFF));
    fillGridCache(sbKdt);
}

const SBLoc* UniCellBKDTGridSBSolver::
findNearest(double lng, double lat) const {
    auto[startIdx, thisLngInc] = thisRowStartIdx[(lat+90.0)/latInc];
    size_t idx = startIdx + static_cast<size_t>((lng+180.0)/thisLngInc);
    auto singleLoc = gridSingleCache[idx];
    return singleLoc ? singleLoc :
           gridTreeCache[idx].kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}


