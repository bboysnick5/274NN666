//
//  UpgradeBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UpgradeBKDTGridSBSolver.hpp"
//#include <omp.h>


UpgradeBKDTGridSBSolver::
UpgradeBKDTGridSBSolver(double alpc) : GridSBSolver(alpc) {}


void UpgradeBKDTGridSBSolver::fillGridCache() {
    gridTreeCache = std::vector<KDTree<3, const SBLoc*, DistType::EUC>>
    (rowSize*colSize);
    gridSingleCache = std::vector<const SBLoc*>(rowSize*colSize);
    std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>> ptLocPairs(numLocs);
    size_t totalTreeSize = 0, singleLocs = 0;
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
    //default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
    //reduction(+:totalTreeSize, singleLocs) collapse(2)
    for (size_t r = 0; r < rowSize; ++r) {
        for (size_t c = 0; c < colSize; ++c) {
            size_t idx = r*colSize+c;
            auto locsEnd = sbKdt.rangeDiffKNNPairs(
                           SBLoc::latLngToCart3DPt((c+0.5) * lngInc - 180.0,
                                                    (r+0.5) * latInc - 90.0),
                           SBLoc::xyzDistFromLngLat(r*latInc-90.0,
                                                    (r+1)*latInc-90.0, lngInc),
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
    size_t multiLocs = rowSize*colSize - singleLocs;
    std::cout << "ave tree size: " << totalTreeSize/(rowSize*colSize)
    << "\nSingle loc cells: " << singleLocs
    << "\nMulti-loc cells:" << multiLocs << std::endl;
}

double UpgradeBKDTGridSBSolver::calcLatIncFromAlpc() {
    double surfaceArea = 4*M_PI*SBLoc::EARTH_RADIUS*SBLoc::EARTH_RADIUS;
    double numCells = numLocs/AVE_LOC_PER_CELL;
    double sideLen = sqrt(surfaceArea/numCells);
    return std::fabs(SBLoc::latFromHavDist(sideLen, 0));
}

void UpgradeBKDTGridSBSolver::build() {
    std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>> kdtData;
    kdtData.reserve(sbData->size());
    std::transform(sbData->begin(), sbData->end(), std::back_inserter(kdtData),
                   [&](const SBLoc& loc){ return
                       std::make_pair(SBLoc::latLngToCart3DPt(loc.lng, loc.lat), &loc);});
    sbKdt = KDTree<3, const SBLoc*, DistType::EUC>(kdtData.begin(), kdtData.end());
    numLocs = sbKdt.size();
    latInc = calcLatIncFromAlpc();
    rowSize = std::ceil(180.0/(latInc -
                               std::numeric_limits<double>::epsilon()*256));
    colSize = rowSize;
    lngInc = 360.0/colSize + std::numeric_limits<double>::epsilon()*256;
    fillGridCache();
}

const SBLoc* UpgradeBKDTGridSBSolver::findNearest(double lng, double lat) const {
    size_t idx = std::floor((lat+90.0)/latInc)*colSize + (lng+180.0)/lngInc;
    auto singleLoc = gridSingleCache[idx];
    return singleLoc ? singleLoc :
           gridTreeCache[idx].kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}






