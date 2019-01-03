//
//  BKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/23/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//


/*
 Build Version 1
 Inaccurate and waste of resource.
 Precalculate the side length in HavDist. Then apply that to each cell.
 */

#include <thread>
//#include <omp.h>
#include "BKDTGridSBSolver.hpp"


BKDTGridSBSolver::BKDTGridSBSolver(double alpc) : GridSBSolver(alpc) {}

std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>>::iterator
BKDTGridSBSolver::cacheAllPossibleLocsOneCell(size_t r0, size_t c0, double diff,
    std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>>::iterator begin) {
    double rowDistCellCtrGridCtr = ((r0+ 0.5 -rowSize/2 )*sideLen),
           colDistCellCtrGridCtr = ((c0 + 0.5 -colSize/2)*sideLen);
    double cellCtrLat = SBLoc::latFromHavDist(rowDistCellCtrGridCtr, midLat),
           cellCtrLng = SBLoc::lngFromHavDist(colDistCellCtrGridCtr,
                                              midLng, cellCtrLat);
    return sbKdt.rangeDiffKNNPairs(SBLoc::latLngToCart3DPt(cellCtrLng,
                                                           cellCtrLat),
                                   diff, begin);
}

void BKDTGridSBSolver::fillGridCache() {
    gridTreeCache = std::vector<KDTree<3, const SBLoc*, DistType::EUC>>
                    (rowSize*colSize);
    gridSingleCache = std::vector<const SBLoc*>(rowSize*colSize);
    std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>> ptLocPairs(numLocs);
    size_t totalTreeSize = 0, singleLocs = 0;
    double diff = xyzDistFromSideLen();
//#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
//default(none) schedule(guided) shared(diff) firstprivate(ptLocPairs) \
//reduction(+:totalTreeSize, singleLocs) collapse(2)
    for (size_t r = 0; r < rowSize; ++r) {
        for (size_t c = 0; c < colSize; ++c) {
            size_t idx = r*colSize+c;
            auto locsEnd = cacheAllPossibleLocsOneCell(r, c, diff,
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

double BKDTGridSBSolver::xyzDistFromSideLen() {
    double lat2 = SBLoc::latFromHavDist(sideLen*sqrt(2), 0);
    return Point<3, DistType::EUC>::dist(SBLoc::latLngToCart3DPt(0, 0),
                                         SBLoc::latLngToCart3DPt(0, lat2)); 
}
 

void BKDTGridSBSolver::build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>> kdtData;
    kdtData.reserve(locData->size());
    std::transform(locData->begin(), locData->end(), std::back_inserter(kdtData),
                   [&](const SBLoc& loc){ return
                       std::make_pair(SBLoc::latLngToCart3DPt(loc.lng, loc.lat), &loc);});
    sbKdt = KDTree<3, const SBLoc*, DistType::EUC>(kdtData.begin(), kdtData.end());
    numLocs = sbKdt.size();
    findKeyLngLat(locData);
    rowSize = sqrt(locData->size()) / AVE_LOC_PER_CELL;
    sideLen = SBLoc::havDist(0, minLat, 0, maxLat) / rowSize;
    double lowestLatCircleRadius = SBLoc::EARTH_RADIUS * cos(minLat < 0 ? 0 : minLat);
    double longestColDistSpan = 2 * M_PI * lowestLatCircleRadius *
                                (std::fabs(maxLng - minLng)/(2*M_PI));
    colSize = longestColDistSpan/sideLen + 1;
    
    fillGridCache();
}

const SBLoc* BKDTGridSBSolver::findNearest(double lng, double lat) const {
    auto idxPr = getIdx(lng, lat);
    size_t idx = idxPr.first*colSize+idxPr.second;
    auto singleLoc = gridSingleCache[idx];
    return singleLoc != nullptr ? singleLoc :
           gridTreeCache[idx].kNNValue(SBLoc::latLngToCart3DPt(lng, lat), 1);
}






