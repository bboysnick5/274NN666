//
//  BKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/23/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "BKDTGridSBSolver.hpp"


BKDTGridSBSolver::BKDTGridSBSolver(double alpc) : GridSBSolver(alpc) {}

std::vector<std::pair<Point<3>, const SBLoc*>>::iterator
BKDTGridSBSolver::cacheAllPossibleLocsOneCell(int r0, int c0, double diff,
    std::vector<std::pair<Point<3>, const SBLoc*>>::iterator begin) {
    double rowDistCellCtrGridCtr = ((r0-rowSize/2)*sideLen + sideLen/2),
           colDistCellCtrGridCtr = ((c0-colSize/2)*sideLen + sideLen/2);
    double cellCtrLat = SBLoc::latFromHavDist(rowDistCellCtrGridCtr, midLat),
           cellCtrLng = SBLoc::lngFromHavDist(colDistCellCtrGridCtr,
                                              midLng, cellCtrLat);
    return sbKdt.rangeDiffKNNPairs(SBLoc::latLngToCart3DXYZ(cellCtrLng,
                                                            cellCtrLat), diff, begin);
}

void BKDTGridSBSolver::fillGridCache() {
    gridTreeCache = std::vector<KDTree<3, const SBLoc*, DistType::EUC>>
                    (rowSize*colSize);
    gridSingleCache = std::vector<const SBLoc*>(rowSize*colSize);
    std::vector<std::pair<Point<3>, const SBLoc*>> ptLocPairs(numLocs);
    int totalTreeSize = 0, singleLocs = 0, multiLocs = 0;
    double diff = xyzDistFromSideLen();
    for (int r = 0; r < rowSize; ++r) {
        for (int c = 0; c < colSize; ++c) {
            int idx = r*colSize+c;
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
    multiLocs = rowSize*colSize - singleLocs;
    std::cout << "ave tree size: " << totalTreeSize/(rowSize*colSize)
              << "\nSingle loc cells: " << singleLocs
              << "\nMulti-loc cells:" << multiLocs << std::endl;
}

double BKDTGridSBSolver::xyzDistFromSideLen() {
    double lat2 = SBLoc::latFromHavDist(sideLen*sqrt(2), 45);
    return Point<3>::euclDist(SBLoc::latLngToCart3DXYZ(0, 45),
                              SBLoc::latLngToCart3DXYZ(0, lat2)); 
}


void BKDTGridSBSolver::build() {
    std::vector<std::pair<Point<3>, const SBLoc*>> kdtData;
    kdtData.reserve(sbData->size());
    std::transform(sbData->begin(), sbData->end(), std::back_inserter(kdtData),
                   [&](const SBLoc& loc){ return
                       std::make_pair(SBLoc::latLngToCart3DXYZ(loc.lng, loc.lat), &loc);});
    sbKdt = KDTree<3, const SBLoc*, DistType::EUC>(kdtData.begin(), kdtData.end());
    numLocs = static_cast<int>(sbKdt.size());
    findKeyLngLat();
    rowSize = sqrt(sbData->size()) / AVE_LOC_PER_CELL;
    sideLen = SBLoc::havDist(0, minLat, 0, maxLat) / rowSize;
    double lowestLatCircleRadius = SBLoc::EARTH_RADIUS * cos(minLat/180*M_PI);
    double longestColDistSpan = 2 * M_PI * lowestLatCircleRadius *
                                (std::fabs(maxLng - minLng)/360);
    colSize = longestColDistSpan/sideLen + 1;
    fillGridCache();
}

const SBLoc* BKDTGridSBSolver::findNearest(double lng, double lat) const {
    auto idxPr = getIdx(lng, lat);
    int idx = idxPr.first*colSize+idxPr.second;
    auto singleLoc = gridSingleCache[idx];
    return singleLoc != nullptr ? singleLoc :
           gridTreeCache[idx].kNNValue(SBLoc::latLngToCart3DXYZ(lng, lat), 1);
}






