//
//  BKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/23/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "BKDTGridSBSolver.hpp"


BKDTGridSBSolver::BKDTGridSBSolver(double alpc) : GridSBSolver(alpc) {}

void BKDTGridSBSolver::checkOneCell(const std::unordered_set<const SBLoc*>& cell,
double cellCtrLng, double cellCtrLat, double& minDist,
std::vector<const SBLoc*>& validLocs) const {
    for (const auto &l : cell) {
        double dist = SBLoc::havDist(cellCtrLng, cellCtrLat, l->lng, l->lat);
        minDist = std::min(minDist, dist);
        validLocs.push_back(l);
    }
}

int totalTreeSize = 0;
std::vector<const SBLoc*>
BKDTGridSBSolver::spiralSearchAllPossibleLocsOneCell(int r0, int c0) {
    double rowDistCellCtrGridCtr = ((r0-rowSize/2)*sideLen + sideLen/2),
           colDistCellCtrGridCtr = ((c0-colSize/2)*sideLen + sideLen/2);
    double cellCtrLat = SBLoc::latFromHavDist(rowDistCellCtrGridCtr, midLat),
           cellCtrLng = SBLoc::lngFromHavDist(colDistCellCtrGridCtr,
                                              midLng, cellCtrLat);
    double minDist = DOUBLE_MAX, thisMinDist = DOUBLE_MAX;
    std::vector<const SBLoc*> validLocs;
    // exact cell check
    checkOneCell(grid[r0][c0], cellCtrLng, cellCtrLat, minDist, validLocs);
    for (int d = 1; validLocs.size() < numLocs; ++d) {
        for (int r = std::max(0, r0-d); r <= std::min(r0+d, rowSize-1); ++r) {
            if (c0-d >= 0)
                checkOneCell(grid[r][c0-d], cellCtrLng, cellCtrLat, thisMinDist, validLocs);
            if (c0+d < colSize)
                checkOneCell(grid[r][c0+d], cellCtrLng, cellCtrLat, thisMinDist, validLocs);
        }
        for (int c = std::max(0, c0-d+1); c < std::min(c0+d, colSize); c++) {
            if (r0 - d >= 0)
                checkOneCell(grid[r0-d][c], cellCtrLng, cellCtrLat, thisMinDist, validLocs);
            if (r0 + d < rowSize)
                checkOneCell(grid[r0+d][c], cellCtrLng, cellCtrLat, thisMinDist, validLocs);
        }
        if (thisMinDist != DOUBLE_MAX &&
            thisMinDist - minDist > sideLen*sqrt(2)*(1+d*0.5))
            break;
        minDist = std::min(minDist, thisMinDist);
        thisMinDist = DOUBLE_MAX;
    }
    validLocs.erase(std::remove_if(validLocs.begin(), validLocs.end(),
        [=](const SBLoc* l) {return SBLoc::havDist(l->lng, l->lat, cellCtrLng,
        cellCtrLat) > minDist + sideLen * sqrt(2);}), validLocs.end());
    totalTreeSize += validLocs.size();
    return validLocs;
}

void BKDTGridSBSolver::fillGridCache() {
    gridTreeCache = std::vector<std::vector<KDTree<3, const SBLoc*>>>(rowSize,
                std::vector<KDTree<3, const SBLoc*>>(colSize));
    for (int r = 0; r < rowSize; ++r) {
        for (int c = 0; c < colSize; ++c) {
            auto locs = spiralSearchAllPossibleLocsOneCell(r, c);
            std::vector<std::pair<Point<3>, const SBLoc*>> locTreePairs;
            std::transform(locs.begin(), locs.end(),
                std::back_inserter(locTreePairs), [](const SBLoc* l) {return
                std::make_pair(SBLoc::latLngToCart3DXYZ(l->lng, l->lat), l);});
            gridTreeCache[r][c] = KDTree<3, const SBLoc*>(locTreePairs.begin(),
                                                          locTreePairs.end());
        }
    }
}


void BKDTGridSBSolver::build() {
    findKeyLngLat();
    constructGrid();
    fillGrid();
    fillGridCache();
    std::cout << "ave tree size: " << totalTreeSize/(rowSize*colSize) << "\n\n";
}

const SBLoc* BKDTGridSBSolver::findNearest(double lng, double lat) const {
    auto idxPr = getIdx(lng, lat);
    return gridTreeCache[idxPr.first][idxPr.second].
           kNNValue(SBLoc::latLngToCart3DXYZ(lng, lat), 1);
}
