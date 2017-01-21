//
//  BKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/23/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "BKDTGridSBSolver.hpp"

/*
double BKDTGridSBSolver::distFromMidPtInCell(double lng, double lat,
            double rowDistCellCtrGridCtr, double colDistCellCtrGridCtr) const {
    double rowDistFromGridCenter = SBLoc::havDist(lng, lat, lng, midLat),
           colDistFromGridCenter = SBLoc::havDist(lng, lat, midLng, lat);
    double rowDiff = (lat < midLat ? -rowDistFromGridCenter
                      : rowDistFromGridCenter) - rowDistCellCtrGridCtr,
           colDiff = (lng < midLng ? -colDistFromGridCenter
                      : colDistFromGridCenter) - colDistCellCtrGridCtr;
    return sqrt(rowDiff*rowDiff + colDiff*colDiff);
}*/

void BKDTGridSBSolver::checkOneCell(const std::unordered_set<SBLoc>& cell,
double cellCtrLng, double cellCtrLat, double& minDist,
std::vector<std::pair<Point<3>, SBLoc>>& validLocs) const {
    for (const auto &l : cell) {
        double dist = SBLoc::havDist(cellCtrLng, cellCtrLat, l.lng, l.lat);
        minDist = std::min(minDist, dist);
        validLocs.emplace_back(SBLoc::latLngToCart3DXYZ(l.lng, l.lat), l);
    }
}

int totalTreeSize = 0;
void BKDTGridSBSolver::fillCacheOneCell(int r0, int c0) {
    double rowDistCellCtrGridCtr = ((r0-rowSize/2)*sideLen + sideLen/2),
           colDistCellCtrGridCtr = ((c0-colSize/2)*sideLen + sideLen/2);
    double cellCtrLat = SBLoc::latFromHavDist(rowDistCellCtrGridCtr, midLat),
           cellCtrLng = SBLoc::lngFromHavDist(colDistCellCtrGridCtr,
                                              midLng, cellCtrLat);
    double minDist = std::numeric_limits<double>::max();
    std::vector<std::pair<Point<3>, SBLoc>> validLocs;
    int maxD = std::numeric_limits<int>::max();
    // exact cell check
    checkOneCell(grid[r0][c0], cellCtrLng, cellCtrLat, minDist, validLocs);
    for (int d = 1; d <= maxD; ++d) {
        double thisMinDist = std::numeric_limits<double>::max();
        std::vector<std::pair<Point<3>, SBLoc>> thisValidLocs;


        for (int r = std::max(0, r0-d); r <= std::min(r0+d, rowSize-1); ++r) {
            if (c0-d >= 0)
                checkOneCell(grid[r][c0-d], cellCtrLng, cellCtrLat, thisMinDist, thisValidLocs);
            if (c0+d < colSize)
                checkOneCell(grid[r][c0+d], cellCtrLng, cellCtrLat, thisMinDist, thisValidLocs);
        }
        for (int c = std::max(0, c0-d+1); c < std::min(c0+d, colSize); c++) {
            if (r0 - d >= 0)
                checkOneCell(grid[r0-d][c], cellCtrLng, cellCtrLat, thisMinDist, thisValidLocs);
            if (r0 + d < rowSize)
                checkOneCell(grid[r0+d][c], cellCtrLng, cellCtrLat, thisMinDist, thisValidLocs);
        }
        if (thisMinDist != std::numeric_limits<double>::max()
            && thisMinDist - minDist > sideLen * sqrt(2)*(1+d*0.35))
            break;
        minDist = std::min(minDist, thisMinDist);
        std::copy(thisValidLocs.begin(), thisValidLocs.end(), std::back_inserter(validLocs));

        
        
        /*
        int prevVLSize = validLocs.size();
        for (int r = std::max(0, r0-d); r <= std::min(r0+d, rowSize-1); ++r) {
            if (c0-d >= 0)
                checkOneCell(grid[r][c0-d], cellCtrLng, cellCtrLat, minDist, validLocs);
            if (c0+d < colSize)
                checkOneCell(grid[r][c0+d], cellCtrLng, cellCtrLat, minDist, validLocs);
        }
        for (int c = std::max(0, c0-d+1); c < std::min(c0+d, colSize); c++) {
            if (r0 - d >= 0)
                checkOneCell(grid[r0-d][c], cellCtrLng, cellCtrLat, minDist, validLocs);
            if (r0 + d < rowSize)
                checkOneCell(grid[r0+d][c], cellCtrLng, cellCtrLat, minDist, validLocs);
        }
        if (minDist != std::numeric_limits<double>::max()
            && maxD == std::numeric_limits<int>::max()) {
            if (minDist < (d-sqrt(2)/2)*sideLen)
                break;
            else
                maxD = d+1;
        }*/
    }
    validLocs.erase(std::remove_if(validLocs.begin(), validLocs.end(), [=](const std::pair<Point<3>, SBLoc> &locPair){return SBLoc::havDist(locPair.second.lng, locPair.second.lat, cellCtrLng, cellCtrLat) > minDist + sideLen * sqrt(2);}), validLocs.end());
    totalTreeSize += validLocs.size();
    gridCache[r0][c0] = KDTree<3, SBLoc>(validLocs.begin(), validLocs.end());
}

void BKDTGridSBSolver::fillGridCache() {
    gridCache = std::vector<std::vector<KDTree<3, SBLoc>>>(rowSize,
                std::vector<KDTree<3, SBLoc>>(colSize));
    for (int r = 0; r < rowSize; ++r) {
        for (int c = 0; c < colSize; ++c) {
            fillCacheOneCell(r, c);
        }
    }
}


void BKDTGridSBSolver::build(const std::vector<SBLoc> &sbData) {
    auto minMaxLngLat = findMaxLngLat(sbData);
    constructGrid(sbData, minMaxLngLat);
    fillGrid(sbData);
    fillGridCache();
    std::cout << "ave tree size: " << totalTreeSize/(rowSize*colSize) << "\n\n";
}

SBLoc BKDTGridSBSolver::findNearest(double lng, double lat) const {
    auto idxPr = getIdx(lng, lat);
    return gridCache[idxPr.first][idxPr.second].
           kNNValue(SBLoc::latLngToCart3DXYZ(lng, lat), 1);
}
