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
std::vector<const SBLoc*>& validLocs, int &vlIdx) const {
    for (const auto &l : cell) {
        double dist = SBLoc::havDist(cellCtrLng, cellCtrLat, l->lng, l->lat);
        minDist = std::min(minDist, dist);
        validLocs[vlIdx++] = l;
    }
}

std::vector<const SBLoc*>::iterator
BKDTGridSBSolver::spiralSearchAllPossibleLocsOneCell(int r0, int c0,
std::vector<const SBLoc*> &validLocs) {
    double rowDistCellCtrGridCtr = ((r0-rowSize/2)*sideLen + sideLen/2),
           colDistCellCtrGridCtr = ((c0-colSize/2)*sideLen + sideLen/2);
    double cellCtrLat = SBLoc::latFromHavDist(rowDistCellCtrGridCtr, midLat),
           cellCtrLng = SBLoc::lngFromHavDist(colDistCellCtrGridCtr,
                                              midLng, cellCtrLat);
    double minDist = DOUBLE_MAX, thisMinDist = DOUBLE_MAX;
    int vlIdx = 0;
    // exact cell check
    checkOneCell(grid[r0][c0], cellCtrLng, cellCtrLat, minDist, validLocs, vlIdx);
    for (int d = 1; vlIdx < numLocs; ++d) {
        for (int r = std::max(0, r0-d); r <= std::min(r0+d, rowSize-1); ++r) {
            if (c0-d >= 0)
                checkOneCell(grid[r][c0-d], cellCtrLng, cellCtrLat, thisMinDist, validLocs, vlIdx);
            if (c0+d < colSize)
                checkOneCell(grid[r][c0+d], cellCtrLng, cellCtrLat, thisMinDist, validLocs, vlIdx);
        }
        for (int c = std::max(0, c0-d+1); c < std::min(c0+d, colSize); c++) {
            if (r0 - d >= 0)
                checkOneCell(grid[r0-d][c], cellCtrLng, cellCtrLat, thisMinDist, validLocs, vlIdx);
            if (r0 + d < rowSize)
                checkOneCell(grid[r0+d][c], cellCtrLng, cellCtrLat, thisMinDist, validLocs, vlIdx);
        }
        if (thisMinDist != DOUBLE_MAX &&
            thisMinDist - minDist > sideLen*sqrt(2)*(1+d*0.5))
            break;
        minDist = std::min(minDist, thisMinDist);
        thisMinDist = DOUBLE_MAX;
    }
    return std::remove_if(validLocs.begin(), validLocs.begin() + vlIdx,
        [=](const SBLoc* l) {return SBLoc::havDist(l->lng, l->lat, cellCtrLng,
        cellCtrLat) > minDist + sideLen * sqrt(2);});
}

void BKDTGridSBSolver::fillGridCache() {
    gridTreeCache = std::vector<std::vector<KDTree<3, const SBLoc*>>>(rowSize,
                    std::vector<KDTree<3, const SBLoc*>>(colSize));
    gridSingleCache = std::vector<std::vector<const SBLoc*>>(rowSize,
                      std::vector<const SBLoc*>(colSize));
    std::vector<const SBLoc*> locs(numLocs, nullptr);
    std::vector<std::pair<Point<3>, const SBLoc*>> locTreePairs;
    int totalTreeSize = 0;
    for (int r = 0; r < rowSize; ++r) {
        for (int c = 0; c < colSize; ++c) {
            auto locsEnd = spiralSearchAllPossibleLocsOneCell(r, c, locs);
            size_t locsSize = locsEnd - locs.begin();
            if (locsSize > 1) {
                std::transform(locs.begin(), locsEnd, std::back_inserter
                    (locTreePairs), [](const SBLoc* l) {return std::make_pair
                    (SBLoc::latLngToCart3DXYZ(l->lng, l->lat), l);});
                gridTreeCache[r][c]= KDTree<3, const SBLoc*>
                    (locTreePairs.begin(), locTreePairs.end());
                locTreePairs.clear();
            } else {
                gridSingleCache[r][c] = locs[0];
            }
            totalTreeSize += locsSize;
        }
    }
    std::cout << "ave tree size: " << totalTreeSize/(rowSize*colSize) << "\n\n";
}


void BKDTGridSBSolver::build() {
    findKeyLngLat();
    constructGrid();
    fillGrid();
    fillGridCache();
}

const SBLoc* BKDTGridSBSolver::findNearest(double lng, double lat) const {
    auto idxPr = getIdx(lng, lat);
    auto singleLoc = gridSingleCache[idxPr.first][idxPr.second];
    return singleLoc != nullptr ? singleLoc :
           gridTreeCache[idxPr.first][idxPr.second].
               kNNValue(SBLoc::latLngToCart3DXYZ(lng, lat), 1);
}
