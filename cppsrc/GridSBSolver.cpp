//
//  GridSBSolver.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/13/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "GridSBSolver.hpp"
#include <utility>
#include <algorithm>
#include <cmath>



std::vector<double> findBoundaryPtsInLngLat(const std::vector<SBLoc> &sbData) {
    double minLng = std::numeric_limits<double>::min(),
           minLat = minLng, maxLng = minLng, maxLat = minLng;
    for (const auto &loc : sbData) {
        minLng = std::min(minLng, loc.lng);
        minLat = std::min(minLat, loc.lat);
        maxLng = std::max(maxLng, loc.lng);
        maxLat = std::max(maxLat, loc.lat);
    }
    return {minLng, maxLng, minLat, maxLat};
}


void GridSBSolver::constructGrid(const std::vector<SBLoc> &sbData,
                                 const std::vector<double>& boundaryPts) {
    midLng = (boundaryPts[0] + boundaryPts[1])/2;
    midLat = (boundaryPts[2] + boundaryPts[3])/2;
    
    ySize = sqrt(sbData.size()) + 2;
    sideLen = SBLoc::havDist(boundaryPts[0], 0, boundaryPts[1], 0)/ySize;
    xSize = SBLoc::havDist(0, 0, 0, 180)/sideLen + 2;
    grid = std::vector<std::vector<std::unordered_set<SBLoc>>>(ySize,
           std::vector<std::unordered_set<SBLoc>>(xSize));
}

std::pair<int, int> GridSBSolver::getIdx(double lng, double lat) const {
    return std::make_pair(SBLoc::havDist(lng, lat, lng, midLat)/sideLen + xSize/2,
                          SBLoc::havDist(lng, lat, midLng, lat)/sideLen + ySize/2);
}

void GridSBSolver::build(const std::vector<SBLoc> &sbData) {
    auto boundaryPts = findBoundaryPtsInLngLat(sbData);
    constructGrid(sbData, boundaryPts);
    for (auto l : sbData) {
        auto idxPr = getIdx(l.lng, l.lat);
        auto &cell = grid[idxPr.first][idxPr.second];
        cell.erase(l);
        cell.insert(l);
    }
}

SBLoc GridSBSolver::findNearest(double lng, double lat) {
    return SBLoc();
}





