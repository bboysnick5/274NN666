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



std::vector<double>
GridSBSolver::findMaxLngLat(const std::vector<SBLoc> &sbData) {
    double minLng = std::numeric_limits<double>::max(), minLat = minLng;
    double maxLng = std::numeric_limits<double>::lowest(), maxLat = maxLng;
    for (const auto &loc : sbData) {
        minLng = std::min(minLng, loc.lng);
        minLat = std::min(minLat, loc.lat);
        maxLng = std::max(maxLng, loc.lng);
        maxLat = std::max(maxLat, loc.lat);
    }
    return {std::floor(minLng), std::ceil(maxLng),
            std::floor(minLat), std::ceil(maxLat)};
}


void GridSBSolver::constructGrid(const std::vector<SBLoc> &sbData,
                                 const std::vector<double>& boundaryPts) {
    midLng = (boundaryPts[0] + boundaryPts[1])/2;
    midLat = (boundaryPts[2] + boundaryPts[3])/2;
    
    rowSize = sqrt(sbData.size());
    rowSize *= 3;
    sideLen = SBLoc::havDist(0, boundaryPts[2], 0, boundaryPts[3])/rowSize;
    double lowestLatCircleRadius = SBLoc::EARTH_RADIUS *
                                   cos(boundaryPts[2]/180*M_PI);
    double longestColDistSpan = 2 * M_PI * lowestLatCircleRadius *
                                (std::fabs(boundaryPts[1]-boundaryPts[0])/360);
    colSize = longestColDistSpan/sideLen + 1;
    grid = std::vector<std::vector<std::unordered_set<SBLoc>>>(rowSize,
           std::vector<std::unordered_set<SBLoc>>(colSize));
}

std::pair<int, int> GridSBSolver::getIdx(double lng, double lat) const {
    double unsignedRowDistFromCenter = SBLoc::havDist(lng, lat, lng, midLat),
           unsignedColDistFromCenter = SBLoc::havDist(lng, lat, midLng, lat);
    //double z = SBLoc::havDist(midLng, midLat, lng, lat);
    //unsignedColDistFromCenter = sqrt(z*z-unsignedRowDistFromCenter*unsignedRowDistFromCenter);
    //double aveLat = (midLat + lat)/2;
    //unsignedColDistFromCenter = SBLoc::havDist(lng, aveLat, midLng, aveLat);

    return std::make_pair((lat < midLat ? -unsignedRowDistFromCenter :
                           unsignedRowDistFromCenter)/sideLen + rowSize/2,
                          (lng < midLng ? -unsignedColDistFromCenter :
                           unsignedColDistFromCenter)/sideLen + colSize/2);
}

void GridSBSolver::fillGrid(const std::vector<SBLoc>& sbData) {
    for (auto l : sbData) {
        auto idxPr = getIdx(l.lng, l.lat);
        auto &cell = grid[idxPr.first][idxPr.second];
        cell.erase(l);
        cell.insert(l);
    }
}

void GridSBSolver::build(const std::vector<SBLoc> &sbData) {
    auto minMaxLngLat = findMaxLngLat(sbData);
    constructGrid(sbData, minMaxLngLat);
    fillGrid(sbData);
}


void GridSBSolver::NNOneCell(const std::unordered_set<SBLoc> &cell, double lng,
                             double lat, double &minDist, SBLoc &best) const {
    for (const auto &l : cell) {
        double dist = SBLoc::havDist(lng, lat, l.lng, l.lat);
        if (dist < minDist) {
            minDist = dist;
            best = l;
        }
    }
}


SBLoc GridSBSolver::findNearest(double lng, double lat) const {
    auto idxPr = getIdx(lng, lat);
    int r0 = idxPr.first, c0 = idxPr.second;
    double minDist = std::numeric_limits<double>::max();
    SBLoc best;
    // exact cell check
    NNOneCell(grid[r0][c0], lng, lat, minDist, best);
    
    for (int d = 1; ; ++d) {
        for (int r = std::max(0, r0-d); r <= std::min(r0+d, rowSize-1); ++r) {
            if (c0-d >= 0)
                NNOneCell(grid[r][c0-d], lng, lat, minDist, best);
            if (c0+d < colSize)
                NNOneCell(grid[r][c0+d], lng, lat, minDist, best);
        }
        for (int c = std::max(0, c0-d+1); c < std::min(c0+d, colSize); c++) {
            if (r0 - d >= 0)
                NNOneCell(grid[r0-d][c], lng, lat, minDist, best);
            if (r0 + d < rowSize)
                NNOneCell(grid[r0+d][c], lng, lat, minDist, best);
        }
        if (minDist < d*sideLen*0.95)
            return best;
    }
    return best;
}





