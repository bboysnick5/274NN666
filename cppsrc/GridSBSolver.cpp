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

GridSBSolver::GridSBSolver(double alpc) :
    AVE_LOC_PER_CELL(alpc),
    maxLng(std::numeric_limits<double>::lowest()),
    maxLat(std::numeric_limits<double>::lowest()),
    minLng(std::numeric_limits<double>::max()),
    minLat(std::numeric_limits<double>::max()) {}

void GridSBSolver::findKeyLngLat(const std::vector<SBLoc> &sbData) {
    for (const auto &loc : sbData) {
        minLng = std::min(minLng, loc.lng);
        minLat = std::min(minLat, loc.lat);
        maxLng = std::max(maxLng, loc.lng);
        maxLat = std::max(maxLat, loc.lat);
    }
    minLng = std::floor(minLng), minLat = std::floor(minLat);
    maxLng = std::ceil(maxLng), maxLat = std::ceil(maxLat);
    midLng = (minLng + maxLng)/2, midLat = (minLat + maxLat)/2;
}

void GridSBSolver::constructGrid(const std::vector<SBLoc> &sbData) {
    rowSize = sqrt(sbData.size()) / AVE_LOC_PER_CELL;
    sideLen = SBLoc::havDist(0, minLat, 0, maxLat) / rowSize;
    double lowestLatCircleRadius = SBLoc::EARTH_RADIUS * cos(minLat/180*M_PI);
    double longestColDistSpan = 2 * M_PI * lowestLatCircleRadius *
                                (std::fabs(maxLng - minLng)/360);
    colSize = longestColDistSpan/sideLen + 1;
    grid = std::vector<std::vector<std::unordered_set<SBLoc>>>(rowSize,
           std::vector<std::unordered_set<SBLoc>>(colSize));
}

std::pair<int, int> GridSBSolver::getIdx(double lng, double lat) const {
    double unsignedRowDistFromCenter = SBLoc::havDist(lng, lat, lng, midLat),
           unsignedColDistFromCenter = SBLoc::havDist(lng, lat, midLng, lat);
    return std::make_pair((lat < midLat ? -unsignedRowDistFromCenter :
                           unsignedRowDistFromCenter)/sideLen + rowSize/2,
                          (lng < midLng ? -unsignedColDistFromCenter :
                           unsignedColDistFromCenter)/sideLen + colSize/2);
}

void GridSBSolver::fillGrid(const std::vector<SBLoc>& sbData) {
    for (auto l : sbData) {
        auto idxPr = getIdx(l.lng, l.lat);
        auto &cell = grid[idxPr.first][idxPr.second];
        numLocs = numLocs - static_cast<int>(cell.erase(l)) + 1;
        cell.insert(l);
    }
}

void GridSBSolver::build(const std::vector<SBLoc> &sbData) {
    findKeyLngLat(sbData);
    constructGrid(sbData);
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
    
    // Spiral Search
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
        if (minDist < d*sideLen*DISTORT_FACTOR)
            return best;
    }
    return best;
}





