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


const double DOUBLE_MAX = std::numeric_limits<double>::max();
const double DOUBLE_MIN = std::numeric_limits<double>::lowest();

GridSBSolver::GridSBSolver(double alpc) :
    AVE_LOC_PER_CELL(alpc),
    minLng(DOUBLE_MAX), maxLng(DOUBLE_MIN),
    minLat(DOUBLE_MAX), maxLat(DOUBLE_MIN) {}

void GridSBSolver::findKeyLngLat(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    for (const auto &loc : *locData) {
        minLng = std::min(minLng, loc.lng);
        minLat = std::min(minLat, loc.lat);
        maxLng = std::max(maxLng, loc.lng);
        maxLat = std::max(maxLat, loc.lat);
    }
    minLng = std::floor(minLng); minLat = std::floor(minLat);
    maxLng = std::ceil(maxLng); maxLat = std::ceil(maxLat);
    midLng = (minLng + maxLng)/2; midLat = (minLat + maxLat)/2;
}

void GridSBSolver::constructGrid(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    rowSize = sqrt(locData->size()) / AVE_LOC_PER_CELL;
    sideLen = SBLoc::havDist(0, minLat, 0, maxLat) / rowSize;
    double lowestLatCircleRadius = SBLoc::EARTH_RADIUS * cos(minLat);
    double longestColDistSpan = lowestLatCircleRadius *
                                std::fabs(maxLng - minLng);
    colSize = longestColDistSpan/sideLen + 1;
    grid = std::vector<std::vector<std::unordered_set<const SBLoc*>>>(rowSize,
           std::vector<std::unordered_set<const SBLoc*>>(colSize));
}

std::pair<size_t, size_t> GridSBSolver::getIdx(double lng, double lat) const {
    double unsignedRowDistFromCenter = SBLoc::havDist(lng, lat, lng, midLat),
           unsignedColDistFromCenter = SBLoc::havDist(lng, lat, midLng, lat);
    return std::make_pair((lat < midLat ? -unsignedRowDistFromCenter :
                           unsignedRowDistFromCenter)/sideLen + rowSize/2,
                          (lng < midLng ? -unsignedColDistFromCenter :
                           unsignedColDistFromCenter)/sideLen + colSize/2);
}

void GridSBSolver::fillGrid(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    for (const auto &l : *locData) {
        auto idxPr = getIdx(l.lng, l.lat);
        auto &cell = grid[idxPr.first][idxPr.second];
        numLocs = numLocs - static_cast<int>(cell.erase(&l)) + 1;
        cell.insert(&l);
    }
}

void GridSBSolver::build(const std::shared_ptr<std::vector<SBLoc>> &locData) {
    findKeyLngLat(locData);
    constructGrid(locData);
    fillGrid(locData);
}

void GridSBSolver::NNOneCell(const std::unordered_set<const SBLoc*> &cell,
double lng, double lat, double &minDist, const SBLoc* &best) const {
    for (const auto &l : cell) {
        double dist = SBLoc::havDist(lng, lat, l->lng, l->lat);
        if (dist < minDist) {
            minDist = dist;
            best = l;
        }
    }
}


const SBLoc* GridSBSolver::findNearest(double lng, double lat) const {
    auto idxPr = getIdx(lng, lat);
    size_t r0 = idxPr.first, c0 = idxPr.second;
    double minDist = DOUBLE_MAX;
    const SBLoc* best = nullptr;
    // exact cell check
    NNOneCell(grid[r0][c0], lng, lat, minDist, best);
    
    // Spiral Search
    for (size_t d = 1; ; ++d) {
        for (int r = std::max(0,static_cast<int>(r0-d)); r <= std::min(static_cast<int>(r0+d), static_cast<int>(rowSize)-1); ++r) {
            if (c0-d >= 0)
                NNOneCell(grid[r][c0-d], lng, lat, minDist, best);
            if (c0+d < colSize)
                NNOneCell(grid[r][c0+d], lng, lat, minDist, best);
        }
        for (int c = std::max(0, static_cast<int>(c0-d+1)); c < std::min(static_cast<int>(c0+d), static_cast<int>(colSize)); c++) {
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





