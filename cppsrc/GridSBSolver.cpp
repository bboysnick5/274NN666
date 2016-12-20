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
    sideLen = SBLoc::havDist(0, boundaryPts[2], 0, boundaryPts[3])/rowSize;
    //rowSize += 200;
    colSize = SBLoc::havDist(boundaryPts[0], boundaryPts[2],
                             boundaryPts[1], boundaryPts[2])/sideLen;
    grid = std::vector<std::vector<std::unordered_set<SBLoc>>>(rowSize,
           std::vector<std::unordered_set<SBLoc>>(colSize));
}

std::pair<int, int> GridSBSolver::getIdx(double lng, double lat) const {
    double unsignedRowDistFromCenter = SBLoc::havDist(lng, lat, lng, midLat),
    unsignedColDistFromCenter = SBLoc::havDist(lng, lat, midLng, lat);
    return std::make_pair((lat < midLat ? -unsignedRowDistFromCenter
                           : unsignedRowDistFromCenter)/sideLen + rowSize/2,
                          (lng < midLng ? -unsignedColDistFromCenter :
                           unsignedColDistFromCenter)/sideLen + colSize/2);
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
    
    // ------------------ DEBUG USE --------------------
    // -------------------------------------------------
    std::unordered_set<size_t> hashes;
    int count = 0;
    for (int x = 0; x < grid.size(); ++x) {
        for (int y = 0; y < grid[0].size(); ++y) {
            auto &set = grid[x][y];
            if (!set.empty()) {
                for (auto &l : set) {
                    hashes.insert(std::hash<SBLoc>()(l));
                }
                std::cout << x << " " << y << "  set size: " <<
                set.size() << std::endl;
                
            }
            count += set.size();
        }
    }
    std::cout << "size: " << count << std::endl;
    std::cout << "hashes size: " << hashes.size() << std::endl;
}


void GridSBSolver::NNInCell(const std::unordered_set<SBLoc> &cell, double lng,
                              double lat, double &minDist, SBLoc &best) {
    for (const auto &l : cell) {
        double dist = SBLoc::havDist(lng, lat, l.lng, l.lat);
        if (dist < minDist) {
            minDist = dist;
            best = l;
        }
    }
}


SBLoc GridSBSolver::findNearest(double lng, double lat) {
    // todo: out of boundary
    
    auto idxPr = getIdx(lng, lat);
    int r0 = idxPr.first, c0 = idxPr.second;
    double minDist = std::numeric_limits<double>::max();
    SBLoc best;
    // exact cell check
    NNInCell(grid[r0][c0], lng, lat, minDist, best);
    
    for (int d = 0; ; ++d) {
        for (int r = std::max(0, r0-d); r <= std::min(r0+d, rowSize-1); ++r) {
            if (c0-d >= 0)
                NNInCell(grid[r][c0-d], lng, lat, minDist, best);
            if (c0+d < colSize)
                NNInCell(grid[r][c0+d], lng, lat, minDist, best);
        }
        for (int c = std::max(0, c0-d+1); c < std::min(c0+d, colSize); c++) {
            if (r0 - d >= 0)
                NNInCell(grid[r0-d][c], lng, lat, minDist, best);
            if (r0 + d < rowSize)
                NNInCell(grid[r0+d][c], lng, lat, minDist, best);
        }
        if (minDist < d*sideLen)
            return best;
    }
    
    return best;
}





