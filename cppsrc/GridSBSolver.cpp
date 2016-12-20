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
    return {minLng, maxLng, minLat, maxLat};
}


void GridSBSolver::constructGrid(const std::vector<SBLoc> &sbData,
                                 const std::vector<double>& boundaryPts) {
    midLng = (boundaryPts[0] + boundaryPts[1])/2;
    midLat = (boundaryPts[2] + boundaryPts[3])/2;
    
    rowSize = sqrt(sbData.size());
    sideLen = SBLoc::havDist(0, boundaryPts[2], 0, boundaryPts[3])/rowSize;
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

SBLoc GridSBSolver::findNearest(double lng, double lat) {
    /*
    // (di, dj) is a vector - direction in which we move right now
    int di = 1;
    int dj = 0;
    // length of current segment
    int segment_length = 1;
    
    // current position (i, j) and how much of current segment we passed
    int i = 0;
    int j = 0;
    int segment_passed = 0;
    for (int k = 0; k < NUMBER_OF_POINTS; ++k) {
        // make a step, add 'direction' vector (di, dj) to current position (i, j)
        i += di;
        j += dj;
        ++segment_passed;
        System.out.println(i + " " + j);
        
        if (segment_passed == segment_length) {
            // done with current segment
            segment_passed = 0;
            
            // 'rotate' directions
            int buffer = di;
            di = -dj;
            dj = buffer;
            
            // increase segment length if necessary
            if (dj == 0) {
                ++segment_length;
            }
        }
    }
    */
    
    return SBLoc();
}





