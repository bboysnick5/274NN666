//
//  GridSBSolver<dist_type>.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/13/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "GridSBSolver.hpp"
#include <utility>
#include <algorithm>
#include <cmath>

template <typename dist_type>
const dist_type dist_type_MAX = std::numeric_limits<dist_type>::max();
template <typename dist_type>
const dist_type dist_type_MIN = std::numeric_limits<dist_type>::lowest();

template <typename dist_type>
GridSBSolver<dist_type>::GridSBSolver(dist_type alpc) :
    AVE_LOC_PER_CELL(alpc),
    minLng(dist_type_MAX<dist_type>), maxLng(dist_type_MIN<dist_type>),
    minLat(dist_type_MAX<dist_type>), maxLat(dist_type_MIN<dist_type>) {}

template <typename dist_type>
void GridSBSolver<dist_type>::findKeyLngLat(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    for (const auto &loc : *locData) {
        minLat = std::min(minLat, loc.geoPt[0]);
        maxLat = std::max(maxLat, loc.geoPt[0]);
        minLng = std::min(minLng, loc.geoPt[1]);
        maxLng = std::max(maxLng, loc.geoPt[1]);
    }
    minLng = std::floor(minLng); minLat = std::floor(minLat);
    maxLng = std::ceil(maxLng); maxLat = std::ceil(maxLat);
    midLng = (minLng + maxLng)/2.0; midLat = (minLat + maxLat)/2.0;
}

template <typename dist_type>
void GridSBSolver<dist_type>::constructGrid(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    rowSize = sqrt(locData->size()) / AVE_LOC_PER_CELL;
    sideLen = SBLoc<dist_type>::havDist({minLat, 0.0}, {maxLat, 0.0}) / rowSize;
    dist_type lowestLatCircleRadius = SBLoc<dist_type>::EARTH_RADIUS * cos(minLat);
    dist_type longestColDistSpan = lowestLatCircleRadius *
                                std::fabs(maxLng - minLng);
    colSize = longestColDistSpan/sideLen + 1;
    grid = std::vector<std::vector<std::unordered_set<const SBLoc<dist_type>*>>>(rowSize,
           std::vector<std::unordered_set<const SBLoc<dist_type>*>>(colSize));
}

template <typename dist_type>
std::pair<size_t, size_t> GridSBSolver<dist_type>::getIdx(dist_type lng, dist_type lat) const {
    dist_type unsignedRowDistFromCenter = SBLoc<dist_type>::havDist({lat, lng}, {midLat, lng}),
           unsignedColDistFromCenter = SBLoc<dist_type>::havDist({lat, lng}, {lat, midLng});
    return std::make_pair((lat < midLat ? -unsignedRowDistFromCenter :
                           unsignedRowDistFromCenter)/sideLen + rowSize/2,
                          (lng < midLng ? -unsignedColDistFromCenter :
                           unsignedColDistFromCenter)/sideLen + colSize/2);
}

template <typename dist_type>
void GridSBSolver<dist_type>::fillGrid(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    for (const auto &l : *locData) {
        auto idxPr = getIdx(l.geoPt[1], l.geoPt[0]);
        auto &cell = grid[idxPr.first][idxPr.second];
        numLocs = numLocs - static_cast<int>(cell.erase(&l)) + 1;
        cell.insert(&l);
    }
}

template <typename dist_type>
void GridSBSolver<dist_type>::build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData) {
    findKeyLngLat(locData);
    constructGrid(locData);
    fillGrid(locData);
}

template <typename dist_type>
void GridSBSolver<dist_type>::NNOneCell(const std::unordered_set<const SBLoc<dist_type>*> &cell,
dist_type lng, dist_type lat, dist_type &minDist, const SBLoc<dist_type>* &best) const {
    for (const auto &l : cell) {
        dist_type dist = SBLoc<dist_type>::havDist({lat, lng}, l->geoPt);
        if (dist < minDist) {
            minDist = dist;
            best = l;
        }
    }
}

template <typename dist_type>
const SBLoc<dist_type>* GridSBSolver<dist_type>::findNearest(const Point<dist_type, 2>& geoSearchPt) const {
    auto idxPr = getIdx(geoSearchPt[1], geoSearchPt[0]);
    size_t r0 = idxPr.first, c0 = idxPr.second;
    dist_type minDist = dist_type_MAX<dist_type>;
    const SBLoc<dist_type>* best = nullptr;
    // exact cell check
    NNOneCell(grid[r0][c0], geoSearchPt[1], geoSearchPt[0], minDist, best);
    
    // Spiral Search
    for (size_t d = 1; ; ++d) {
        for (int r = std::max(0,static_cast<int>(r0-d)); r <= std::min(static_cast<int>(r0+d), static_cast<int>(rowSize)-1); ++r) {
            if (c0-d >= 0)
                NNOneCell(grid[r][c0-d], geoSearchPt[1], geoSearchPt[0], minDist, best);
            if (c0+d < colSize)
                NNOneCell(grid[r][c0+d], geoSearchPt[1], geoSearchPt[0], minDist, best);
        }
        for (int c = std::max(0, static_cast<int>(c0-d+1)); c < std::min(static_cast<int>(c0+d), static_cast<int>(colSize)); c++) {
            if (r0 - d >= 0)
                NNOneCell(grid[r0-d][c], geoSearchPt[1], geoSearchPt[0], minDist, best);
            if (r0 + d < rowSize)
                NNOneCell(grid[r0+d][c], geoSearchPt[1], geoSearchPt[0], minDist, best);
        }
        if (minDist < d*sideLen*DISTORT_FACTOR)
            return best;
    }
    return best;
}

template class GridSBSolver<double>;
template class GridSBSolver<float>;




