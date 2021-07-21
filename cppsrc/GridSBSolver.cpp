//
//  GridSBSolver<FPType, policy>.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/13/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "GridSBSolver.hpp"
#include <utility>
#include <algorithm>
#include <cmath>

template <typename FPType>
const FPType FPType_MAX = std::numeric_limits<FPType>::max();
template <typename FPType>
const FPType FPType_MIN = std::numeric_limits<FPType>::lowest();

template <typename FPType, def::ThreadingPolicy policy>
GridSBSolver<FPType, policy>::GridSBSolver(FPType alpc) :
    AVE_LOC_PER_CELL(alpc),
    minLng(FPType_MAX<FPType>), maxLng(FPType_MIN<FPType>),
    minLat(FPType_MAX<FPType>), maxLat(FPType_MIN<FPType>) {}

template <typename FPType, def::ThreadingPolicy policy>
void GridSBSolver<FPType, policy>::findKeyLngLat(std::span<const SBLoc<FPType>> loc_data_span) {
    for (const auto &loc : loc_data_span) {
        minLat = std::min(minLat, loc.geo_pt[0]);
        maxLat = std::max(maxLat, loc.geo_pt[0]);
        minLng = std::min(minLng, loc.geo_pt[1]);
        maxLng = std::max(maxLng, loc.geo_pt[1]);
    }
    minLng = std::floor(minLng); minLat = std::floor(minLat);
    maxLng = std::ceil(maxLng); maxLat = std::ceil(maxLat);
    midLng = (minLng + maxLng)/2.0; midLat = (minLat + maxLat)/2.0;
}

template <typename FPType, def::ThreadingPolicy policy>
void GridSBSolver<FPType, policy>::constructGrid(std::span<const SBLoc<FPType>> loc_data_span) {
    row_size_ = sqrt(loc_data_span.size()) / AVE_LOC_PER_CELL;
    side_len_ = SBLoc<FPType>::havDist({minLat, 0.0}, {maxLat, 0.0}) / row_size_;
    FPType lowestLatCircleRadius = SBLoc<FPType>::EARTH_RADIUS * std::cos(minLat);
    FPType longestColDistSpan = lowestLatCircleRadius *
                                std::fabs(maxLng - minLng);
    col_size_ = longestColDistSpan/side_len_ + 1;
    grid_ = std::vector<std::vector<std::unordered_set<const SBLoc<FPType>*>>>(row_size_,
           std::vector<std::unordered_set<const SBLoc<FPType>*>>(col_size_));
}

template <typename FPType, def::ThreadingPolicy policy>
std::pair<std::size_t, std::size_t> GridSBSolver<FPType, policy>::getIdx(FPType lng, FPType lat) const {
    FPType unsignedRowDistFromCenter = SBLoc<FPType>::havDist({lat, lng}, {midLat, lng}),
           unsignedColDistFromCenter = SBLoc<FPType>::havDist({lat, lng}, {lat, midLng});
    return std::make_pair((lat < midLat ? -unsignedRowDistFromCenter :
                           unsignedRowDistFromCenter)/side_len_ + row_size_/2,
                          (lng < midLng ? -unsignedColDistFromCenter :
                           unsignedColDistFromCenter)/side_len_ + col_size_/2);
}

template <typename FPType, def::ThreadingPolicy policy>
void GridSBSolver<FPType, policy>::fillGrid(std::span<const SBLoc<FPType>> loc_data_span) {
    for (const auto &l : loc_data_span) {
        auto idxPr = getIdx(l.geo_pt[1], l.geo_pt[0]);
        auto &cell = grid_[idxPr.first][idxPr.second];
        numLocs = numLocs - static_cast<int>(cell.erase(&l)) + 1;
        cell.insert(&l);
    }
}

template <typename FPType, def::ThreadingPolicy policy>
void GridSBSolver<FPType, policy>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    findKeyLngLat(loc_data_span);
    constructGrid(loc_data_span);
    fillGrid(loc_data_span);
}

template <typename FPType, def::ThreadingPolicy policy>
void GridSBSolver<FPType, policy>::NNOneCell(const std::unordered_set<const SBLoc<FPType>*> &cell,
FPType lng, FPType lat, FPType &minDist, const SBLoc<FPType>* &best) const {
    for (const auto &l : cell) {
        FPType dist = SBLoc<FPType>::havDist({lat, lng}, l->geo_pt);
        if (dist < minDist) {
            minDist = dist;
            best = l;
        }
    }
}

template <typename FPType, def::ThreadingPolicy policy>
const SBLoc<FPType>* GridSBSolver<FPType, policy>::FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    auto idxPr = getIdx(geo_search_pt[1], geo_search_pt[0]);
    std::size_t r0 = idxPr.first, c0 = idxPr.second;
    FPType minDist = FPType_MAX<FPType>;
    const SBLoc<FPType>* best = nullptr;
    // exact cell check
    NNOneCell(grid_[r0][c0], geo_search_pt[1], geo_search_pt[0], minDist, best);
    
    // Spiral Search
    for (std::size_t d = 1; ; ++d) {
        for (int r = std::max(0,static_cast<int>(r0-d)); r <= std::min(static_cast<int>(r0+d), static_cast<int>(row_size_)-1); ++r) {
            if (c0-d >= 0)
                NNOneCell(grid_[r][c0-d], geo_search_pt[1], geo_search_pt[0], minDist, best);
            if (c0+d < col_size_)
                NNOneCell(grid_[r][c0+d], geo_search_pt[1], geo_search_pt[0], minDist, best);
        }
        for (int c = std::max(0, static_cast<int>(c0-d+1)); c < std::min(static_cast<int>(c0+d), static_cast<int>(col_size_)); c++) {
            if (r0 - d >= 0)
                NNOneCell(grid_[r0-d][c], geo_search_pt[1], geo_search_pt[0], minDist, best);
            if (r0 + d < row_size_)
                NNOneCell(grid_[r0+d][c], geo_search_pt[1], geo_search_pt[0], minDist, best);
        }
        if (minDist < d*side_len_*DISTORT_FACTOR)
            return best;
    }
    return best;
}

template class GridSBSolver<double, def::ThreadingPolicy::kSingle>;
template class GridSBSolver<float, def::ThreadingPolicy::kSingle>;




