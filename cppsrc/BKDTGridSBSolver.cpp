//
//  BKDTGridSBSolver<KDTType, FPType, Policy>.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/23/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

/*
 Build Version 1
 Inaccurate and waste of resource.
 Precalculate the side length in HavDist. Then apply that to each cell.
 */

#include <thread>
// #include <omp.h>
#include "BKDTGridSBSolver.hpp"
#include "SBSolver.hpp"

template <SolverConfig Config>
BKDTGridSBSolver<Config>::BKDTGridSBSolver(FPType alpc)
    : GridSBSolver<Config>(alpc) {}

template <SolverConfig Config>
typename std::vector<typename decltype(Config)::KDTType::node_type>::iterator
BKDTGridSBSolver<Config>::cacheAllPossibleLocsOneCell(
    std::size_t r0, std::size_t c0, FPType diffSq,
    typename std::vector<
        typename KDTree<FPType, 3, const SBLoc<FPType> *,
                        def::DistType::kEuc>::node_type>::iterator begin) {
    FPType rowDistCellCtrGridCtr =
               ((r0 + 0.5 - this->row_size_ / 2) * this->side_len_),
           colDistCellCtrGridCtr =
               ((c0 + 0.5 - this->col_size_ / 2) * this->side_len_);
    FPType cellCtrLat =
               this->midLat + SBLoc<FPType>::deltaLatOnSameLngFromHavDist(
                                  rowDistCellCtrGridCtr),
           cellCtrLng = SBLoc<FPType>::lngFromSameLatHavDist(
               colDistCellCtrGridCtr, this->midLng, cellCtrLat);
    return sbKdt.NNsWithFence(
        SBLoc<FPType>::GeoPtToCartPt({cellCtrLat, cellCtrLng}), diffSq, begin);
}

template <SolverConfig Config>
void BKDTGridSBSolver<Config>::FillGridCache() {
    gridTreeCache = std::vector<
        KDTree<FPType, 3, const SBLoc<FPType> *, def::DistType::kEuc>>(
        this->row_size_ * this->col_size_);
    gridSingleCache =
        std::vector<const SBLoc<FPType> *>(this->row_size_ * this->col_size_);
    std::vector<typename KDTree<FPType, 3, const SBLoc<FPType> *,
                                def::DistType::kEuc>::node_type>
        pt_loc_vec(this->numLocs);
    std::size_t totalTreeSize = 0, singleLocs = 0;
    FPType diffSq = xyzDistSqFromSideLen();
    //#pragma omp parallel for num_threads(std::thread::hardware_concurrency()) \
//default(none) schedule(guided) shared(diffSq) firstprivate(pt_loc_vec) \
//reduction(+:totalTreeSize, singleLocs) collapse(2)
    for (std::size_t r = 0; r < this->row_size_; ++r) {
        for (std::size_t c = 0; c < this->col_size_; ++c) {
            std::size_t idx = r * this->col_size_ + c;
            auto locsEnd =
                cacheAllPossibleLocsOneCell(r, c, diffSq, pt_loc_vec.begin());
            std::size_t locsSize = locsEnd - pt_loc_vec.begin();
            if (locsSize > 1) {
                gridTreeCache[idx] =
                    KDTree<FPType, 3, const SBLoc<FPType> *,
                           def::DistType::kEuc>(pt_loc_vec.begin(), locsEnd);
            } else {
                gridSingleCache[idx] = pt_loc_vec[0].value;
                singleLocs++;
            }
            totalTreeSize += locsSize;
        }
    }
    std::size_t multiLocs = this->row_size_ * this->col_size_ - singleLocs;
    std::cout << "ave tree size: "
              << totalTreeSize / (this->row_size_ * this->col_size_)
              << "\nSingle loc cells: " << singleLocs
              << "\nMulti-loc cells:" << multiLocs << std::endl;
}

template <SolverConfig Config>
typename decltype(Config)::FPType
BKDTGridSBSolver<Config>::xyzDistSqFromSideLen() {
    FPType lat2 =
        SBLoc<FPType>::deltaLatOnSameLngFromHavDist(this->side_len_ * sqrt(2));
    return SBLoc<FPType>::CartPtType::template dist<def::DistType::kEucSq>(
        SBLoc<FPType>::GeoPtToCartPt({0.0, 0.0}),
        SBLoc<FPType>::GeoPtToCartPt({lat2, 0.0}));
}

template <SolverConfig Config>
void BKDTGridSBSolver<Config>::Build(
    std::span<const SBLoc<FPType>> loc_data_span) {
    std::vector<typename KDTree<FPType, 3, const SBLoc<FPType> *,
                                def::DistType::kEuc>::node_type>
        kdt_data_vec;
    kdt_data_vec.reserve(loc_data_span.size());
    std::transform(loc_data_span.rbegin(), loc_data_span.rend(),
                   std::back_inserter(kdt_data_vec),
                   [&](const SBLoc<FPType> &loc) ->
                   typename KDTree<FPType, 3, const SBLoc<FPType> *,
                                   def::DistType::kEuc>::node_type {
                       return {SBLoc<FPType>::GeoPtToCartPt(loc.geo_pt), &loc};
                   });
    sbKdt = KDTree<FPType, 3, const SBLoc<FPType> *, def::DistType::kEuc>(
        kdt_data_vec.begin(), kdt_data_vec.end());
    this->numLocs = sbKdt.size();
    this->findKeyLngLat(loc_data_span);
    this->row_size_ = sqrt(loc_data_span.size()) / this->AVE_LOC_PER_CELL;
    this->side_len_ =
        SBLoc<FPType>::havDist({this->minLat, 0.0}, {this->maxLat, 0.0}) /
        this->row_size_;
    FPType lowestLatCircleRadius =
        SBLoc<FPType>::EARTH_RADIUS *
        std::cos(this->minLat < 0 ? 0 : this->minLat);
    FPType longestColDistSpan =
        2 * def::kMathPi<FPType> * lowestLatCircleRadius *
        (std::fabs(this->maxLng - this->minLng) / (2 * def::kMathPi<FPType>));
    this->col_size_ = longestColDistSpan / this->side_len_ + 1;

    FillGridCache();
}

template <SolverConfig Config>
const SBLoc<typename decltype(Config)::FPType>
    *BKDTGridSBSolver<Config>::FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const {
    auto idxPr = this->getIdx(geo_search_pt[1], geo_search_pt[0]);
    std::size_t idx = idxPr.first * this->col_size_ + idxPr.second;
    auto singleLoc = gridSingleCache[idx];
    return singleLoc != nullptr
               ? singleLoc
               : gridTreeCache[idx].kNNValue(
                     SBLoc<FPType>::GeoPtToCartPt(geo_search_pt), 1);
}

template class BKDTGridSBSolver<kConfigStSoaNoSimdKDTree<double>>;
template class BKDTGridSBSolver<kConfigStSoaNoSimdKDTree<float>>;

/*
template class BKDTGridSBSolver<KDTreeCusMem, double>;
template class BKDTGridSBSolver<KDTreeCusMem, float>;

template class BKDTGridSBSolver<KDTreeExpandLongest, double>;
template class BKDTGridSBSolver<KDTreeExpandLongest, float>;

template class BKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class BKDTGridSBSolver<KDTreeExpandLongestVec, float>;
*/
