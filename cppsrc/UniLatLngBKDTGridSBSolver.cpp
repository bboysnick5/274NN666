//
//  UniLatLngBKDTGridSBSolver.cpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#include "UniLatLngBKDTGridSBSolver.hpp"
#include "Utility.hpp"
//#include <omp.h>

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
UniLatLngBKDTGridSBSolver<KDTType, FPType>::
UniLatLngBKDTGridSBSolver(FPType alpc, std::size_t maxCacheCellVecSize)
: BKDTSBSolver<KDTType, FPType>(), AVE_LOC_PER_CELL(alpc),
  kMaxCacheCellVecSize_(maxCacheCellVecSize) {}


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void UniLatLngBKDTGridSBSolver<KDTType, FPType>::PrintSolverInfo() const {
    std::cout << "Total cache locs: " << totalNodeSize
    << "\nAve tree size: " << totalNodeSize/grid_cache_.size()
    << "\nAve tree height: "
    << static_cast<std::size_t>(log2(totalNodeSize/grid_cache_.size() + 1)) + 1
    << "\nRatio of cache locs over actual num locs: "
    << totalNodeSize/totalLocSize
    << "\nTotal num of loc cells: " << grid_cache_.size()
    << "\nSingle loc cells: " << singleLocs
    << "\nVector loc cells: " << vecLocs
    << "\nKd-tree loc cells: " << grid_cache_.size() -singleLocs -vecLocs << "\n";
}


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void UniLatLngBKDTGridSBSolver<KDTType, FPType>::
FillCacheCell(FPType thisCtrLng, FPType thisCtrLat, FPType diagonalDistSq3DEUC,
              std::vector<typename KDT<KDTType, FPType>::node_type>& pt_loc_vec) {
    this->loc_kdt_.NNsWithFence(SBLoc<FPType>::GeoPtTo3dEucPt({thisCtrLat, thisCtrLng}),
                                   diagonalDistSq3DEUC, std::back_inserter(pt_loc_vec));
    std::size_t locsSize = pt_loc_vec.size();
    this->totalNodeSize += locsSize;
    if (locsSize == 1) {
        this->grid_cache_.emplace_back(pt_loc_vec[0].value);
        this->singleLocs++;
    } else if (locsSize < kMaxCacheCellVecSize_) {
        this->grid_cache_.emplace_back(pt_loc_vec);
        this->vecLocs++;
    } else {
        this->grid_cache_.emplace_back(std::in_place_type<KDT<KDTType, FPType>>,
                                     pt_loc_vec.begin(), pt_loc_vec.end());
    }
    pt_loc_vec.clear();
}


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void UniLatLngBKDTGridSBSolver<KDTType, FPType>::FillGridCache() {
    col_size_ = row_size_;
    lng_inc_ = 2.0*def::kMathPi<FPType>/col_size_ + 2.0*def::kMathPi<FPType>/(col_size_*65536);
    grid_cache_.reserve(this->loc_kdt_.size()*1.2/AVE_LOC_PER_CELL);
    std::vector<typename KDT<KDTType, FPType>::node_type> pt_loc_vec;
    pt_loc_vec.reserve(kMaxCacheCellVecSize_);
    
    FPType thisCtrLat = 0.5 * (lat_inc_ - def::kMathPi<FPType>);
    for (std::size_t r = 0; r < row_size_; ++r, thisCtrLat += lat_inc_) {
        FPType thisCtrLng = 0.5 * lng_inc_ - def::kMathPi<FPType>;
        FPType lat1 = r*this->lat_inc_ - 0.5*def::kMathPi<FPType>;
        FPType diagonalDistSq3DEUC = SBLoc<FPType>::EUC3DDistSqFromLatDeltaLng(lat1, lat1 + lat_inc_, lng_inc_);
        for (std::size_t c = 0; c < col_size_; ++c, thisCtrLng += lng_inc_) {
            FillCacheCell(thisCtrLng, thisCtrLat, diagonalDistSq3DEUC, pt_loc_vec);
        }
    }
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void UniLatLngBKDTGridSBSolver<KDTType, FPType>::calcSideLenFromAlpc() {
    FPType surfaceArea = 4*def::kMathPi<FPType>*SBLoc<FPType>::EARTH_RADIUS*SBLoc<FPType>::EARTH_RADIUS;
    FPType numCells = this->loc_kdt_.size()/AVE_LOC_PER_CELL;
    side_len_ = sqrt(surfaceArea/numCells);
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
void UniLatLngBKDTGridSBSolver<KDTType, FPType>::Build(std::span<const SBLoc<FPType>> loc_data_span) {
    totalLocSize = loc_data_span.size();
    BKDTSBSolver<KDTType, FPType>::GenerateKDT(loc_data_span);
    calcSideLenFromAlpc();
    lat_inc_ = std::fabs(SBLoc<FPType>::deltaLatOnSameLngFromHavDist(side_len_));
    lat_inc_inverse_ = 1.0/lat_inc_;
    row_size_ = std::ceil(def::kMathPi<FPType>/lat_inc_);
    FillGridCache();
    this->loc_kdt_.Clear();
}

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
const SBLoc<FPType>* UniLatLngBKDTGridSBSolver<KDTType, FPType>::
ReturnNNLocFromCacheVariant(const PointND<FPType, 2>& geo_search_pt,
                            const std::variant<std::vector<typename KDT<KDTType, FPType>::node_type>, const SBLoc<FPType>*, KDT<KDTType, FPType>>& v) const {
    switch (v.index()) {
        case 0: {
            const auto p = SBLoc<FPType>::GeoPtTo3dEucPt(geo_search_pt);
            const auto &vec = std::get<0>(v);
            return Utility::MinElementGivenDistFunc(vec.cbegin(), vec.cend(),
                                               [&p](const auto& nh){return p.template
                                                    dist<PointND<FPType, 3>::DistType::EUCSQ>(nh.key);}, std::less()
                                               )->value;
        }
        case 1:
            return std::get<1>(v);
        default:
            return std::get<2>(v).kNNValue(SBLoc<FPType>::GeoPtTo3dEucPt(geo_search_pt), 1);
    }
}


template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
const SBLoc<FPType>* UniLatLngBKDTGridSBSolver<KDTType, FPType>::
FindNearestLoc(const PointND<FPType, 2>& geo_search_pt) const {
    return ReturnNNLocFromCacheVariant(geo_search_pt, grid_cache_[static_cast<std::size_t>
    ((geo_search_pt[0]+0.5*def::kMathPi<FPType>)*lat_inc_inverse_)*col_size_+ static_cast<std::size_t>((geo_search_pt[1]+def::kMathPi<FPType>)/lng_inc_)]);
}



template class UniLatLngBKDTGridSBSolver<KDTree, double>;
template class UniLatLngBKDTGridSBSolver<KDTree, float>;

template class UniLatLngBKDTGridSBSolver<KDTreeCusMem, double>;
template class UniLatLngBKDTGridSBSolver<KDTreeCusMem, float>;

template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongest, double>;
template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongest, float>;

template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, double>;
template class UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, float>;

