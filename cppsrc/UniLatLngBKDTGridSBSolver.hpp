//
//  UniLatLngBKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#ifndef UniLatLngBKDTGridSBSolver_hpp
#define UniLatLngBKDTGridSBSolver_hpp


#include "BKDTSBSolver.hpp"
#include "KDTree.hpp"
#include <stdio.h>
#include <vector>
#include <iterator>
#include <variant>


template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
class UniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType, FPType, policy> {
public:
    UniLatLngBKDTGridSBSolver(FPType = 1, std::size_t = 1500);
    void Build(std::span<const SBLoc<FPType>>) override final;
    const SBLoc<FPType>* FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    virtual void PrintSolverInfo() const override;
    virtual ~UniLatLngBKDTGridSBSolver() override {}

    
protected:
    const FPType AVE_LOC_PER_CELL;
    const std::size_t kMaxCacheCellVecSize_;
    FPType lng_inc_, lat_inc_, lat_inc_inverse_, side_len_;
    std::size_t totalLocSize, totalNodeSize = 0, singleLocs = 0,
           vecLocs = 0, row_size_, col_size_;
    std::vector<std::variant<std::vector<typename KDT<KDTType, FPType>::node_type>,
                             const SBLoc<FPType>*, KDT<KDTType, FPType>>> grid_cache_;
    void calcSideLenFromAlpc();
    void FillCacheCell(FPType, FPType, FPType,
                       std::vector<typename KDT<KDTType, FPType>::node_type>&);
    const SBLoc<FPType>* ReturnNNLocFromCacheVariant(const typename SBLoc<FPType>::GeoPtType&,
          const std::variant<std::vector<typename KDT<KDTType, FPType>::node_type>,
          const SBLoc<FPType>*, KDT<KDTType, FPType>>&) const;
    
private:
    virtual void FillGridCache();
};

#endif /* UniLatLngBKDTGridSBSolver_hpp */
