//
//  UnionUnionUniCellBKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 2/20/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef UnionUniCellBKDTGridSBSolver_hpp
#define UnionUniCellBKDTGridSBSolver_hpp

#include "UnionUniLatLngBKDTGridSBSolver.hpp"
#include "KDTree.hpp"
#include <stdio.h>
#include <vector>
#include <memory>
#include <iterator>


template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy Policy>
class UnionUniCellBKDTGridSBSolver final : public UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy> {
public:
    UnionUniCellBKDTGridSBSolver(FPType = 1.0, std::size_t = 1500);
    const SBLoc<FPType>* FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    virtual ~UnionUniCellBKDTGridSBSolver() override {}

    
private:
    std::vector<std::pair<std::size_t, FPType>> thisRowStartIdxThisLngIncInverseVec;
    std::vector<std::pair<std::size_t, FPType>> col_size_CosLngIncEachRowVec;
    std::size_t totalCacheCells;
    virtual void FillGridCache() override final;
    virtual void LoopBody(typename def::ThreadingPolicyTag<def::ThreadingPolicy::kSingle>) override;
    virtual void LoopBody(typename def::ThreadingPolicyTag<def::ThreadingPolicy::kMultiOmp>) override;
};


#endif /* UnionUniCellBKDTGridSBSolver_hpp */
