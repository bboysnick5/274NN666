/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-09 10:55:41
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/UnionUniCellBKDTGridSBSolver.hpp
 */
//
//  UnionUnionUniCellBKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 2/20/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef UnionUniCellBKDTGridSBSolver_hpp
#define UnionUniCellBKDTGridSBSolver_hpp

#include <stdio.h>

#include <iterator>
#include <memory>
#include <vector>

#include "KDTree.hpp"
#include "UnionUniLatLngBKDTGridSBSolver.hpp"

template <SolverConfig Config>
class UnionUniCellBKDTGridSBSolver final
    : public UnionUniLatLngBKDTGridSBSolver<Config> {
    using FPType = typename decltype(Config)::FPType;
    using KDTType = typename decltype(Config)::KDTType;

   public:
    UnionUniCellBKDTGridSBSolver(FPType = 1.0, std::size_t = 1500);
    const SBLoc<FPType> *FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;

   private:
    std::vector<std::pair<std::size_t, FPType>>
        thisRowStartIdxThisLngIncInverseVec;
    std::vector<std::pair<std::size_t, FPType>> col_size_CosLngIncEachRowVec;
    std::size_t totalCacheCells;
    virtual void FillGridCache() override final;
    virtual void LoopBody(
        typename def::ThreadingPolicyTag<def::ThreadingPolicy::kSingle>)
        override;
    virtual void LoopBody(
        typename def::ThreadingPolicyTag<def::ThreadingPolicy::kMultiOmp>)
        override;
};

#endif /* UnionUniCellBKDTGridSBSolver_hpp */
