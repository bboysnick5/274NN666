/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-08 20:54:30
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/KDTSBSolver.hpp
 */

#ifndef KDTSBSolver_hpp
#define KDTSBSolver_hpp

#include <stdio.h>

#include <span>

#include "KDTree.hpp"
#include "KDTreeCusMem.hpp"
#include "KDTreeExpandLongest.hpp"
#include "KDTreeExpandLongestVec.hpp"
#include "SBSolver.hpp"

template <SolverConfig Config>
class KDTSBSolver : public SBSolver<Config> {
    using FPType = typename decltype(Config)::FPType;
    using KDTType = typename decltype(Config)::KDTType;

   public:
    void Build(std::span<const SBLoc<FPType>>) override;
    const SBLoc<FPType> *FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    virtual void PrintSolverInfo() const override;

   protected:
    KDTType loc_kdt_;
    virtual void GenerateKDT(std::span<const SBLoc<FPType>>);
};

#endif /* KDTSBSolver_hpp */
