/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-08 17:05:23
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/BFCartPtSolver.hpp
 */

#ifndef BFCartPtSolver_hpp
#define BFCartPtSolver_hpp

#include <stdio.h>

#include <vector>

#include "BFSBSolver.hpp"
#include "Point.hpp"
#include "SBLoc.hpp"

template <SolverConfig Config>
class BFCartPtSBSolver final : public BFSBSolver<Config> {
    using FPType = typename decltype(Config)::FPType;

   public:
    virtual void Build(std::span<const SBLoc<FPType>>) override;
    const SBLoc<FPType> *FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    void PrintSolverInfo() const override;

   private:
    std::vector<typename SBLoc<FPType>::CartPtType> loc_cart_pt_vec_;
};

#endif /* BFCartPtSolver_hpp */
