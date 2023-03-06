/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-08 20:39:04
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/BFSBSolver.hpp
 */

#ifndef BFSBSolver_hpp
#define BFSBSolver_hpp

#include <concepts>
#include <memory>
#include <span>

#include "Definition.hpp"
#include "SBLoc.hpp"
#include "SBSolver.hpp"

template <SolverConfig Config>
class BFSBSolver : public SBSolver<Config> {
    using FPType = typename decltype(Config)::FPType;

   public:
    virtual void Build(std::span<const SBLoc<FPType>>) override;
    const SBLoc<FPType> *FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    void PrintSolverInfo() const override;

   protected:
    std::span<const SBLoc<FPType>> loc_span_;
    std::vector<typename SBLoc<FPType>::GeoPtType> loc_geo_pt_vec_;
};

#endif /* BFSBSolver_hpp */
