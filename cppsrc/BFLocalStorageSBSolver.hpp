/*
 * @Author: Nick Liu
 * @Date: 2021-07-09 13:33:01
 * @LastEditTime: 2022-08-08 17:04:40
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/BFLocalStorageSBSolver.hpp
 */

#ifndef BFLocalStorageSBSolver_hpp
#define BFLocalStorageSBSolver_hpp

#include "Definition.hpp"
#include "SBLoc.hpp"
#include "SBSolver.hpp"

// #include <Vc/Vc>

#include <concepts>
#include <memory>
#include <span>
#include <vector>

template <SolverConfig Config>
class BFLocalStorageSBSolver : public SBSolver<Config> {
    using FPType = typename decltype(Config)::FPType;

   public:
    virtual void Build(std::span<const SBLoc<FPType>>) override final;
    const SBLoc<FPType> *FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    void PrintSolverInfo() const override;

   protected:
    std::vector<typename SBLoc<FPType>::GeoPtType> geo_pt_vec_;
    std::vector<const SBLoc<FPType> *> loc_address_vec_;
};

#endif /* BFLocalStorageSBSolver_hpp */
