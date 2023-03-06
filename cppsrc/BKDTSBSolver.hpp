/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-09 10:31:33
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/BKDTSBSolver.hpp
 */

#ifndef BKDTSBSolver_hpp
#define BKDTSBSolver_hpp

#include <stdio.h>

#include "KDTSBSolver.hpp"

template <SolverConfig Config>
class BKDTSBSolver : public KDTSBSolver<Config> {
    using FPType = typename decltype(Config)::FPType;
    using KDTType = typename decltype(Config)::KDTType;

   public:
    virtual void Build(std::span<const SBLoc<FPType>>) override;
    virtual void PrintSolverInfo() const override;

   protected:
    virtual void GenerateKDT(std::span<const SBLoc<FPType>>) override final;
};

#endif /* BKDTSBSolver_hpp */
