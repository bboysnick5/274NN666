/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-11 11:34:36
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/SBSolver.hpp
 */

#ifndef SBSolver_hpp
#define SBSolver_hpp

#include <stdio.h>

#include <boost/mp11.hpp>
#include <memory>
#include <span>
#include <type_traits>
#include <vector>

#include "Definition.hpp"
#include "KDTree.hpp"
#include "KDTreeCusMem.hpp"
#include "KDTreeExpandLongest.hpp"
#include "KDTreeExpandLongestVec.hpp"
#include "Point.hpp"
#include "SBLoc.hpp"
#include "SolverConfig.hpp"

template <SolverConfig Config>
class SBSolver {
   public:
    using FPType = typename decltype(Config)::FPType;
    using KDTType = typename decltype(Config)::KDTType;

    virtual void Build(std::span<const SBLoc<FPType>>) = 0;
    virtual const SBLoc<FPType>* FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType) const = 0;
    virtual void PrintSolverInfo() const = 0;
    virtual ~SBSolver() = default;
};

template <SolverConfig Config>
concept SBSolverC = requires {
    typename SBSolver<Config>;
};

    // TODO: add a description string to solver with its name and config

#endif /* SBSolver_hpp */
