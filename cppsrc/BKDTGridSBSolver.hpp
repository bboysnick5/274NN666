/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-08 20:49:47
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/BKDTGridSBSolver.hpp
 */

#ifndef BKDTGridSBSolver_hpp
#define BKDTGridSBSolver_hpp

#include <stdio.h>

#include <iterator>
#include <vector>

#include "GridSBSolver.hpp"
#include "KDTree.hpp"

template <SolverConfig Config>
class BKDTGridSBSolver final : public GridSBSolver<Config> {
    using FPType = typename decltype(Config)::FPType;
    using KDTType = typename decltype(Config)::KDTType;

   public:
    BKDTGridSBSolver(FPType aveLocPerCell = 1);
    virtual void Build(std::span<const SBLoc<FPType>>) override;
    virtual const SBLoc<FPType> *FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    virtual void PrintSolverInfo() const override {}

   private:
    std::vector<KDTType> gridTreeCache;
    std::vector<const SBLoc<FPType> *> gridSingleCache;
    typename std::vector<typename KDTType::node_type>::iterator
        cacheAllPossibleLocsOneCell(
            std::size_t, std::size_t, FPType,
            typename std::vector<
                typename KDTree<FPType, 3, const SBLoc<FPType> *,
                                def::DistType::kEuc>::node_type>::iterator);
    void FillGridCache();
    FPType xyzDistSqFromSideLen();

    KDTree<FPType, 3, const SBLoc<FPType> *, def::DistType::kEuc> sbKdt;
};

#endif /* BKDTGridSBSolver_hpp */
