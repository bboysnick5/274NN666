/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-10 00:27:39
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/UniCellBKDTGridSBSolver.hpp
 */

#ifndef UniCellBKDTGridSBSolver_hpp
#define UniCellBKDTGridSBSolver_hpp

#include <stdio.h>

#include <iterator>
#include <memory>
#include <vector>

#include "KDTree.hpp"
#include "UniLatLngBKDTGridSBSolver.hpp"

template <SolverConfig Config>
class UniCellBKDTGridSBSolver final : public UniLatLngBKDTGridSBSolver<Config> {
    using FPType = typename decltype(Config)::FPType;
    using KDTType = typename decltype(Config)::KDTType;

   public:
    UniCellBKDTGridSBSolver(FPType = 1.0, std::size_t = 1500);
    // void Build(std::span<const SBLoc<FPType>>) override final;
    const SBLoc<FPType> *FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    // virtual void PrintSolverInfo() const override;

   private:
    std::vector<std::pair<std::size_t, FPType>> thisRowStartIdx;
    void FillGridCache() override final;
};

#endif /* UniCellBKDTGridSBSolver_hpp */
