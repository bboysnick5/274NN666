/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-08-09 10:35:25
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/UniLatLngBKDTGridSBSolver.hpp
 */

#ifndef UniLatLngBKDTGridSBSolver_hpp
#define UniLatLngBKDTGridSBSolver_hpp

#include <stdio.h>

#include <iterator>
#include <variant>
#include <vector>

#include "BKDTSBSolver.hpp"
#include "KDTree.hpp"

template <SolverConfig Config>
class UniLatLngBKDTGridSBSolver : public BKDTSBSolver<Config> {
    using FPType = typename decltype(Config)::FPType;
    using KDTType = typename decltype(Config)::KDTType;

   public:
    UniLatLngBKDTGridSBSolver(FPType = 1, std::size_t = 1500);
    void Build(std::span<const SBLoc<FPType>>) override final;
    const SBLoc<FPType> *FindNearestLoc(
        typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    virtual void PrintSolverInfo() const override;

   protected:
    const FPType AVE_LOC_PER_CELL;
    const std::size_t kMaxCacheCellVecSize_;
    FPType lng_inc_, lat_inc_, lat_inc_inverse_, side_len_;
    std::size_t totalLocSize, totalNodeSize = 0, singleLocs = 0, vecLocs = 0,
                              row_size_, col_size_;
    std::vector<std::variant<std::vector<typename KDTType::node_type>,
                             const SBLoc<FPType> *, KDTType>>
        grid_cache_;
    void calcSideLenFromAlpc();
    void FillCacheCell(FPType, FPType, FPType,
                       std::vector<typename KDTType::node_type> &);
    const SBLoc<FPType> *ReturnNNLocFromCacheVariant(
        const typename SBLoc<FPType>::GeoPtType &,
        const std::variant<std::vector<typename KDTType::node_type>,
                           const SBLoc<FPType> *, KDTType> &) const;

   private:
    virtual void FillGridCache();
};

#endif /* UniLatLngBKDTGridSBSolver_hpp */
