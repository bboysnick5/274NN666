//
//  BKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/23/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef BKDTGridSBSolver_hpp
#define BKDTGridSBSolver_hpp

#include <stdio.h>
#include <vector>
#include <iterator>
#include "GridSBSolver.hpp"
#include "KDTree.hpp"

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
class BKDTGridSBSolver final : public GridSBSolver<FPType, policy> {
public:
    BKDTGridSBSolver(FPType aveLocPerCell = 1);
    virtual void Build(std::span<const SBLoc<FPType>>) override;
    virtual const SBLoc<FPType>* FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    virtual void PrintSolverInfo() const override {}
    virtual ~BKDTGridSBSolver() override {}

    
private:
    std::vector<KDTType<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::kEuc>> gridTreeCache;
    std::vector<const SBLoc<FPType>*> gridSingleCache;
    typename std::vector<typename KDTree<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::kEuc>::node_type>::iterator
        cacheAllPossibleLocsOneCell(std::size_t, std::size_t, FPType,
                                    typename std::vector<typename KDTree<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::kEuc>::node_type>::iterator);
    void FillGridCache();
    FPType xyzDistSqFromSideLen();
    
    KDTree<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::kEuc> sbKdt;
};


#endif /* BKDTGridSBSolver_hpp */
