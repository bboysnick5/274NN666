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

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
class BKDTGridSBSolver : public GridSBSolver<FPType> {
public:
    BKDTGridSBSolver(FPType aveLocPerCell = 1);
    void Build(const std::shared_ptr<std::vector<SBLoc<FPType>>>&) override;
    const SBLoc<FPType>* FindNearestLoc(const PointND<FPType, 2>&) const override;
    virtual ~BKDTGridSBSolver() override {}

    
private:
    std::vector<KDTType<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::EUC>> gridTreeCache;
    std::vector<const SBLoc<FPType>*> gridSingleCache;
    typename std::vector<typename KDTree<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::EUC>::node_type>::iterator
        cacheAllPossibleLocsOneCell(std::size_t, std::size_t, FPType,
                                    typename std::vector<typename KDTree<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::EUC>::node_type>::iterator);
    void FillGridCache();
    FPType xyzDistSqFromSideLen();
    
    KDTree<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::EUC> sbKdt;
};


#endif /* BKDTGridSBSolver_hpp */
