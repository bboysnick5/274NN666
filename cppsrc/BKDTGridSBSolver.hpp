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

template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type>
class BKDTGridSBSolver : public GridSBSolver<dist_type> {
public:
    BKDTGridSBSolver(dist_type aveLocPerCell = 1);
    void Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
    const SBLoc<dist_type>* FindNearestLoc(const PointND<dist_type, 2>&) const override;
    virtual ~BKDTGridSBSolver() override {}

    
private:
    std::vector<KDTType<dist_type, 3, const SBLoc<dist_type>*, PointND<dist_type, 3>::DistType::EUC>> gridTreeCache;
    std::vector<const SBLoc<dist_type>*> gridSingleCache;
    typename std::vector<typename KDTree<dist_type, 3, const SBLoc<dist_type>*, PointND<dist_type, 3>::DistType::EUC>::node_type>::iterator
        cacheAllPossibleLocsOneCell(std::size_t, std::size_t, dist_type,
                                    typename std::vector<typename KDTree<dist_type, 3, const SBLoc<dist_type>*, PointND<dist_type, 3>::DistType::EUC>::node_type>::iterator);
    void FillGridCache();
    dist_type xyzDistSqFromSideLen();
    
    KDTree<dist_type, 3, const SBLoc<dist_type>*, PointND<dist_type, 3>::DistType::EUC> sbKdt;
};


#endif /* BKDTGridSBSolver_hpp */
