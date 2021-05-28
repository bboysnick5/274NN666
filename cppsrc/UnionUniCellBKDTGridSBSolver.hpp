//
//  UnionUnionUniCellBKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 2/20/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef UnionUniCellBKDTGridSBSolver_hpp
#define UnionUniCellBKDTGridSBSolver_hpp

#include "UnionUniLatLngBKDTGridSBSolver.hpp"
#include "KDTree.hpp"
#include <stdio.h>
#include <vector>
#include <memory>
#include <iterator>


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType,
          class dist_type, Def::Threading_Policy policy>
class UnionUniCellBKDTGridSBSolver : public UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy> {
public:
    UnionUniCellBKDTGridSBSolver(dist_type = 1.0, size_t = 1500);
    const SBLoc<dist_type>* findNearest(const Point<dist_type, 2>&) const override final;
    
    
private:
    std::vector<std::pair<size_t, dist_type>> thisRowStartIdxThisLngIncInverseVec;
    std::vector<std::size_t> colSizeInEachRowVec;
    std::size_t totalCacheCells;
    virtual void fillGridCache() override final;
    virtual void loopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs,
                          typename UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::template Policy_Tag<Def::Threading_Policy::SINGLE>) override final;
    virtual void loopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs,
                          typename UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::template Policy_Tag<Def::Threading_Policy::MULTI_OMP>) override final;
    void ompLoopCollapsePrep();
};


#endif /* UnionUniCellBKDTGridSBSolver_hpp */
