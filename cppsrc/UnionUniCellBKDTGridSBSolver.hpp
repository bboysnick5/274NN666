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


template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
class UnionUniCellBKDTGridSBSolver : public UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy> {
public:
    UnionUniCellBKDTGridSBSolver(dist_type = 1.0, std::size_t = 1500);
    const SBLoc<dist_type>* FindNearestLoc(const PointND<dist_type, 2>&) const override final;
    virtual ~UnionUniCellBKDTGridSBSolver() override {}

    
private:
    std::vector<std::pair<std::size_t, dist_type>> thisRowStartIdxThisLngIncInverseVec;
    std::vector<std::pair<std::size_t, dist_type>> colSizeCosLngIncEachRowVec;
    std::size_t totalCacheCells;
    virtual void FillGridCache() override final;
    virtual void LoopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs,
                          typename def::Policy_Tag<def::ThreadingPolicy::kSingle>) override final;
    virtual void LoopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs,
                          typename def::Policy_Tag<def::ThreadingPolicy::kMultiOmp>) override final;
    void ompLoopCollapsePrep();
};


#endif /* UnionUniCellBKDTGridSBSolver_hpp */
