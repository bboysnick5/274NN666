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


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
class UnionUniCellBKDTGridSBSolver : public UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type> {
public:
    UnionUniCellBKDTGridSBSolver(dist_type = 1, size_t = 1500);
    const SBLoc<dist_type>* findNearest(dist_type lat, dist_type lng) const override final;
    
    
private:
    std::vector<std::pair<size_t, dist_type>> thisRowStartIdx;
    void fillGridCache() override final;
};


#endif /* UnionUniCellBKDTGridSBSolver_hpp */
