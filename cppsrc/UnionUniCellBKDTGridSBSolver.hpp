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


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
class UnionUniCellBKDTGridSBSolver : public UnionUniLatLngBKDTGridSBSolver<KDTType> {
public:
    UnionUniCellBKDTGridSBSolver(double = 1, size_t = 1500);
    //void build(const std::shared_ptr<std::vector<SBLoc>>&) override final;
    const SBLoc* findNearest(double lng, double lat) const override final;
    //virtual void printSolverInfo() const override;
    
    
private:
    std::vector<std::pair<size_t, double>> thisRowStartIdx;
    void fillGridCache() override final;
};


#endif /* UnionUniCellBKDTGridSBSolver_hpp */
