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


class BKDTGridSBSolver : public GridSBSolver {
public:
    BKDTGridSBSolver(double aveLocPerCell = 1);
    void build();
    const SBLoc* findNearest(double lng, double lat) const;
    
    
private:
    std::vector<KDTree<3, const SBLoc*, DistType::EUC>> gridTreeCache;
    std::vector<const SBLoc*> gridSingleCache;
    std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>>::iterator
        cacheAllPossibleLocsOneCell(size_t, size_t, double,
        std::vector<std::pair<Point<3, DistType::EUC>, const SBLoc*>>::iterator);
    void fillGridCache();
    double xyzDistFromSideLen();
    
    KDTree<3, const SBLoc*, DistType::EUC> sbKdt;
};


#endif /* BKDTGridSBSolver_hpp */
