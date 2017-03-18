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
#include "KDTreeVec.hpp"


class BKDTGridSBSolver : public GridSBSolver {
public:
    BKDTGridSBSolver(double aveLocPerCell = 1);
    void build();
    const SBLoc* findNearest(double lng, double lat) const;
    
    
private:
    std::vector<KDTreeVec<3, const SBLoc*, DistType::EUC>> gridTreeCache;
    std::vector<const SBLoc*> gridSingleCache;
    std::vector<std::pair<Point<3>, const SBLoc*>>::iterator
        cacheAllPossibleLocsOneCell(int, int, double,
        std::vector<std::pair<Point<3>, const SBLoc*>>::iterator);
    void fillGridCache();
    double xyzDistFromSideLen();
    
    KDTreeVec<3, const SBLoc*, DistType::EUC> sbKdt;
};


#endif /* BKDTGridSBSolver_hpp */
