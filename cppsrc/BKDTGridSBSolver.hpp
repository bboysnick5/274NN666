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
    std::vector<std::vector<KDTree<3, const SBLoc*>>> gridTreeCache;
    std::vector<std::vector<const SBLoc*>> gridSingleCache;
    std::vector<const SBLoc*>::iterator spiralSearchAllPossibleLocsOneCell(int,
        int, std::vector<const SBLoc*>&);
    void fillGridCache();
    inline void checkOneCell(const std::unordered_set<const SBLoc*>&, double,
        double, double&, std::vector<const SBLoc*>&, int&) const;
};


#endif /* BKDTGridSBSolver_hpp */
