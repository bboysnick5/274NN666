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
#include "GridSBSolver.hpp"
#include "KDTree.hpp"


class BKDTGridSBSolver : public GridSBSolver {
public:
    void build(const std::vector<SBLoc> &sbData);
    SBLoc findNearest(double lng, double lat) const;
    
private:

    std::vector<std::vector<KDTree<3, SBLoc>>> gridCache;
    
    //double distFromMidPtInCell(double, double, double, double) const;
    void fillCacheOneCell(int r0, int c0);
    void fillGridCache();
    void checkOneCell(const std::unordered_set<SBLoc>&, double, double,
                      double&, std::vector<std::pair<Point<3>, SBLoc>>&) const;
};


#endif /* BKDTGridSBSolver_hpp */
