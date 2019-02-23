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
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lat, double lng) const override;
    
    
private:
    std::vector<KDTree<3, const SBLoc*, Point<3>::DistType::EUC>> gridTreeCache;
    std::vector<const SBLoc*> gridSingleCache;
    std::vector<std::pair<Point<3>, const SBLoc*>>::iterator
        cacheAllPossibleLocsOneCell(size_t, size_t, double,
        std::vector<std::pair<Point<3>, const SBLoc*>>::iterator);
    void fillGridCache();
    double xyzDistFromSideLen();
    
    KDTree<3, const SBLoc*, Point<3>::DistType::EUC> sbKdt;
};


#endif /* BKDTGridSBSolver_hpp */
