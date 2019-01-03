//
//  UpgradeBKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#ifndef UpgradeBKDTGridSBSolver_hpp
#define UpgradeBKDTGridSBSolver_hpp

#include <stdio.h>
#include <vector>
#include <iterator>
#include "GridSBSolver.hpp"
#include "KDTree.hpp"


class UpgradeBKDTGridSBSolver : public GridSBSolver {
public:
    UpgradeBKDTGridSBSolver(double aveLocPerCell = 1);
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lng, double lat) const override;
    
    
private:
    std::vector<KDTree<3, const SBLoc*, DistType::EUC>> gridTreeCache;
    std::vector<const SBLoc*> gridSingleCache;
    double latInc, lngInc;
    KDTree<3, const SBLoc*, DistType::EUC> sbKdt;
    
    void fillGridCache();
    double calcLatIncFromAlpc();
};


#endif /* UpgradeBKDTGridSBSolver_hpp */
