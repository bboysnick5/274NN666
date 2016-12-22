//
//  BKDTSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef BKDTSBSolver_hpp
#define BKDTSBSolver_hpp

#include <stdio.h>
#include "SBSolver.hpp"
#include "KDTree.hpp"

class BKDTSBSolver : public SBSolver {
public:
    void build(const std::vector<SBLoc> &sbData);
    SBLoc findNearest(double lng, double lat) const;
    
private:
    KDTree<3, SBLoc> kdt;
};


#endif /* BKDTSBSolver_hpp */
