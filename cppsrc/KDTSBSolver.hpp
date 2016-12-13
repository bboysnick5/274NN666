//
//  KDTSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef KDTSBSolver_hpp
#define KDTSBSolver_hpp

#include <stdio.h>
#include "SBSolver.hpp"
#include "KDTree.hpp"

class KDTSBSolver : public SBSolver {
public:
    void build(const std::vector<SBLoc> &sbData);
    SBLoc findNearest(double lng, double lat);
    
private:
    KDTree<3, SBLoc> kdt;
};



#endif /* KDTSBSolver_hpp */
