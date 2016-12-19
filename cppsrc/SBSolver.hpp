//
//  SBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef SBSolver_hpp
#define SBSolver_hpp

#include <stdio.h>
#include <vector>
#include "Point.hpp"
#include "SBLoc.hpp"


class SBSolver {
public:
    virtual void build(const std::vector<SBLoc> &sbData) = 0;
    virtual SBLoc findNearest(double lng, double lat) = 0;
};




#endif /* SBSolver_hpp */
