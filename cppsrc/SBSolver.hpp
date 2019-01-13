//
//  SBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef SBSolver_hpp
#define SBSolver_hpp

#include "SBLoc.hpp"
#include "Point.hpp"
#include <stdio.h>
#include <vector>
#include "memory"


class SBSolver {
public:
    virtual void build(const std::shared_ptr<std::vector<SBLoc>>&) = 0;
    virtual const SBLoc* findNearest(double lng, double lat) const = 0;
    virtual void printSolverInfo() const = 0;
};


#endif /* SBSolver_hpp */
