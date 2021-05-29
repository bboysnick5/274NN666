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
#include <memory>

template <typename dist_type>
class SBSolver {
public:
    virtual void build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) = 0;
    virtual const SBLoc<dist_type>* findNearest(const Point<dist_type, 2>&) const = 0;
    virtual void printSolverInfo() const = 0;
    virtual ~SBSolver() {}
};


#endif /* SBSolver_hpp */
