//
//  BFEUCPtSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 1/16/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef BFEUCPtSolver_hpp
#define BFEUCPtSolver_hpp

#include "BFSBSolver.hpp"
#include <stdio.h>

class BFEUCPtSBSolver : public BFSBSolver {
public:
    const SBLoc* findNearest(double lat, double lng) const override;
    void printSolverInfo() const override;
};

#endif /* BFEUCPtSolver_hpp */
