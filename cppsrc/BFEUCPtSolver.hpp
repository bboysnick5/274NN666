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

template <typename dist_type>
class BFEUCPtSBSolver : public BFSBSolver<dist_type> {
public:
    const SBLoc<dist_type>* FindNearestLoc(const Point<dist_type, 2>&) const override;
    void PrintSolverInfo() const override;
    virtual ~BFEUCPtSBSolver() override {}
};

#endif /* BFEUCPtSolver_hpp */
