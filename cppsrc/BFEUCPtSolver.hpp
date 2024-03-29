//
//  BFEUCPtSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 1/16/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef BFEUCPtSolver_hpp
#define BFEUCPtSolver_hpp

#include "Point.hpp"
#include "SBLoc.hpp"
#include "BFSBSolver.hpp"
#include <stdio.h>

template <typename FPType, def::ThreadingPolicy Policy>
class BFEUCPtSBSolver final : public BFSBSolver<FPType, Policy> {
public:
    const SBLoc<FPType>* FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    void PrintSolverInfo() const override;
    virtual ~BFEUCPtSBSolver() override {}
};

#endif /* BFEUCPtSolver_hpp */
