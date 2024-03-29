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
#include <span>

template <typename FPType>
class SBSolverWrapper {
public:
    virtual void Build(std::span<const SBLoc<FPType>>) = 0;
    virtual const SBLoc<FPType>* FindNearestLoc(typename SBLoc<FPType>::GeoPtType) const = 0;
    virtual void PrintSolverInfo() const = 0;
    virtual ~SBSolverWrapper() = default;
};

template <typename FPType, def::ThreadingPolicy>
class SBSolver : public SBSolverWrapper<FPType> {
public:
    virtual ~SBSolver() {}
};


#endif /* SBSolver_hpp */
