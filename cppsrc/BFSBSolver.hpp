//
//  BFSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/11/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef BFSBSolver_hpp
#define BFSBSolver_hpp

#include "SBSolver.hpp"
#include "SBLoc.hpp"
#include <memory>
#include <span>
#include <concepts>


template <typename FPType, def::ThreadingPolicy Policy>
class BFSBSolver : public SBSolver<FPType, Policy> {
public:
    virtual void Build(std::span<const SBLoc<FPType>>) override final;
    const SBLoc<FPType>* FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    void PrintSolverInfo() const override;
    
    virtual ~BFSBSolver() {}
    
protected:
    std::span<const SBLoc<FPType>> loc_data_span_;
};


#endif /* BFSBSolver_hpp */
