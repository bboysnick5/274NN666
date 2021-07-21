//
//  BFLocalStorageSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 7/9/21.
//  Copyright © 2021 Yunlong Liu. All rights reserved.
//

#ifndef BFLocalStorageSBSolver_hpp
#define BFLocalStorageSBSolver_hpp

#include "Definition.hpp"
#include "SBSolver.hpp"
#include "SBLoc.hpp"

#include <Vc/Vc>

#include <vector>
#include <memory>
#include <span>
#include <concepts>


template <typename FPType, def::ThreadingPolicy policy>
class BFLocalStorageSBSolver : public SBSolver<FPType, policy> {
public:
    virtual void Build(std::span<const SBLoc<FPType>>) override final;
    const SBLoc<FPType>* FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    void PrintSolverInfo() const override;
    virtual ~BFLocalStorageSBSolver() {}
    
protected:
    std::vector<typename SBLoc<FPType>::GeoPtType> geo_pt_vec_;
    std::vector<const SBLoc<FPType>*> loc_address_vec_;
};

#endif /* BFLocalStorageSBSolver_hpp */
