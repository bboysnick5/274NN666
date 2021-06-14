//
//  BFSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/11/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef BFSBSolver_hpp
#define BFSBSolver_hpp

#include <memory>
#include <vector>
#include "SBSolver.hpp"


template <typename FPType>
class BFSBSolver : public SBSolver<FPType> {
public:
    const SBLoc<FPType>* FindNearestLoc(const PointND<FPType, 2>&) const override;
    void Build(const std::shared_ptr<std::vector<SBLoc<FPType>>> &sbData) override;
    void PrintSolverInfo() const override;
    
    virtual ~BFSBSolver() {}
    
protected:
    std::shared_ptr<std::vector<SBLoc<FPType>>> locData;
};


#endif /* BFSBSolver_hpp */
