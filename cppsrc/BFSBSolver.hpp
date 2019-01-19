//
//  BFSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/11/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef BFSBSolver_hpp
#define BFSBSolver_hpp

#include <stdio.h>
#include <vector>
#include "SBSolver.hpp"



class BFSBSolver : public SBSolver {
public:
    const SBLoc* findNearest(double lng, double lat) const override;
    void build(const std::shared_ptr<std::vector<SBLoc>> &sbData) override;
    void printSolverInfo() const override;
    
    virtual ~BFSBSolver() {}
    
protected:
    std::shared_ptr<std::vector<SBLoc>> locData;
};


#endif /* BFSBSolver_hpp */
