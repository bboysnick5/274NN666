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


template <typename dist_type>
class BFSBSolver : public SBSolver<dist_type> {
public:
    const SBLoc<dist_type>* findNearest(dist_type, dist_type) const override;
    void build(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &sbData) override;
    void printSolverInfo() const override;
    
    virtual ~BFSBSolver() {}
    
protected:
    std::shared_ptr<std::vector<SBLoc<dist_type>>> locData;
};


#endif /* BFSBSolver_hpp */
