//
//  BKDTSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef BKDTSBSolver_hpp
#define BKDTSBSolver_hpp

#include <stdio.h>
#include "KDTSBSolver.hpp"

template <template <size_t, typename elem, typename Point<3>::DistType> class Tree>
class BKDTSBSolver : public KDTSBSolver<Tree> {
    
public:
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    virtual void printSolverInfo() const override;
    virtual ~BKDTSBSolver() {}
    
protected:
    void generateKDT(const std::shared_ptr<std::vector<SBLoc>>&) override;
};


#endif /* BKDTSBSolver_hpp */
