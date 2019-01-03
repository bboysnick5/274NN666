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

class BKDTSBSolver : public KDTSBSolver {
public:
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    
protected:
    void generateKDT(const std::shared_ptr<std::vector<SBLoc>>&) override;
};


#endif /* BKDTSBSolver_hpp */
