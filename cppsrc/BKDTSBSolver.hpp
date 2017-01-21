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
    void build(const std::vector<SBLoc> &sbData);
};


#endif /* BKDTSBSolver_hpp */
