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
    void build(const std::vector<SBLoc> &sbData);
    SBLoc findNearest(double lng, double lat) const;
    
private:
    std::vector<SBLoc> sbData;
};


#endif /* BFSBSolver_hpp */
