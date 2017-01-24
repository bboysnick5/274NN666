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
    const SBLoc* findNearest(double lng, double lat) const;
    
private:
    void build();
};


#endif /* BFSBSolver_hpp */
