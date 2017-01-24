//
//  KDTSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef KDTSBSolver_hpp
#define KDTSBSolver_hpp

#include <stdio.h>
#include "SBSolver.hpp"
#include "KDTree.hpp"

class KDTSBSolver : public SBSolver {
public:
    virtual void build();
    const SBLoc* findNearest(double lng, double lat) const;
    
protected:
    KDTree<3, const SBLoc*> kdt;
};



#endif /* KDTSBSolver_hpp */
