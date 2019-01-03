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
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lng, double lat) const override;
    
protected:
    KDTree<3, const SBLoc*, DistType::EUC> locKdt;
    
    virtual void generateKDT(const std::shared_ptr<std::vector<SBLoc>>&);
};



#endif /* KDTSBSolver_hpp */
