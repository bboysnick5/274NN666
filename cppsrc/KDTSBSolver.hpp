//
//  KDTSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef KDTSBSolver_hpp
#define KDTSBSolver_hpp

#include "SBSolver.hpp"
#include "KDTree.hpp"
#include "KDTreeCusMem.hpp"
#include "KDTreeExpandLongest.hpp"
#include <stdio.h>

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
using KDT = KDTType<3, const SBLoc*, Point<3>::DistType::EUC>;

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
class KDTSBSolver : public SBSolver {
public:
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lng, double lat) const override;
    virtual void printSolverInfo() const override;
    
protected:
    KDT<KDTType> locKdt;
    virtual void generateKDT(const std::shared_ptr<std::vector<SBLoc>>&);
};



#endif /* KDTSBSolver_hpp */
