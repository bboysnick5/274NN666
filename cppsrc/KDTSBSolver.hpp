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
#include "KDTreeExpandLongestVec.hpp"
#include <stdio.h>


template <template <class value_type, size_t, class, typename Point<value_type, 3>::DistType> class KDTType>
using KDT = KDTType<double, 3, const SBLoc*, Point<double, 3>::DistType::EUC>;

template <template <class value_type, size_t, class, typename Point<value_type, 3>::DistType> class KDTType>
class KDTSBSolver : public SBSolver {
public:
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lat, double lng) const override;
    virtual void printSolverInfo() const override;
    
protected:
    KDT<KDTType> locKdt;
    virtual void generateKDT(const std::shared_ptr<std::vector<SBLoc>>&);
};



#endif /* KDTSBSolver_hpp */
