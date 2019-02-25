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

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
using KDT = KDTType<dist_type, 3, const SBLoc<dist_type>*, Point<dist_type, 3>::DistType::EUC>;

template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
class KDTSBSolver : public SBSolver<dist_type> {
public:
    void build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
    const SBLoc<dist_type>* findNearest(dist_type lat, dist_type lng) const override;
    virtual void printSolverInfo() const override;
    
protected:
    KDT<KDTType, dist_type> locKdt;
    virtual void generateKDT(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&);
};



#endif /* KDTSBSolver_hpp */
