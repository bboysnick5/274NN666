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

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type>
using KDT = KDTType<dist_type, 3, const SBLoc<dist_type>*, PointND<dist_type, 3>::DistType::EUC>;

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type>
class KDTSBSolver : public SBSolver<dist_type> {
public:
    void Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
    const SBLoc<dist_type>* FindNearestLoc(const PointND<dist_type, 2>&) const override;
    virtual void PrintSolverInfo() const override;
    virtual ~KDTSBSolver() override {}


protected:
    KDT<KDTType, dist_type> locKdt;
    virtual void GenerateKDT(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&);
};



#endif /* KDTSBSolver_hpp */
