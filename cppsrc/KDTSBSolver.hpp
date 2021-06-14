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

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
using KDT = KDTType<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::EUC>;

template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
class KDTSBSolver : public SBSolver<FPType> {
public:
    void Build(const std::shared_ptr<std::vector<SBLoc<FPType>>>&) override;
    const SBLoc<FPType>* FindNearestLoc(const PointND<FPType, 2>&) const override;
    virtual void PrintSolverInfo() const override;
    virtual ~KDTSBSolver() override {}


protected:
    KDT<KDTType, FPType> locKdt;
    virtual void GenerateKDT(const std::shared_ptr<std::vector<SBLoc<FPType>>>&);
};



#endif /* KDTSBSolver_hpp */
