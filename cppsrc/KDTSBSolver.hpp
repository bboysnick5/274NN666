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
#include <span>

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
using KDT = KDTType<FPType, 3, const SBLoc<FPType>*, PointND<FPType, 3>::DistType::kEuc>;

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
class KDTSBSolver : public SBSolver<FPType, policy> {
public:
    void Build(std::span<const SBLoc<FPType>>) override;
    const SBLoc<FPType>* FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    virtual void PrintSolverInfo() const override;
    virtual ~KDTSBSolver() override {}


protected:
    KDT<KDTType, FPType> loc_kdt_;
    virtual void GenerateKDT(std::span<const SBLoc<FPType>>);
};



#endif /* KDTSBSolver_hpp */
