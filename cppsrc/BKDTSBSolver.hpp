//
//  BKDTSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/12/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef BKDTSBSolver_hpp
#define BKDTSBSolver_hpp

#include <stdio.h>
#include "KDTSBSolver.hpp"



template <template <typename FPType, std::size_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
class BKDTSBSolver : public KDTSBSolver<KDTType, FPType> {
    
public:
    virtual void Build(std::span<const SBLoc<FPType>>) override;
    virtual void PrintSolverInfo() const override;
    virtual ~BKDTSBSolver() {}
    
protected:
    virtual void GenerateKDT(std::span<const SBLoc<FPType>>) override final;
};


#endif /* BKDTSBSolver_hpp */
