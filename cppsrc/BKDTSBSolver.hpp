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



template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
class BKDTSBSolver : public KDTSBSolver<KDTType, FPType, Policy> {
    
public:
    virtual void Build(std::span<const SBLoc<FPType>>) override;
    virtual void PrintSolverInfo() const override;
    virtual ~BKDTSBSolver() {}
    
protected:
    virtual void GenerateKDT(std::span<const SBLoc<FPType>>) override final;
};


#endif /* BKDTSBSolver_hpp */
