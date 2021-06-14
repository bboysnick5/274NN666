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
    void Build(const std::shared_ptr<std::vector<SBLoc<FPType>>>&) override;
    virtual void PrintSolverInfo() const override;
    virtual ~BKDTSBSolver() {}
    
protected:
    void GenerateKDT(const std::shared_ptr<std::vector<SBLoc<FPType>>>&) override final;
};


#endif /* BKDTSBSolver_hpp */
