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



template <template <class DT, std::size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
class BKDTSBSolver : public KDTSBSolver<KDTType, dist_type> {
    
public:
    void Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
    virtual void PrintSolverInfo() const override;
    virtual ~BKDTSBSolver() {}
    
protected:
    void GenerateKDT(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override final;
};


#endif /* BKDTSBSolver_hpp */
