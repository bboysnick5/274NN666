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



template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
class BKDTSBSolver : public KDTSBSolver<KDTType, dist_type> {
    
public:
    void build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
    virtual void printSolverInfo() const override;
    virtual ~BKDTSBSolver() {}
    
protected:
    void generateKDT(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
};


#endif /* BKDTSBSolver_hpp */
