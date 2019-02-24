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



template <template <class value_type, size_t, class, typename Point<value_type, 3>::DistType> class KDTType>
class BKDTSBSolver : public KDTSBSolver<KDTType> {
    
public:
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    virtual void printSolverInfo() const override;
    virtual ~BKDTSBSolver() {}
    
protected:
    void generateKDT(const std::shared_ptr<std::vector<SBLoc>>&) override;
};


#endif /* BKDTSBSolver_hpp */
