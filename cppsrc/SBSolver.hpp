//
//  SBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef SBSolver_hpp
#define SBSolver_hpp

#include <stdio.h>
#include <vector>
#include "memory"
#include "Point.hpp"
#include "SBLoc.hpp"


class SBSolver {
public:
    inline void build(const std::shared_ptr<std::vector<SBLoc>> &sbData);
    virtual const SBLoc* findNearest(double lng, double lat) const = 0;
protected:
    std::shared_ptr<std::vector<SBLoc>> sbData;
private:
    virtual void build() = 0;
};

inline void SBSolver::build(const std::shared_ptr<std::vector<SBLoc>> &sbData) {
    this->sbData = sbData;
    build();
}

#endif /* SBSolver_hpp */
