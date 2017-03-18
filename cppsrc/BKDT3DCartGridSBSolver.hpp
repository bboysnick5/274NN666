//
//  BKDT3DCartGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 1/27/17.
//  Copyright © 2017 Yunlong Liu. All rights reserved.
//

#ifndef BKDT3DCartGridSBSolver_hpp
#define BKDT3DCartGridSBSolver_hpp

#include <stdio.h>
#include "SBSolver.hpp"
#include "KDTree.hpp"

class BKDT3DCartGridSBSolver : public SBSolver {
public:
    void build();
    const SBLoc* findNearest(double lng, double lat) const;
    BKDT3DCartGridSBSolver(double aveLocsPerCell = 1);
    
private:
    size_t getIdx(size_t x, size_t y, size_t z) const;
    size_t getIdx(double x, double y, double z) const;
    
    std::vector<KDTree<3, const SBLoc*, DistType::EUC>> gridCache;
    KDTree<3, const SBLoc*, DistType::EUC> sbKdt;
    size_t gridXSize, gridYSize, gridZSize;
    double minX, minY, minZ, maxX, maxY, maxZ, midX, midY, midZ, sideLen;
    const double AVE_LOC_PER_CELL;
};


#endif /* BKDT3DCartGridSBSolver_hpp */
