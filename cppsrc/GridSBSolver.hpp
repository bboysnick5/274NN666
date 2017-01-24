//
//  GridSBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/13/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef GridSBSolver_hpp
#define GridSBSolver_hpp

#include <stdio.h>
#include <unordered_set>
#include <vector>
#include "SBSolver.hpp"

const double DOUBLE_MAX = std::numeric_limits<double>::max();
const double DOUBLE_MIN = std::numeric_limits<double>::lowest();

class GridSBSolver : public SBSolver {
    
public:
    virtual void build();
    virtual const SBLoc* findNearest(double lng, double lat) const;
    GridSBSolver(double aveLocPerCell = 1);
    
protected:
    void findKeyLngLat();
    std::pair<int, int> getIdx(double lng, double lat) const;

    void constructGrid();
    void fillGrid();
    void NNOneCell(const std::unordered_set<const SBLoc*>&,
                  double, double, double&, const SBLoc*&) const;
    
    const double AVE_LOC_PER_CELL;
    std::vector<std::vector<std::unordered_set<const SBLoc*>>> grid;
    double sideLen, minLng, maxLng, minLat, maxLat, midLng, midLat;
    int rowSize, colSize, numLocs = 0;
    
private:
    static constexpr double DISTORT_FACTOR = 0.95;
    
};

#endif /* GridSBSolver_hpp */
