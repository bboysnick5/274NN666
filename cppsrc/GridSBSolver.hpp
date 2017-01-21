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

class GridSBSolver : public SBSolver {
    
public:
    virtual void build(const std::vector<SBLoc> &sbData);
    virtual SBLoc findNearest(double lng, double lat) const;
    
    
protected:
    static constexpr double DISTORT_FACTOR = 0.9;
    
    static std::vector<double> findMaxLngLat(const std::vector<SBLoc>&);
    std::pair<int, int> getIdx(double lng, double lat) const;

    void constructGrid(const std::vector<SBLoc>&, const std::vector<double>&);
    void fillGrid(const std::vector<SBLoc>&);
    void NNOneCell(const std::unordered_set<SBLoc>&,
                  double, double, double&, SBLoc&) const;
    
    std::vector<std::vector<std::unordered_set<SBLoc>>> grid;
    double sideLen, midLng, midLat;
    int rowSize, colSize;
    
};

#endif /* GridSBSolver_hpp */
