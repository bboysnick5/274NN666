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
        
    void build(const std::vector<SBLoc> &sbData);
    SBLoc findNearest(double lng, double lat);
    
    
private:
    
    void constructGrid(const std::vector<SBLoc> &sbData,
                       const std::vector<double>& boundaryPts);
    std::pair<int, int> getIdx(double lng, double lat) const;
    void NNInCell(const std::unordered_set<SBLoc> &,
                  double, double, double&, SBLoc&);
    
    std::vector<std::vector<std::unordered_set<SBLoc>>> grid;
    double sideLen, midLng, midLat;
    int rowSize, colSize;
    
};

#endif /* GridSBSolver_hpp */
