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
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lat, double lng) const override;
    GridSBSolver(double aveLocPerCell = 1);
    
protected:
    void findKeyLngLat(const std::shared_ptr<std::vector<SBLoc>>&);
    std::pair<size_t, size_t> getIdx(double lng, double lat) const;

    void constructGrid(const std::shared_ptr<std::vector<SBLoc>>&);
    void fillGrid(const std::shared_ptr<std::vector<SBLoc>>&);
    void NNOneCell(const std::unordered_set<const SBLoc*>&,
                  double, double, double&, const SBLoc*&) const;
    
    const double AVE_LOC_PER_CELL;
    std::vector<std::vector<std::unordered_set<const SBLoc*>>> grid;
    double sideLen, minLng, maxLng, minLat, maxLat, midLng, midLat;
    size_t rowSize, colSize, numLocs = 0;
    
private:
    static constexpr double DISTORT_FACTOR = 0.95;
    
};

#endif /* GridSBSolver_hpp */

