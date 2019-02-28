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


template <typename dist_type>
class GridSBSolver : public SBSolver<dist_type> {
    
public:
    void build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
    const SBLoc<dist_type>* findNearest(const Point<dist_type, 2>&) const override;
    GridSBSolver(dist_type aveLocPerCell = 1);
    
protected:
    void findKeyLngLat(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&);
    std::pair<size_t, size_t> getIdx(dist_type lng, dist_type lat) const;

    void constructGrid(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&);
    void fillGrid(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&);
    void NNOneCell(const std::unordered_set<const SBLoc<dist_type>*>&,
                  dist_type, dist_type, dist_type&, const SBLoc<dist_type>*&) const;
    
    const dist_type AVE_LOC_PER_CELL;
    std::vector<std::vector<std::unordered_set<const SBLoc<dist_type>*>>> grid;
    dist_type sideLen, minLng, maxLng, minLat, maxLat, midLng, midLat;
    size_t rowSize, colSize, numLocs = 0;
    
private:
    static constexpr dist_type DISTORT_FACTOR = 0.95;
    
};

#endif /* GridSBSolver_hpp */

