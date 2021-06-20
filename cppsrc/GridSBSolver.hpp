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


template <typename FPType>
class GridSBSolver : public SBSolver<FPType> {
    
public:
    void Build(std::span<const SBLoc<FPType>>) override;
    const SBLoc<FPType>* FindNearestLoc(const PointND<FPType, 2>&) const override;
    GridSBSolver(FPType aveLocPerCell = 1);
    virtual ~GridSBSolver() override {}
    
protected:
    void findKeyLngLat(std::span<const SBLoc<FPType>>);
    std::pair<std::size_t, std::size_t> getIdx(FPType lng, FPType lat) const;

    void constructGrid(std::span<const SBLoc<FPType>>);
    void fillGrid(std::span<const SBLoc<FPType>>);
    void NNOneCell(const std::unordered_set<const SBLoc<FPType>*>&,
                  FPType, FPType, FPType&, const SBLoc<FPType>*&) const;
    
    const FPType AVE_LOC_PER_CELL;
    std::vector<std::vector<std::unordered_set<const SBLoc<FPType>*>>> grid_;
    FPType side_len_, minLng, maxLng, minLat, maxLat, midLng, midLat;
    std::size_t row_size_, col_size_, numLocs = 0;
    
private:
    static constexpr FPType DISTORT_FACTOR = 0.95;
    
};

#endif /* GridSBSolver_hpp */

