//
//  UniCellBKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#ifndef UniCellBKDTGridSBSolver_hpp
#define UniCellBKDTGridSBSolver_hpp

#include <stdio.h>
#include <vector>
#include <memory>
#include <iterator>
#include "BKDTSBSolver.hpp"
#include "KDTree.hpp"


class UniCellBKDTGridSBSolver : public BKDTSBSolver {
public:
    UniCellBKDTGridSBSolver(double aveLocPerCell = 1);
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lng, double lat) const override;
    
    
private:
    const double AVE_LOC_PER_CELL = 1;
    
    std::vector<std::pair<std::unique_ptr<KDTree<3,const SBLoc*,DistType::EUC>>,
                          const SBLoc*>> gridCache;
    std::vector<std::pair<size_t, double>> thisRowStartIdx;
    double latInc;
    
    void fillGridCache(size_t, double);
    double calcSideLenFromAlpc();
};



#endif /* UniCellBKDTGridSBSolver_hpp */