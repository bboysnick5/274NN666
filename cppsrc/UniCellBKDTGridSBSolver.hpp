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
#include <iterator>
#include "GridSBSolver.hpp"
#include "KDTree.hpp"


class UniCellBKDTGridSBSolver : public GridSBSolver {
public:
    UniCellBKDTGridSBSolver(double aveLocPerCell = 1);
    void build();
    const SBLoc* findNearest(double lng, double lat) const;
    
    
private:
    std::vector<KDTree<3, const SBLoc*, DistType::EUC>> gridTreeCache;
    std::vector<const SBLoc*> gridSingleCache;
    std::vector<std::pair<size_t, double>> thisRowStartIdx;
    double latInc;
    
    void fillGridCache(const KDTree<3, const SBLoc*, DistType::EUC>&);
    double calcSideLenFromAlpc();
};



#endif /* UniCellBKDTGridSBSolver_hpp */
