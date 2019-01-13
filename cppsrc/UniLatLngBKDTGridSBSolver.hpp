//
//  UniLatLngBKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 12/29/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#ifndef UniLatLngBKDTGridSBSolver_hpp
#define UniLatLngBKDTGridSBSolver_hpp


#include "BKDTSBSolver.hpp"
#include "KDTree.hpp"
#include <stdio.h>
#include <vector>
#include <iterator>


template <template <size_t, typename, typename Point<3>::DistType> class Tree>
class UniLatLngBKDTGridSBSolver : public BKDTSBSolver<Tree> {
public:
    UniLatLngBKDTGridSBSolver(double aveLocPerCell = 1);
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lng, double lat) const override;
    virtual void printSolverInfo() const override final;
    
protected:
    const double AVE_LOC_PER_CELL = 1;
    double lngInc, latInc, sideLen;
    size_t totalNodeSize = 0, singleLocs = 0, rowSize, colSize;
    std::vector<std::pair<Tree<3,const SBLoc*, Point<3>::DistType::EUC>,
                const SBLoc*>> gridCache;
    
    double calcSideLenFromAlpc();
    
private:
    void fillGridCache();
};


#endif /* UniLatLngBKDTGridSBSolver_hpp */
