//
//  UniCellBKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 12/30/18.
//  Copyright © 2018 Yunlong Liu. All rights reserved.
//

#ifndef UniCellBKDTGridSBSolver_hpp
#define UniCellBKDTGridSBSolver_hpp

#include "UniLatLngBKDTGridSBSolver.hpp"
#include "KDTree.hpp"
#include <stdio.h>
#include <vector>
#include <memory>
#include <iterator>


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
class UniCellBKDTGridSBSolver : public UniLatLngBKDTGridSBSolver<KDTType> {
public:
    UniCellBKDTGridSBSolver(double = 1, size_t = 40);
    //void build(const std::shared_ptr<std::vector<SBLoc>>&) override final;
    const SBLoc* findNearest(double lng, double lat) const override final;
    //virtual void printSolverInfo() const override;
    

private:
    std::vector<std::pair<size_t, double>> thisRowStartIdx;
    void fillGridCache() override final;
};



#endif /* UniCellBKDTGridSBSolver_hpp */
