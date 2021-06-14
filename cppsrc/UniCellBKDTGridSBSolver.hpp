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


template <template <class DT, std::size_t N, class, typename PointND<DT, N>::DistType> class KDTType, class dist_type>
class UniCellBKDTGridSBSolver : public UniLatLngBKDTGridSBSolver<KDTType, dist_type> {
public:
    UniCellBKDTGridSBSolver(dist_type = 1.0, std::size_t = 1500);
    //void Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override final;
    const SBLoc<dist_type>* FindNearestLoc(const PointND<dist_type, 2>&) const override final;
    //virtual void PrintSolverInfo() const override;
    virtual ~UniCellBKDTGridSBSolver() override {}


private:
    std::vector<std::pair<std::size_t, dist_type>> thisRowStartIdx;
    void FillGridCache() override final;
};



#endif /* UniCellBKDTGridSBSolver_hpp */
