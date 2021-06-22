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


template <template <typename FPType, std::uint_fast8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>
class UniCellBKDTGridSBSolver final : public UniLatLngBKDTGridSBSolver<KDTType, FPType> {
public:
    UniCellBKDTGridSBSolver(FPType = 1.0, std::size_t = 1500);
    //void Build(std::span<const SBLoc<FPType>>) override final;
    const SBLoc<FPType>* FindNearestLoc(PointND<FPType, 2> geo_search_pt) const override;
    //virtual void PrintSolverInfo() const override;
    virtual ~UniCellBKDTGridSBSolver() override {}


private:
    std::vector<std::pair<std::size_t, FPType>> thisRowStartIdx;
    void FillGridCache() override final;
};



#endif /* UniCellBKDTGridSBSolver_hpp */
