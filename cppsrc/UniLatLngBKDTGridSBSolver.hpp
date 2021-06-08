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
#include <variant>


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
class UniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType, dist_type> {
public:
    UniLatLngBKDTGridSBSolver(dist_type = 1, size_t = 1500);
    void build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
    const SBLoc<dist_type>* findNearest(const Point<dist_type, 2>&) const override;
    virtual void printSolverInfo() const override final;
    virtual ~UniLatLngBKDTGridSBSolver() override {}

    
protected:
    const dist_type AVE_LOC_PER_CELL;
    const size_t MAX_CACHE_CELL_VEC_SIZE;
    dist_type lngInc, latInc, latIncInverse, sideLen;
    size_t totalLocSize, totalNodeSize = 0, singleLocs = 0,
           vecLocs = 0, rowSize, colSize;
    std::vector<std::variant<std::vector<typename KDT<KDTType, dist_type>::node_type>,
                             const SBLoc<dist_type>*, KDT<KDTType, dist_type>>> gridCache;
    void calcSideLenFromAlpc();
    void fillCacheCell(dist_type, dist_type, dist_type,
                       std::vector<typename KDT<KDTType, dist_type>::node_type>&);
    const SBLoc<dist_type>* returnNNLocFromCacheVariant(const Point<dist_type, 2>&,
          const std::variant<std::vector<typename KDT<KDTType, dist_type>::node_type>,
          const SBLoc<dist_type>*, KDT<KDTType, dist_type>>&) const;
    
private:
    virtual void fillGridCache();
};

#endif /* UniLatLngBKDTGridSBSolver_hpp */
