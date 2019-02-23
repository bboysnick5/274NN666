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


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
class UniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType> {
public:
    UniLatLngBKDTGridSBSolver(double = 1, size_t = 1500);
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lat, double lng) const override;
    virtual void printSolverInfo() const override final;
    
protected:
    const double AVE_LOC_PER_CELL;
    const size_t MAX_CACHE_CELL_VEC_SIZE;
    double lngInc, latInc, latIncInverse, sideLen;
    size_t totalLocSize, totalNodeSize = 0, singleLocs = 0,
           vecLocs = 0, rowSize, colSize;
    std::vector<std::variant<std::vector<std::pair<Point<3>, const SBLoc*>>,
                             const SBLoc*, KDT<KDTType>>> gridCache;
    void calcSideLenFromAlpc();
    void fillCacheCell(double, double, double,
                       std::vector<std::pair<Point<3>, const SBLoc*>>&);
    const SBLoc* returnNNLocFromCacheVariant(double, double,
          const std::variant<std::vector<std::pair<Point<3>, const SBLoc*>>,
          const SBLoc*, KDT<KDTType>>&) const;
    
private:
    virtual void fillGridCache();
};

/*

template <template <size_t, class, typename Point<3>::DistType> class KDTType>
class UniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType> {
public:
    UniLatLngBKDTGridSBSolver(double = 1, size_t = 1500);
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lng, double lat) const override;
    virtual void printSolverInfo() const override final;
    
protected:
    const double AVE_LOC_PER_CELL;
    const size_t MAX_CACHE_CELL_VEC_SIZE;
    double lngInc, latInc, sideLen;
    size_t totalLocSize, totalNodeSize = 0, singleLocs = 0,
    vecLocs = 0, rowSize, colSize;
    std::vector<std::variant<std::tuple<Point<3>*, const SBLoc**, size_t>,
    const SBLoc*, const KDT<KDTType>>> gridCache;
    void calcSideLenFromAlpc();
    void fillCacheCell(double, double, double,
                       std::vector<std::pair<Point<3>, const SBLoc*>>&);
    static const SBLoc* returnNNLocFromCacheVariant(double, double,
                                                    const std::variant<std::tuple<Point<3>*, const SBLoc**, size_t>,
                                                    const SBLoc*, const KDT<KDTType>>&);
    ~UniLatLngBKDTGridSBSolver();
private:
    virtual void fillGridCache();
}; */

#endif /* UniLatLngBKDTGridSBSolver_hpp */
