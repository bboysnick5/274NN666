//
//  UnionUnionUniLatLngBKDTGridSBSolver.hpp
//  274F16NearestSB
//
//  Created by nick on 2/19/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef UnionUniLatLngBKDTGridSBSolver_hpp
#define UnionUniLatLngBKDTGridSBSolver_hpp

#include "BKDTSBSolver.hpp"
#include "KDTree.hpp"
#include <stdio.h>
#include <vector>
#include <iterator>


template <template <size_t, class, typename Point<3>::DistType> class KDTType>
class UnionUniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType> {
public:
    UnionUniLatLngBKDTGridSBSolver(double = 1, size_t = 1500);
    void build(const std::shared_ptr<std::vector<SBLoc>>&) override;
    const SBLoc* findNearest(double lat, double lng) const override;
    virtual void printSolverInfo() const override final;
    
    
    
    
protected:
    
    struct UnionCell {
        union {
            const SBLoc *cacheLoc;
            std::pair<Point<3>, const SBLoc*> *cacheLocs;
            const KDT<KDTType> *cacheTree;
        };
        
        UnionCell(size_t, const std::vector<std::pair<Point<3>, const SBLoc*>>&);
        ~UnionCell();
        size_t size() const;
        
    private:
        size_t _size;
    };
    
    struct BitCell {
        
        BitCell(size_t, const std::vector<std::pair<Point<3>, const SBLoc*>>&);
        ~BitCell();
        size_t size() const;
        size_t rawSize() const;
        
        const SBLoc* getSingleLoc() const;
        const std::pair<Point<3>, const SBLoc*>* getLocPairs() const;
        const KDT<KDTType>* getCacheTree() const;
        
    private:
        uintptr_t ptr;
    };
    
    /*
    struct CellAlloc {
        typedef Cell value_type;
        value_type* allocate(size_t n) { return static_cast<value_type*>(::operator new(sizeof(value_type) * n)); }
        void deallocate(value_type* p, size_t n) { return ::operator delete(static_cast<void*>(p)); }
        template<class U, class... Args>
        void construct(U* p, Args&&... args) { ::new(static_cast<void*>(p)) U{ std::forward<Args>(args)... }; }
    }; */
    
    const double AVE_LOC_PER_CELL;
    const size_t MAX_CACHE_CELL_VEC_SIZE;
    double lngInc, lngIncInverse, latInc, latIncInverse, sideLen;
    size_t totalLocSize, totalNodeSize = 0, singleLocs = 0,
    vecLocs = 0, rowSize, colSize;
    std::vector<BitCell> gridCache;
    
    void calcSideLenFromAlpc();
    void fillCacheCell(double, double, double,
                       std::vector<std::pair<Point<3>, const SBLoc*>>&);
    const SBLoc* returnNNLocFromCacheVariant(double, double, const BitCell&) const;
    
private:
    virtual void fillGridCache();
};

/*
 
 template <template <size_t, class, typename Point<3>::DistType> class KDTType>
 class UnionUniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType> {
 public:
 UnionUniLatLngBKDTGridSBSolver(double = 1, size_t = 1500);
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
 ~UnionUniLatLngBKDTGridSBSolver();
 private:
 virtual void fillGridCache();
 }; */

#endif /* UnionUniLatLngBKDTGridSBSolver_hpp */
