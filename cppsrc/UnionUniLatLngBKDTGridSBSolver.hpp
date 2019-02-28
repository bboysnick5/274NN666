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


template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>
class UnionUniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType, dist_type> {
public:
    UnionUniLatLngBKDTGridSBSolver(dist_type = 1, size_t = 1500);
    void build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
    const SBLoc<dist_type>* findNearest(const Point<dist_type, 2>&) const override;
    virtual void printSolverInfo() const override final;
    
    
    
    
protected:
    
    struct UnionCell {
        union {
            const SBLoc<dist_type> *cacheLoc;
            typename KDT<KDTType, dist_type>::node_type *cacheLocs;
            const KDT<KDTType, dist_type> *cacheTree;
        };
        
        UnionCell(const std::vector<typename KDT<KDTType, dist_type>::node_type>&, size_t);
        ~UnionCell();
        size_t size() const;
        
    private:
        size_t _size;
    };
    
    struct BitCell {
        
        BitCell(const std::vector<typename KDT<KDTType, dist_type>::node_type>&,
                size_t, const BitCell&, const BitCell&);
        BitCell(const std::vector<typename KDT<KDTType, dist_type>::node_type>&,
                size_t, const BitCell&);
        BitCell(const std::vector<typename KDT<KDTType, dist_type>::node_type>&,
                size_t);
        ~BitCell();
        size_t size() const;
        size_t rawSize() const;
        
        const SBLoc<dist_type>* getSingleLoc() const;
        const typename KDT<KDTType, dist_type>::node_type* getLocPairs() const;
        const KDT<KDTType, dist_type>* getCacheTree() const;
        
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
    
    const dist_type AVE_LOC_PER_CELL;
    const size_t MAX_CACHE_CELL_VEC_SIZE;
    dist_type lngInc, lngIncInverse, latInc, latIncInverse, sideLen;
    size_t totalLocSize, totalNodeSize = 0, singleLocs = 0,
    vecLocs = 0, rowSize, colSize;
    std::vector<BitCell> gridCache;
    
    void calcSideLenFromAlpc();
    void fillCacheCell(const Point<dist_type, 2>&, dist_type, size_t,
                       std::vector<typename KDT<KDTType, dist_type>::node_type>&);
    const SBLoc<dist_type>* returnNNLocFromCacheVariant(const Point<dist_type, 2>&, const BitCell&) const;
    
private:
    virtual void fillGridCache();
};

/*
 
 template <template <class DT, size_t, class, typename Point<DT, 3>::DistType> class KDTType, class dist_type>

 class UnionUniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType> {
 public:
 UnionUniLatLngBKDTGridSBSolver(dist_type = 1, size_t = 1500);
 void build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
 const SBLoc<dist_type>* findNearest(dist_type lng, dist_type lat) const override;
 virtual void printSolverInfo() const override final;
 
 protected:
 const dist_type AVE_LOC_PER_CELL;
 const size_t MAX_CACHE_CELL_VEC_SIZE;
 dist_type lngInc, latInc, sideLen;
 size_t totalLocSize, totalNodeSize = 0, singleLocs = 0,
 vecLocs = 0, rowSize, colSize;
 std::vector<std::variant<std::tuple<Point<dist_type, 3>*, const SBLoc<dist_type>**, size_t>,
 const SBLoc<dist_type>*, const KDT<KDTType>>> gridCache;
 void calcSideLenFromAlpc();
 void fillCacheCell(dist_type, dist_type, dist_type,
 std::vector<typename KDT<KDTType, dist_type>::node_type>&);
 static const SBLoc<dist_type>* returnNNLocFromCacheVariant(dist_type, dist_type,
 const std::variant<std::tuple<Point<dist_type, 3>*, const SBLoc<dist_type>**, size_t>,
 const SBLoc<dist_type>*, const KDT<KDTType>>&);
 ~UnionUniLatLngBKDTGridSBSolver();
 private:
 virtual void fillGridCache();
 }; */

#endif /* UnionUniLatLngBKDTGridSBSolver_hpp */
