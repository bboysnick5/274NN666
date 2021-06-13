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
#include "Definition.hpp"

#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <iterator>
#include <type_traits>
#include <new>


template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
class UnionUniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType, dist_type> {
    
public:
    UnionUniLatLngBKDTGridSBSolver(dist_type = 1, std::size_t = 1500);
    void Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
    const SBLoc<dist_type>* FindNearestLoc(const PointND<dist_type, 2>&) const override;
    virtual void PrintSolverInfo() const override final;
    virtual ~UnionUniLatLngBKDTGridSBSolver() override {}
    
protected:
    
    struct BitCell;
    
    inline static dist_type EUC3DDistSqFromLatCosDeltaLng(dist_type lat1, dist_type lat2, dist_type cosDeltaLng);
    
    void calcSideLenFromAlpc();
    
    void LoopBodyThreadingPolicyDispatch(std::vector<typename KDT<KDTType, dist_type>::node_type>&);
    
    inline void FillCacheCell(const PointND<dist_type, 2>&, dist_type, std::size_t,
                       std::vector<typename KDT<KDTType, dist_type>::node_type>&);
    
    inline void FillCacheCell(std::size_t idxToFill, const PointND<dist_type, 2>&, dist_type, std::size_t,
                       std::vector<typename KDT<KDTType, dist_type>::node_type>&);

    const SBLoc<dist_type>* ReturnNNLocFromCacheVariant(const PointND<dist_type, 2>&, const BitCell&) const;
    
    const dist_type AVE_LOC_PER_CELL;
    const std::size_t MAX_CACHE_CELL_VEC_SIZE;
    dist_type lngInc, lngIncInverse, latInc, latIncInverse, sideLen;
    std::size_t rowSize, colSize, numGivenLocs;
    
    
    std::vector<BitCell> grid_cache_;
    
    /*
    struct CellAlloc {
        typedef Cell value_type;
        value_type* allocate(std::size_t n) { return static_cast<value_type*>(::operator new(sizeof(value_type) * n)); }
        void deallocate(value_type* p, std::size_t n) { return ::operator delete(static_cast<void*>(p)); }
        template<class U, class... Args>
        void construct(U* p, Args&&... args) { ::new(static_cast<void*>(p)) U{ std::forward<Args>(args)... }; }
    }; */
    



    
private:
    
    virtual void FillGridCache();
    
    virtual void LoopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs, def::Policy_Tag<def::ThreadingPolicy::kSingle>);
    
    virtual void LoopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs, def::Policy_Tag<def::ThreadingPolicy::kMultiOmp>);
    
    virtual void LoopBody(std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs, def::Policy_Tag<def::ThreadingPolicy::kMultiHand>);
    
    
    
    
protected:
    
    struct BitCell {
        
        using DataPtr = std::conditional_t<policy == def::ThreadingPolicy::kSingle, uintptr_t, std::atomic<uintptr_t>>;

        inline BitCell(uintptr_t otherPtr);
        
        BitCell(const BitCell&);
        BitCell(BitCell&&);

        inline BitCell& operator=(BitCell&&) noexcept;
        BitCell& operator=(const BitCell&) noexcept;

        inline BitCell(std::vector<typename KDT<KDTType, dist_type>::node_type>&,
                       std::size_t, std::initializer_list<const BitCell*> prevCells);
        ~BitCell();
        
        template <typename T = DataPtr,
        typename std::enable_if_t<std::is_same<T, uintptr_t>::value, bool> = true>
        inline void setPtr(uintptr_t);
        
        template <typename T = DataPtr,
        typename std::enable_if_t<std::is_same<T, std::atomic<uintptr_t>>::value, bool> = true>
        inline void setPtr(uintptr_t);
        
        template <typename T = DataPtr,
        typename std::enable_if_t<std::is_same<T, uintptr_t>::value, bool> = true>
        inline uintptr_t getPtr() const;
        
        template <typename T = DataPtr,
        typename std::enable_if_t<std::is_same<T, std::atomic<uintptr_t>>::value, bool> = true>
        inline uintptr_t getPtr() const;
        
        inline static std::size_t size(uintptr_t);
        inline static std::size_t RawSizeBits(uintptr_t);
        inline static bool IsUniqueVecLoc(uintptr_t);

        
        inline static const SBLoc<dist_type>* GetSingleLoc(uintptr_t);
        inline static const typename KDT<KDTType, dist_type>::node_type* GetLocPairs(uintptr_t);
        inline static const KDT<KDTType, dist_type>* GetCacheTree(uintptr_t);
        
    private:
        void deAlloc();
        
        inline constexpr static uintptr_t MASK_OUT_16TH_BIT = ~(1ULL << 48);
        inline constexpr static uintptr_t MASK_OUT_LEAST_SIG_BIT = ~1ull;
        
        // magic
        std::conditional_t<policy == def::ThreadingPolicy::kSingle, uintptr_t, std::atomic<uintptr_t>> ptr;
    };
    
};


template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::BitCell(uintptr_t ptrVal) {
    setPtr(ptrVal);
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::BitCell(BitCell&& rhs) {
    setPtr(rhs.getPtr());
    rhs.setPtr(0);
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::BitCell(const BitCell& rhs) {
    setPtr(rhs.getPtr());
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline typename UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell&
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::operator=(BitCell&& rhs) noexcept {
    if (this != &rhs) [[likely]]{
        deAlloc();
        setPtr(rhs.getPtr());
        rhs.setPtr(0);
    }
    return *this;
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
typename UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell&
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::operator=(const BitCell& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        deAlloc();
        setPtr(rhs.getPtr());
    }
    return *this;
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::
BitCell(std::vector<typename KDT<KDTType, dist_type>::node_type> &bufVec,
        std::size_t maxCacheVecSize, std::initializer_list<const BitCell*> prevCells) {
    if (std::size_t size = bufVec.size(); size == 1) {
        setPtr((reinterpret_cast<std::uintptr_t>(bufVec[0].value) & MASK_OUT_16TH_BIT) | (1ull << 48));
    } else if (size < maxCacheVecSize) {
        for (const BitCell* prevCell : prevCells) {
            if (prevCell) {
                if (uintptr_t cellPtrVal = prevCell->getPtr();
                    cellPtrVal && size == (cellPtrVal >> 48)
                    && std::is_permutation(bufVec.cbegin(), bufVec.cend(), reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(cellPtrVal << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT),
                              [](const auto &nh1, const auto &nh2){return nh1.value == nh2.value;})) {
                    setPtr(cellPtrVal & MASK_OUT_LEAST_SIG_BIT);
                    return;
                }
            }
        }
        //auto *cacheLocs = static_cast<typename KDT<KDTType, dist_type>::node_type*>(
          //                ::operator new(size*sizeof(typename KDT<KDTType, dist_type>::node_type), std::nothrow));
        auto *cacheLocs = static_cast<typename KDT<KDTType, dist_type>::node_type*>(::operator new(size*sizeof(typename KDT<KDTType, dist_type>::node_type), std::align_val_t{32}));
        std::uninitialized_move(std::make_move_iterator(bufVec.begin()), std::make_move_iterator(bufVec.end()), cacheLocs);
        setPtr((reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK_OUT_16TH_BIT) | (size << 48) | 1ull);
    } else [[unlikely]] {
        setPtr(reinterpret_cast<std::uintptr_t>(new KDT<KDTType, dist_type>(bufVec.begin(), bufVec.end())) & MASK_OUT_16TH_BIT);
    }
}


template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::deAlloc() {
    if (uintptr_t ptrValue = getPtr();
        ptrValue != 0) [[likely]] {
        std::size_t size = ptrValue >> 48;
        if (size > 1 && (ptrValue & 1ull)) [[likely]] {
            typename KDT<KDTType, dist_type>::node_type *rawCacheLocPtr = reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(ptrValue << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT);
            std::destroy_n(rawCacheLocPtr, size);
            ::operator delete(rawCacheLocPtr, std::align_val_t{32});
        } else if (size == 0) [[unlikely]] {
            delete reinterpret_cast<KDT<KDTType, dist_type>*>(static_cast<intptr_t>(ptrValue << 16) >> 16);
        }
    }
}


template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::~BitCell() {
    deAlloc();
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
template <typename T, typename std::enable_if_t<std::is_same<T, uintptr_t>::value, bool>>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::getPtr() const {
    return ptr;
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
template <typename T, typename std::enable_if_t<std::is_same<T, std::atomic<uintptr_t>>::value, bool>>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::getPtr() const {
    return ptr.load(std::memory_order_relaxed);
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
template <typename T, typename std::enable_if_t<std::is_same<T, uintptr_t>::value, bool>>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::setPtr(uintptr_t ptrToSet) {
    ptr = ptrToSet;
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
template <typename T, typename std::enable_if_t<std::is_same<T, std::atomic<uintptr_t>>::value, bool>>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::setPtr(uintptr_t ptrToSet) {
    ptr.store(ptrToSet, std::memory_order_relaxed);
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::size(uintptr_t ptr) {
    if (ptr) [[likely]] {
        if (std::uintptr_t size = ptr >> 48;
            size != 0) [[likely]] {
            return size;
        } else {
            return GetCacheTree(ptr)->size();
        }
    } else {
        return 0;
    }
    // return ptr ? (ptr >> 48 ? ptr >> 48 : getCacheTree(ptr)->size()) : 0; attribute not working with ternary operator?
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::RawSizeBits(uintptr_t ptr) {
    return ptr >> 48;
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline bool UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::IsUniqueVecLoc(uintptr_t ptr) {
    return ptr & 1ull;
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline const SBLoc<dist_type>* UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::GetSingleLoc(uintptr_t ptr) {
    return reinterpret_cast<const SBLoc<dist_type>*>(static_cast<intptr_t>(ptr << 16) >> 16);
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline const typename KDT<KDTType, dist_type>::node_type*
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::GetLocPairs(uintptr_t ptr) {
    return reinterpret_cast<typename KDT<KDTType, dist_type>::node_type*>((static_cast<intptr_t>(ptr << 16) >> 16) & std::numeric_limits<uintptr_t>::max() - 1);
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type, def::ThreadingPolicy policy>
inline const KDT<KDTType, dist_type>*
UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::BitCell::GetCacheTree(uintptr_t ptr) {
    return reinterpret_cast<const KDT<KDTType, dist_type>*>(static_cast<intptr_t>(ptr << 16) >> 16);
}




template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
FillCacheCell(std::size_t idxToFill,
              const PointND<dist_type, 2>& thisCtrGeoPt, dist_type diagonalDistSq3DEUC, std::size_t thisColSize,
              std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt(thisCtrGeoPt),
        diagonalDistSq3DEUC, std::back_inserter(ptLocPairs));
    grid_cache_[idxToFill] = {ptLocPairs, MAX_CACHE_CELL_VEC_SIZE,
                            {idxToFill > 0 ? &grid_cache_[idxToFill - 1] : nullptr,
                             idxToFill >= colSize ? &grid_cache_[idxToFill - colSize] : nullptr,
                             idxToFill > colSize ? &grid_cache_[idxToFill - colSize - 1] : nullptr}};
    ptLocPairs.clear();
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
FillCacheCell(const PointND<dist_type, 2>& thisCtrGeoPt, dist_type diagonalDistSq3DEUC, std::size_t thisColSize,
              std::vector<typename KDT<KDTType, dist_type>::node_type>& ptLocPairs) {
    this->locKdt.rangeDiffKNNPairs(SBLoc<dist_type>::geoPtToCart3DPt(thisCtrGeoPt), diagonalDistSq3DEUC, std::back_inserter(ptLocPairs));
    auto size = this->grid_cache_.size();
    this->grid_cache_.emplace_back(ptLocPairs, MAX_CACHE_CELL_VEC_SIZE,
                                 std::initializer_list<const BitCell*>{
                                    size > 0 ? &this->grid_cache_.back() : nullptr,
                                    size >= thisColSize ? &this->grid_cache_[size - thisColSize] : nullptr,
                                    size > thisColSize ? &this->grid_cache_[size - thisColSize - 1] : nullptr});
    ptLocPairs.clear();
}

template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType,
          class dist_type, def::ThreadingPolicy policy>
inline dist_type UnionUniLatLngBKDTGridSBSolver<KDTType, dist_type, policy>::
EUC3DDistSqFromLatCosDeltaLng(dist_type lat1, dist_type lat2, dist_type cosDeltaLng) {
    return 2.0 - 2.0 * (sin(lat1)*sin(lat2) + std::cos(lat1)*std::cos(lat2)*cosDeltaLng);
}
/*
 
 template <template <class DT, std::size_t, class, typename PointND<DT, 3>::DistType> class KDTType, class dist_type>

 class UnionUniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType> {
 public:
 UnionUniLatLngBKDTGridSBSolver(dist_type = 1, std::size_t = 1500);
 void Build(const std::shared_ptr<std::vector<SBLoc<dist_type>>>&) override;
 const SBLoc<dist_type>* FindNearestLoc(dist_type lng, dist_type lat) const override;
 virtual void PrintSolverInfo() const override final;
 
 protected:
 const dist_type AVE_LOC_PER_CELL;
 const std::size_t MAX_CACHE_CELL_VEC_SIZE;
 dist_type lngInc, latInc, sideLen;
 std::size_t totalLocSize, totalNodeSize = 0, singleLocs = 0,
 vecLocs = 0, rowSize, colSize;
 std::vector<std::variant<std::tuple<PointND<dist_type, 3>*, const SBLoc<dist_type>**, std::size_t>,
 const SBLoc<dist_type>*, const KDT<KDTType>>> grid_cache_;
 void calcSideLenFromAlpc();
 void FillCacheCell(dist_type, dist_type, dist_type,
 std::vector<typename KDT<KDTType, dist_type>::node_type>&);
 static const SBLoc<dist_type>* ReturnNNLocFromCacheVariant(dist_type, dist_type,
 const std::variant<std::tuple<PointND<dist_type, 3>*, const SBLoc<dist_type>**, std::size_t>,
 const SBLoc<dist_type>*, const KDT<KDTType>>&);
 ~UnionUniLatLngBKDTGridSBSolver();
 private:
 virtual void FillGridCache();
 }; */

#endif /* UnionUniLatLngBKDTGridSBSolver_hpp */
