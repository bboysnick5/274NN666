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


template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
class UnionUniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType, FPType> {
    
public:
    UnionUniLatLngBKDTGridSBSolver(FPType = 1, std::size_t = 1500);
    virtual void Build(std::span<const SBLoc<FPType>>) override;
    virtual const SBLoc<FPType>* FindNearestLoc(PointND<FPType, 2> geo_search_pt) const override;
    virtual void PrintSolverInfo() const override final;
    virtual ~UnionUniLatLngBKDTGridSBSolver() override {}
    
protected:
    
    struct BitCell;
    
    inline static FPType EUC3DDistSqFromLatCosDeltaLng(FPType lat1, FPType lat2, FPType cosDeltaLng);
    
    void calcSideLenFromAlpc();
    
    void LoopBodyThreadingPolicyDispatch();
    
    inline void FillCacheCell(const PointND<FPType, 2>&, FPType, std::size_t,
                              std::vector<typename KDT<KDTType, FPType>::node_type>&);
    
    inline void FillCacheCell(std::size_t idx_to_fill, const PointND<FPType, 2>&, FPType, std::size_t,
                              std::vector<typename KDT<KDTType, FPType>::node_type>&);

    const SBLoc<FPType>* ReturnNNLocFromCacheVariant(const PointND<FPType, 2>&, const BitCell&) const;
    
    const FPType kAveActualLocsPerCell_;
    const std::size_t kMaxCacheCellVecSize_;
    FPType lng_inc_, lng_inc_inverse_, lat_inc_, lat_inc_inverse_, side_len_;
    std::size_t row_size_, col_size_, num_actual_locs_;
    
    
    std::vector<BitCell> grid_cache_;
    
    /*
    struct CellAlloc {
        typedef Cell value_type;
        value_type* allocate(std::uint8_t N) { return static_cast<value_type*>(::operator new(sizeof(value_type) * n)); }
        void deallocate(value_type* p, std::uint8_t N) { return ::operator delete(static_cast<void*>(p)); }
        template<class U, class... Args>
        void construct(U* p, Args&&... args) { ::new(static_cast<void*>(p)) U{ std::forward<Args>(args)... }; }
    }; */
    



    
private:
    
    virtual void FillGridCache();
    
    virtual void LoopBody(def::Policy_Tag<def::ThreadingPolicy::kSingle>);
    
    virtual void LoopBody(def::Policy_Tag<def::ThreadingPolicy::kMultiOmp>);
    
    virtual void LoopBody(def::Policy_Tag<def::ThreadingPolicy::kMultiHand>);
    
    
    
    
protected:
    
    struct BitCell {
        
        using DataPtr = std::conditional_t<policy == def::ThreadingPolicy::kSingle, uintptr_t, std::atomic<uintptr_t>>;

        inline BitCell(uintptr_t otherPtr);
        
        BitCell(const BitCell&);
        BitCell(BitCell&&);

        inline BitCell& operator=(BitCell&&) noexcept;
        BitCell& operator=(const BitCell&) noexcept;

        inline BitCell(std::vector<typename KDT<KDTType, FPType>::node_type>&,
                       std::size_t, std::initializer_list<const BitCell*> prevCells);
        ~BitCell();
        
        template <typename T = DataPtr,
        typename std::enable_if_t<std::is_same<T, uintptr_t>::value, bool> = true>
        inline void SetPtr(uintptr_t);
        
        template <typename T = DataPtr,
        typename std::enable_if_t<std::is_same<T, std::atomic<uintptr_t>>::value, bool> = true>
        inline void SetPtr(uintptr_t);
        
        template <typename T = DataPtr,
        typename std::enable_if_t<std::is_same<T, uintptr_t>::value, bool> = true>
        inline uintptr_t GetPtr() const;
        
        template <typename T = DataPtr,
        typename std::enable_if_t<std::is_same<T, std::atomic<uintptr_t>>::value, bool> = true>
        inline uintptr_t GetPtr() const;
        
        inline static std::size_t size(uintptr_t);
        inline static std::size_t RawSizeBits(uintptr_t);
        inline static bool IsUniqueVecLoc(uintptr_t);

        
        inline static const SBLoc<FPType>* GetSingleLoc(uintptr_t);
        inline static const typename KDT<KDTType, FPType>::node_type* GetLocPairs(uintptr_t);
        inline static const KDT<KDTType, FPType>* GetCacheTree(uintptr_t);
        
    private:
        void DeAlloc();
        
        inline constexpr static uintptr_t MASK_OUT_16TH_BIT = ~(1ULL << 48);
        inline constexpr static uintptr_t MASK_OUT_LEAST_SIG_BIT = ~1ull;
        
        // magic
        std::conditional_t<policy == def::ThreadingPolicy::kSingle, uintptr_t, std::atomic<uintptr_t>> ptr;
    };
    
};


template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::BitCell(uintptr_t ptrVal) {
    SetPtr(ptrVal);
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::BitCell(BitCell&& rhs) {
    SetPtr(rhs.GetPtr());
    rhs.SetPtr(0);
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::BitCell(const BitCell& rhs) {
    SetPtr(rhs.GetPtr());
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline typename UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell&
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::operator=(BitCell&& rhs) noexcept {
    if (this != &rhs) [[likely]]{
        DeAlloc();
        SetPtr(rhs.GetPtr());
        rhs.SetPtr(0);
    }
    return *this;
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
typename UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell&
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::operator=(const BitCell& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        DeAlloc();
        SetPtr(rhs.GetPtr());
    }
    return *this;
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::
BitCell(std::vector<typename KDT<KDTType, FPType>::node_type> &pt_loc_vec,
        std::size_t max_cache_vec_size, std::initializer_list<const BitCell*> prev_cells) {
    if (std::size_t size = pt_loc_vec.size(); size == 1) {
        SetPtr((reinterpret_cast<std::uintptr_t>(pt_loc_vec[0].value) & MASK_OUT_16TH_BIT) | (1ull << 48));
    } else if (size < max_cache_vec_size) {
        for (const BitCell* prev_cell : prev_cells) {
            if (prev_cell) {
                if (uintptr_t cell_ptr_val = prev_cell->GetPtr();
                    cell_ptr_val && size == (cell_ptr_val >> 48)
                    && std::is_permutation(pt_loc_vec.cbegin(), pt_loc_vec.cend(), reinterpret_cast<typename KDT<KDTType, FPType>::node_type*>((static_cast<intptr_t>(cell_ptr_val << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT),
                              [](const auto &nh1, const auto &nh2){return nh1.value == nh2.value;})) {
                    SetPtr(cell_ptr_val & MASK_OUT_LEAST_SIG_BIT);
                    return;
                }
            }
        }
        //auto *cacheLocs = static_cast<typename KDT<KDTType, FPType>::node_type*>(
          //                ::operator new(size*sizeof(typename KDT<KDTType, FPType>::node_type), std::nothrow));
        auto *cacheLocs = static_cast<typename KDT<KDTType, FPType>::node_type*>(::operator new(size*sizeof(typename KDT<KDTType, FPType>::node_type)));
        std::uninitialized_move(std::make_move_iterator(pt_loc_vec.begin()), std::make_move_iterator(pt_loc_vec.end()), cacheLocs);
        SetPtr((reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK_OUT_16TH_BIT) | (size << 48) | 1ull);
    } else [[unlikely]] {
        SetPtr(reinterpret_cast<std::uintptr_t>(new KDT<KDTType, FPType>(pt_loc_vec.begin(), pt_loc_vec.end())) & MASK_OUT_16TH_BIT);
    }
}


template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::DeAlloc() {
    if (uintptr_t cell_ptr_val = GetPtr();
        cell_ptr_val != 0) [[likely]] {
        std::size_t size = cell_ptr_val >> 48;
        if (size > 1 && (cell_ptr_val & 1ull)) [[likely]] {
            typename KDT<KDTType, FPType>::node_type *rawCacheLocPtr = reinterpret_cast<typename KDT<KDTType, FPType>::node_type*>((static_cast<intptr_t>(cell_ptr_val << 16) >> 16) & MASK_OUT_LEAST_SIG_BIT);
            std::destroy_n(rawCacheLocPtr, size);
            ::operator delete(rawCacheLocPtr);
        } else if (size == 0) [[unlikely]] {
            delete reinterpret_cast<KDT<KDTType, FPType>*>(static_cast<intptr_t>(cell_ptr_val << 16) >> 16);
        }
    }
}


template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::~BitCell() {
    DeAlloc();
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
template <typename T, typename std::enable_if_t<std::is_same<T, uintptr_t>::value, bool>>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::GetPtr() const {
    return ptr;
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
template <typename T, typename std::enable_if_t<std::is_same<T, std::atomic<uintptr_t>>::value, bool>>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::GetPtr() const {
    return ptr.load(std::memory_order_relaxed);
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
template <typename T, typename std::enable_if_t<std::is_same<T, uintptr_t>::value, bool>>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::SetPtr(uintptr_t ptr_to_set) {
    ptr = ptr_to_set;
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
template <typename T, typename std::enable_if_t<std::is_same<T, std::atomic<uintptr_t>>::value, bool>>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::SetPtr(uintptr_t ptr_to_set) {
    ptr.store(ptr_to_set, std::memory_order_relaxed);
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::size(uintptr_t ptr) {
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

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::RawSizeBits(uintptr_t ptr) {
    return ptr >> 48;
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline bool UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::IsUniqueVecLoc(uintptr_t ptr) {
    return ptr & 1ull;
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline const SBLoc<FPType>* UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::GetSingleLoc(uintptr_t ptr) {
    return reinterpret_cast<const SBLoc<FPType>*>(static_cast<intptr_t>(ptr << 16) >> 16);
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline const typename KDT<KDTType, FPType>::node_type*
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::GetLocPairs(uintptr_t ptr) {
    return reinterpret_cast<typename KDT<KDTType, FPType>::node_type*>((static_cast<intptr_t>(ptr << 16) >> 16) & (std::numeric_limits<uintptr_t>::max() - 1));
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType, def::ThreadingPolicy policy>
inline const KDT<KDTType, FPType>*
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::BitCell::GetCacheTree(uintptr_t ptr) {
    return reinterpret_cast<const KDT<KDTType, FPType>*>(static_cast<intptr_t>(ptr << 16) >> 16);
}




template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
FillCacheCell(std::size_t idx_to_fill,
              const PointND<FPType, 2>& thisCtrGeoPt, FPType diagonalDistSq3DEUC, 
              std::size_t this_col_size,
              std::vector<typename KDT<KDTType, FPType>::node_type>& pt_loc_vec) {
    this->loc_kdt_.NNsWithFence(SBLoc<FPType>::GeoPtTo3dEucPt(thisCtrGeoPt),
        diagonalDistSq3DEUC, std::back_inserter(pt_loc_vec));
    grid_cache_[idx_to_fill] = {pt_loc_vec, kMaxCacheCellVecSize_,
                            {idx_to_fill > 0 ? &grid_cache_[idx_to_fill - 1] : nullptr,
                             idx_to_fill >= col_size_ ? &grid_cache_[idx_to_fill - col_size_] : nullptr,
                             idx_to_fill > col_size_ ? &grid_cache_[idx_to_fill - col_size_ - 1] : nullptr}};
    pt_loc_vec.clear();
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
FillCacheCell(const PointND<FPType, 2>& thisCtrGeoPt, FPType diagonalDistSq3DEUC, std::size_t this_col_size,
              std::vector<typename KDT<KDTType, FPType>::node_type>& pt_loc_vec) {
    this->loc_kdt_.NNsWithFence(SBLoc<FPType>::GeoPtTo3dEucPt(thisCtrGeoPt), diagonalDistSq3DEUC, std::back_inserter(pt_loc_vec));
    auto size = this->grid_cache_.size();
    this->grid_cache_.emplace_back(pt_loc_vec, kMaxCacheCellVecSize_,
                                 std::initializer_list<const BitCell*>{
                                    size > 0 ? &this->grid_cache_.back() : nullptr,
                                    size >= this_col_size ? &this->grid_cache_[size - this_col_size] : nullptr,
                                    size > this_col_size ? &this->grid_cache_[size - this_col_size - 1] : nullptr});
    pt_loc_vec.clear();
}

template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy policy>
inline FPType UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, policy>::
EUC3DDistSqFromLatCosDeltaLng(FPType lat1, FPType lat2, FPType cosDeltaLng) {
    return 2.0 - 2.0 * (sin(lat1)*sin(lat2) + std::cos(lat1)*std::cos(lat2)*cosDeltaLng);
}
/*
 
 template <template <typename FPType, std::uint8_t N, class, typename PointND<FPType, N>::DistType> class KDTType, typename FPType>

 class UnionUniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType> {
 public:
 UnionUniLatLngBKDTGridSBSolver(FPType = 1, std::size_t = 1500);
 void Build(std::span<const SBLoc<FPType>>) override;
 const SBLoc<FPType>* FindNearestLoc(FPType lng, FPType lat) const override;
 virtual void PrintSolverInfo() const override final;
 
 protected:
 const FPType AVE_LOC_PER_CELL;
 const std::size_t kMaxCacheCellVecSize_;
 FPType lng_inc_, lat_inc_, side_len_;
 std::size_t totalLocSize, totalNodeSize = 0, singleLocs = 0,
 vecLocs = 0, row_size_, col_size_;
 std::vector<std::variant<std::tuple<PointND<FPType, 3>*, const SBLoc<FPType>**, std::size_t>,
 const SBLoc<FPType>*, const KDT<KDTType>>> grid_cache_;
 void calcSideLenFromAlpc();
 void FillCacheCell(FPType, FPType, FPType,
 std::vector<typename KDT<KDTType, FPType>::node_type>&);
 static const SBLoc<FPType>* ReturnNNLocFromCacheVariant(FPType, FPType,
 const std::variant<std::tuple<PointND<FPType, 3>*, const SBLoc<FPType>**, std::size_t>,
 const SBLoc<FPType>*, const KDT<KDTType>>&);
 ~UnionUniLatLngBKDTGridSBSolver();
 private:
 virtual void FillGridCache();
 }; */

#endif /* UnionUniLatLngBKDTGridSBSolver_hpp */
