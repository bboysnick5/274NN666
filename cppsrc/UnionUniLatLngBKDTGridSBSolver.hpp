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


template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy Policy>
class UnionUniLatLngBKDTGridSBSolver : public BKDTSBSolver<KDTType, FPType, Policy> {
    
public:
    UnionUniLatLngBKDTGridSBSolver(FPType = 1, std::size_t = 1500);
    virtual void Build(std::span<const SBLoc<FPType>>) override;
    virtual const SBLoc<FPType>* FindNearestLoc(typename SBLoc<FPType>::GeoPtType geo_search_pt) const override;
    virtual void PrintSolverInfo() const override final;
    virtual ~UnionUniLatLngBKDTGridSBSolver() override {}
    
protected:
    
    struct BitCell;
    
    inline static FPType EUC3DDistSqFromLatCosDeltaLng(FPType lat1, FPType lat2, FPType cosDeltaLng);
    
    void calcSideLenFromAlpc();
    
    void LoopBodyThreadingPolicyDispatch();
    
    inline void FillCacheCell(const typename SBLoc<FPType>::GeoPtType&, FPType, std::size_t,
                              std::vector<typename KDT<KDTType, FPType>::node_type>&);
    
    inline void FillCacheCell(std::size_t idx_to_fill, const typename SBLoc<FPType>::GeoPtType&, FPType, std::size_t,
                              std::vector<typename KDT<KDTType, FPType>::node_type>&);

    const SBLoc<FPType>* ReturnNNLocFromCacheVariant(const typename SBLoc<FPType>::GeoPtType&, const BitCell&) const;
    
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
    
    virtual void LoopBody(def::ThreadingPolicyTag<def::ThreadingPolicy::kSingle>);
    
    virtual void LoopBody(def::ThreadingPolicyTag<def::ThreadingPolicy::kMultiOmp>);
    
    virtual void LoopBody(def::ThreadingPolicyTag<def::ThreadingPolicy::kMultiHand>);
    
    
    
    
protected:
    
    struct BitCell {
        
        using DataPtr = std::conditional_t<Policy == def::ThreadingPolicy::kSingle, std::uintptr_t, std::atomic<std::uintptr_t>>;

        inline BitCell(std::uintptr_t);
        
        BitCell(const BitCell&);
        BitCell(BitCell&&);

        inline BitCell& operator=(BitCell&&) noexcept;
        BitCell& operator=(const BitCell&) noexcept;

        inline BitCell(std::vector<typename KDT<KDTType, FPType>::node_type>&,
                       std::size_t, std::initializer_list<const BitCell*> prevCells);
        ~BitCell();
        
        inline void SetRawPtr(std::uintptr_t);
        inline std::uintptr_t GetRawPtr() const;
        
        inline static std::size_t size(std::uintptr_t);
        inline static std::size_t RawSizeBits(std::uintptr_t);
        inline static bool IsUniqueVecLoc(std::uintptr_t);

        
        inline static const SBLoc<FPType>* GetSingleLoc(std::uintptr_t);
        inline static const typename KDT<KDTType, FPType>::node_type* GetLocPairs(std::uintptr_t);
        inline static const KDT<KDTType, FPType>* GetCacheTree(std::uintptr_t);
        
    private:
        void DeAlloc();
        
        inline constexpr static std::uintptr_t MASK_OUT_16TH_BIT = ~(1ULL << 48);
        inline constexpr static std::uintptr_t MASK_OUT_LEAST_SIG_BIT = ~1ull;
        
        // magic
        std::conditional_t<Policy == def::ThreadingPolicy::kSingle, std::uintptr_t, std::atomic<std::uintptr_t>> ptr_;
    };
    
};


template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::BitCell(std::uintptr_t ptrVal) {
    SetRawPtr(ptrVal);
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::BitCell(BitCell&& rhs) {
    SetRawPtr(rhs.GetRawPtr());
    rhs.SetRawPtr(0);
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::BitCell(const BitCell& rhs) {
    SetRawPtr(rhs.GetRawPtr());
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline typename UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell&
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::operator=(BitCell&& rhs) noexcept {
    if (this != &rhs) [[likely]]{
        DeAlloc();
        SetRawPtr(rhs.GetRawPtr());
        rhs.SetRawPtr(0);
    }
    return *this;
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
typename UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell&
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::operator=(const BitCell& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        DeAlloc();
        SetRawPtr(rhs.GetRawPtr());
    }
    return *this;
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::
BitCell(std::vector<typename KDT<KDTType, FPType>::node_type> &pt_loc_vec,
        std::size_t max_cache_vec_size, std::initializer_list<const BitCell*> prev_cells) {
    if (std::size_t size = pt_loc_vec.size(); size == 1) {
        SetRawPtr((reinterpret_cast<std::uintptr_t>(pt_loc_vec[0].value) & MASK_OUT_16TH_BIT) | (1ull << 48));
    } else if (size < max_cache_vec_size) {
        for (const BitCell* prev_cell : prev_cells) {
            if (prev_cell) {
                if (std::uintptr_t cell_ptr_val = prev_cell->GetRawPtr();
                    cell_ptr_val && size == (cell_ptr_val >> 48)
                    && std::is_permutation(pt_loc_vec.cbegin(), pt_loc_vec.cend(), reinterpret_cast<typename KDT<KDTType,
                                           FPType>::node_type*>((static_cast<intptr_t>(cell_ptr_val << 16) >> 16)
                                                                & MASK_OUT_LEAST_SIG_BIT),
                                           [](const auto &nh1, const auto &nh2){return nh1.value == nh2.value;})) {
                    SetRawPtr(cell_ptr_val & MASK_OUT_LEAST_SIG_BIT);
                    return;
                }
            }
        }
        //auto *cacheLocs = static_cast<typename KDT<KDTType, FPType>::node_type*>(
          //                ::operator new(size*sizeof(typename KDT<KDTType, FPType>::node_type), std::nothrow));
        auto *cacheLocs = static_cast<typename KDT<KDTType, FPType>::node_type*>(::operator new(size*sizeof(typename KDT<KDTType, FPType>::node_type)));
        std::uninitialized_move(std::make_move_iterator(pt_loc_vec.begin()), std::make_move_iterator(pt_loc_vec.end()), cacheLocs);
        SetRawPtr((reinterpret_cast<std::uintptr_t>(cacheLocs) & MASK_OUT_16TH_BIT) | (size << 48) | 1ull);
    } else [[unlikely]] {
        SetRawPtr(reinterpret_cast<std::uintptr_t>(new KDT<KDTType, FPType>(pt_loc_vec.begin(), pt_loc_vec.end())) & MASK_OUT_16TH_BIT);
    }
}


template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::DeAlloc() {
    if (std::uintptr_t cell_ptr_val = GetRawPtr();
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


template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::~BitCell() {
    DeAlloc();
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::GetRawPtr() const {
    if constexpr (std::is_same_v<decltype(ptr_), std::uintptr_t>) {
        return ptr_;
    } else {
        return ptr_.load(std::memory_order_relaxed);
    }
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::SetRawPtr(std::uintptr_t ptr_to_set) {
    if constexpr (std::is_same_v<decltype(ptr_), std::uintptr_t>) {
        ptr_ = ptr_to_set;
    } else {
        ptr_.store(ptr_to_set, std::memory_order_relaxed);
    }
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::size(std::uintptr_t ptr) {
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

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline std::size_t UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::RawSizeBits(std::uintptr_t ptr) {
    return ptr >> 48;
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline bool UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::IsUniqueVecLoc(std::uintptr_t ptr) {
    return ptr & 1ull;
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline const SBLoc<FPType>* UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::GetSingleLoc(std::uintptr_t ptr) {
    return reinterpret_cast<const SBLoc<FPType>*>(static_cast<intptr_t>(ptr << 16) >> 16);
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline const typename KDT<KDTType, FPType>::node_type*
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::GetLocPairs(std::uintptr_t ptr) {
    return reinterpret_cast<typename KDT<KDTType, FPType>::node_type*>((static_cast<intptr_t>(ptr << 16) >> 16) & (std::numeric_limits<std::uintptr_t>::max() - 1));
}

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType, def::ThreadingPolicy Policy>
inline const KDT<KDTType, FPType>*
UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::BitCell::GetCacheTree(std::uintptr_t ptr) {
    return reinterpret_cast<const KDT<KDTType, FPType>*>(static_cast<intptr_t>(ptr << 16) >> 16);
}




template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy Policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::
FillCacheCell(std::size_t idx_to_fill,
              const typename SBLoc<FPType>::GeoPtType& thisCtrGeoPt, FPType diagonalDistSq3DEUC, 
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

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy Policy>
inline void UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::
FillCacheCell(const typename SBLoc<FPType>::GeoPtType& thisCtrGeoPt, FPType diagonalDistSq3DEUC, std::size_t this_col_size,
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

template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType,
          typename FPType, def::ThreadingPolicy Policy>
inline FPType UnionUniLatLngBKDTGridSBSolver<KDTType, FPType, Policy>::
EUC3DDistSqFromLatCosDeltaLng(FPType lat1, FPType lat2, FPType cosDeltaLng) {
    return 2.0 - 2.0 * (sin(lat1)*sin(lat2) + std::cos(lat1)*std::cos(lat2)*cosDeltaLng);
}
/*
 
 template <template <typename FPType, std::uint8_t N, class, typename def::DistType> class KDTType, typename FPType>

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
 std::vector<std::variant<std::tuple<SBLoc<FPType>::GeoCart3DPtType*, const SBLoc<FPType>**, std::size_t>,
 const SBLoc<FPType>*, const KDT<KDTType>>> grid_cache_;
 void calcSideLenFromAlpc();
 void FillCacheCell(FPType, FPType, FPType,
 std::vector<typename KDT<KDTType, FPType>::node_type>&);
 static const SBLoc<FPType>* ReturnNNLocFromCacheVariant(FPType, FPType,
 const std::variant<std::tuple<SBLoc<FPType>::GeoCart3DPtType*, const SBLoc<FPType>**, std::size_t>,
 const SBLoc<FPType>*, const KDT<KDTType>>&);
 ~UnionUniLatLngBKDTGridSBSolver();
 private:
 virtual void FillGridCache();
 }; */

#endif /* UnionUniLatLngBKDTGridSBSolver_hpp */
