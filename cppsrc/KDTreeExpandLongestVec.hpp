//
//  KDTreeExpandLongestVec.hpp
//  274F16NearestSB
//
//  Created by nick on 2/12/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef KDTreeExpandLongestVec_hpp
#define KDTreeExpandLongestVec_hpp

#include "Point.hpp"
#include "Utility.hpp"
//#include "PoolAllocator.hpp"



#include "BoundedPQueue.hpp"

#include <immintrin.h>
//#include <Eigen/Core>
//#include "xsimd/xsimd.hpp"
//#include <llama/llama.hpp>
//#include <absl/container/fixed_array.h>
//#include <oneapi/dpl/execution>
//#include <oneapi/dpl/algorithm>
#include <stdexcept>
#include <unordered_map>
#include <map>
#include <utility>
#include <set>
#include <new>
#include <tuple>
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>
#include <functional>
#include <iostream>
#include <memory>
#include <array>
#include <concepts>
#include <cstdlib>
#include <bit>
#include <array>

//#include <ranges>

namespace ns {
    struct NodePtr {};
    struct PtNodePtrArray {};
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
class KDTreeExpandLongestVec {
public:
    
    typedef ElemType*                             value_iterator;
    typedef const ElemType*                       const_value_iterator;

    
    struct node_type {
        PointND<FPType, N> key;
        ElemType value;
    };
    
    
    // Constructor: KDTreeExpandLongestVec();
    // Usage: KDTreeExpandLongestVec<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTreeExpandLongestVec.
    KDTreeExpandLongestVec();
    
    // Constructor: KDTreeExpandLongestVec(FwdItType begin, FwdItType end);
    // Usage: KDTreeExpandLongestVec<3, int> myTree(vec.begin(), bec.end());
    // ----------------------------------------------------
    // Constructs a KDTreeExpandLongestVec from a collection. The tree will
    // be balanced using median constructing method
    // NOTE: The tree will not eliminate duplicates and the
    //       intended behavior will not be comprimised, tho
    //       less efficient with extra wasteful space.
    //template <std::random_access_iterator RAI> requires def::non_const_iterator<RAI>
    //requires std::same_as<typename std::iterator_traits<RAI>::value_type, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
    template <typename RAI, std::enable_if_t<!std::is_const_v<typename std::remove_pointer_t<typename std::iterator_traits<RAI>::pointer>>, int> = 0>
    KDTreeExpandLongestVec(RAI, RAI);
    
    //template <std::random_access_iterator ConstRAI> requires def::const_iterator<ConstRAI>
    //&& std::same_as<typename std::iter_value_t<ConstRAI>, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
    template <typename ConstRAI, std::enable_if_t<std::is_const_v<typename std::remove_pointer_t<typename std::iterator_traits<ConstRAI>::pointer>>, int> = 0>
    KDTreeExpandLongestVec(ConstRAI, ConstRAI);
    
    // Destructor: ~KDTreeExpandLongestVec()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTreeExpandLongestVec.
    ~KDTreeExpandLongestVec();
    
    // KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs);
    // KDTreeExpandLongestVec& operator=(const KDTreeExpandLongestVec& rhs);
    // Usage: KDTreeExpandLongestVec<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Copy constructor and copy assignment operator.
    KDTreeExpandLongestVec(const KDTreeExpandLongestVec&);
    KDTreeExpandLongestVec& operator=(const KDTreeExpandLongestVec&) noexcept;
    
    // KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs);
    // KDTreeExpandLongestVec& operator=(const KDTreeExpandLongestVec& rhs);
    // Usage: KDTreeExpandLongestVec<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Move constructor and move assignment operator.
    KDTreeExpandLongestVec(KDTreeExpandLongestVec&&) noexcept;
    KDTreeExpandLongestVec& operator=(KDTreeExpandLongestVec&&) noexcept;
    
    // std::uint8_t Dim() const;
    // Usage: std::uint8_t dim = kd.Dim();
    // ----------------------------------------------------
    // Returns the Dim of the tree.
    constexpr std::uint8_t dimension() const;
    typename def::DistType distType() const;
    
    // std::uint32_t size() const;
    // std::uint32_t height() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree, the max height
    // and whether the tree is empty.
    std::uint32_t size() const;
    std::uint32_t cap() const;
    std::uint32_t height() const;
    bool Empty() const;
    
    bool operator==(const KDTreeExpandLongestVec& rhs) const;
    bool operator!=(const KDTreeExpandLongestVec& rhs) const;

    
    void Clear();
    void ShrinkToFit();
    
    void PrintTreeInfo() const;
    
    // bool contains(const PointND<FPType, N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTreeExpandLongestVec.
    bool contains(const PointND<FPType, N>&) const;
    
    // void insert(const PointND<FPType, N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTreeExpandLongestVec, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const PointND<FPType, N>&, const ElemType&);
    
    // ElemType& operator[](const PointND<FPType, N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTreeExpandLongestVec.
    // If the point does not exist, then it is added to the KDTreeExpandLongestVec using the
    // default value of ElemType as its key.
    ElemType& operator[](const PointND<FPType, N>& pt);
    
    // ElemType& at(const PointND<FPType, N>& pt);
    // const ElemType& at(const PointND<FPType, N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function //throws an out_of_range exception.
    ElemType& at(const PointND<FPType, N>& pt);
    const ElemType& at(const PointND<FPType, N>& pt) const;
    
    // ElemType kNNValue(const PointND<FPType, N>& key, std::size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTreeExpandLonge*stVec
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const PointND<FPType, N>& key, std::size_t k) const;
    
    // NdTypeOutIt NNsWithFence(const PointND<FPType, N>&, FPType, NdTypeOutIt) const
    // Usage: NdTypeOutIt end = kd.NNsWithFence(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a FPType fence in SQUARE, fill the output iterator result
    // with a collection of all the possible node_type structs whose points as key
    // are at least fence distance closer to p than the rest of the points in the tree,
    // and return the output iterator with final value +1 position.
    // The output iterator is passed in and filled with point-element structs and the end will be returned.
    template <typename NdTypeOutIt> 
        //requires std::output_iterator<NdTypeOutIt, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
    NdTypeOutIt NNsWithFence(const PointND<FPType, N>& p, FPType fenseSq, NdTypeOutIt result) const;
    
private:
    
    struct MetaNode {
        std::uint32_t right_idx;
        std::uint8_t dim_to_expand;
        PointND<FPType, N> key;
        bool operator==(const MetaNode&) const = default;
        bool operator!=(const MetaNode&) const = default;
    };
    
    MetaNode *nd_arr_;
    ElemType *elem_arr_;
    std::uint32_t size_;
    std::uint32_t cap_;
    std::uint32_t height_;
    
    constexpr static std::uint8_t kMaxBalancedTreeHeight = 32;
        
    void DeAlloc();
    
    // ----------------------------------------------------
    // Helper method for range constructor
    //template <std::random_access_iterator RAI>
    template <typename RAI>
    void RangeCtorHelper(RAI, RAI);

    template <class RAI>
    void RangeCtorRecursion(MetaNode*&, ElemType*&, RAI, RAI, std::array<FPType, N>&,
                                  std::array<FPType, N>&, std::array<FPType, N>&);
    
    template <typename RAI>
    void MvConstructOneNdIncIter(MetaNode*& cur_nd, ElemType*& cur_elem,
                                 std::uint32_t right_idx, std::uint8_t dim, RAI nh_it);
    
    template <class ConstRAI>
    static std::tuple<std::array<FPType, N>, std::array<FPType,N>, std::array<FPType, N>>
           ComputeInitBBoxHlSpread(ConstRAI, ConstRAI);
    
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(MetaNode *cur, std::uint8_t dim, const PointND<FPType, N> &pt,
                        BoundedPQueue<ElemType, FPType> &bpq) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const PointND<FPType, N>& key) const;
    
    template <typename def::DistType thisDt = DT,
    typename std::enable_if<thisDt == def::DistType::kEuc, int>::type = 0>
    void NNValueHelper(MetaNode*, std::uint8_t, const PointND<FPType, N>&,
                       const ElemType *&, FPType&) const;
    
    template <typename def::DistType thisDt = DT,
    typename std::enable_if<thisDt != def::DistType::kEuc, int>::type = 0>
    void NNValueHelper(MetaNode*, std::uint8_t, const PointND<FPType, N>&,
                       const ElemType*&, FPType&) const;
    
    
    // MetaNode** findNodePtr(const PointND<FPType, N>& pt);
    // MetaNode*const* findNodePtr(const PointND<FPType, N>& pt) const;
    // Usage: MetaNode **nodePtr = findNodePtr(pt);
    // ----------------------------------------------------
    // Returns the pointer pointing to the node address
    // corresponding to the given PointND. In this FPType pointing
    // fashion, we can construct a node at that location.
    MetaNode** findNodePtr(const PointND<FPType, N>& pt);
    MetaNode*const* findNodePtr(const PointND<FPType, N>& pt) const;
    
    FPType branchMin(const PointND<FPType, N>&, const PointND<FPType, N>&, std::uint8_t) const;
    
};

/** KDTreeExpandLongestVec class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::KDTreeExpandLongestVec() {
    size_ = cap_ = height_ = 0;
    nd_arr_ = nullptr;
    elem_arr_ = nullptr;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs) : size_(rhs.size_), cap_(rhs.cap_), height_(rhs.height_) {
    nd_arr_ = new MetaNode[size_];
    std::copy_n(rhs.nd_arr_, size_, nd_arr_);
    elem_arr_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType), std::nothrow));
    std::uninitialized_copy_n(rhs.elem_arr_, size_, elem_arr_);
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>&
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator=(const KDTreeExpandLongestVec& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        DeAlloc();
        size_ = cap_ = rhs.size_;
        height_ = rhs.height_;
        nd_arr_ = new MetaNode[size_];
        std::copy_n(rhs.nd_arr_, size_, nd_arr_);
        elem_arr_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType)));
        std::uninitialized_copy_n(rhs.elem_arr_, size_, elem_arr_);
    }
    return *this;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
KDTreeExpandLongestVec(KDTreeExpandLongestVec&& rhs) noexcept
: nd_arr_(rhs.nd_arr_), elem_arr_(rhs.elem_arr_), size_(rhs.size_), cap_(rhs.cap_), height_(rhs.height_) {
    rhs.nd_arr_ = nullptr;
    rhs.elem_arr_ = nullptr;
    rhs.size_ = rhs.cap_ = 0;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
operator=(KDTreeExpandLongestVec&& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        DeAlloc();
        nd_arr_ = std::exchange(rhs.nd_arr_, nullptr);
        elem_arr_ = std::exchange(rhs.elem_arr_, nullptr);
        size_ = std::exchange(rhs.size_, 0);
        cap_ = std::exchange(rhs.cap_, 0);
        height_ = rhs.height_;
    }
    return *this;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::~KDTreeExpandLongestVec() {
    DeAlloc();
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <typename RAI>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
MvConstructOneNdIncIter(MetaNode*& cur_nd, ElemType*& cur_elem,
                        std::uint32_t right_idx, std::uint8_t dim, RAI nh_it) {
    *cur_nd++ = {right_idx, dim, std::move(nh_it->key)};
    new (cur_elem++) ElemType (std::move(nh_it->value));
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::DeAlloc() {
    delete[] nd_arr_;
    std::destroy_n(elem_arr_, size_);
    ::operator delete(static_cast<void*>(elem_arr_));
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
//template <std::random_access_iterator ConstRAI> requires def::const_iterator<ConstRAI>
//&& std::same_as<typename std::iter_value_t<ConstRAI>, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
template <typename ConstRAI, std::enable_if_t<std::is_const_v<typename std::remove_pointer_t<typename std::iterator_traits<ConstRAI>::pointer>>, int>>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::KDTreeExpandLongestVec(ConstRAI cbegin, ConstRAI cend) 
    : size_(static_cast<std::uint32_t>(cend - cbegin)), cap_(size_), height_(static_cast<std::uint32_t>(std::log2(size_)+1)) {
    std::vector<node_type> constructData(cbegin, cend);
    RangeCtorHelper(constructData.begin(), constructData.end());
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
//template <std::random_access_iterator RAI> requires def::non_const_iterator<RAI>
//requires std::same_as<typename std::iterator_traits<RAI>::value_type, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
template <typename RAI, std::enable_if_t<!std::is_const_v<typename std::remove_pointer_t<typename std::iterator_traits<RAI>::pointer>>, int>>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::KDTreeExpandLongestVec(RAI begin, RAI end) 
    : size_(static_cast<std::uint32_t>(end - begin)), cap_(size_), height_(static_cast<std::uint32_t>(std::log2(size_)+1)) {
    RangeCtorHelper(begin, end);
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
//template <std::random_access_iterator RAI>
template <typename RAI>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::RangeCtorHelper(RAI data_begin, RAI data_end) {
    if (size_ == 0) [[unlikely]]
        return;
    nd_arr_ = new MetaNode[size_];
    elem_arr_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType), std::nothrow));
    if (size_ == 1) [[unlikely]] {
        *nd_arr_ = {0, N, std::move(data_begin->key)};
        new (elem_arr_) ElemType(std::move(data_begin->value));
        return;
    }
    auto [lows, highs, hl_spread] = ComputeInitBBoxHlSpread(std::as_const(data_begin), std::as_const(data_end));
    MetaNode* cur_nd = nd_arr_;
    ElemType* cur_elem = elem_arr_;
    RangeCtorRecursion(cur_nd, cur_elem, data_begin, data_end, lows, highs, hl_spread);

    /*
    struct ActRecord {
        std::uint8_t dim;
        std::uint32_t end_idx_or_zero_as_traced; // use 0 to mark as traced
        FPType prev_high_low_on_dim;
    };
    std::array<ActRecord, kMaxBalancedTreeHeight> ar_stack;
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator ar_it = ar_stack.begin();
    RAI this_begin = data_begin, this_end = data_end;

    auto RevertOneNode = [&ar_it](auto& highs_or_lows, auto& hl_spread, auto& highs, auto& lows) {
        auto& [dim_on_depth, end_idx_or_traced, prev_high_low_on_dim] = *ar_it;
        end_idx_or_traced = 0;
        highs_or_lows[dim_on_depth] = prev_high_low_on_dim;
        hl_spread[dim_on_depth] = highs[dim_on_depth] - lows[dim_on_depth];
    };

    while (true) {
        std::uint8_t dim = static_cast<std::uint8_t>(std::max_element(hl_spread.cbegin(), hl_spread.cend())
            - hl_spread.cbegin());
        std::ptrdiff_t left_branch_size = (this_end - this_begin) / 2;
        MetaNode* nd_ptr_this_iter = cur_nd;
        RAI median = this_begin + left_branch_size;
        std::nth_element(this_begin, median, this_end,
            [dim](const auto& nh1, const auto& nh2) {return nh1.key[dim] < nh2.key[dim]; });
        MvConstructOneNdIncIter(cur_nd, cur_elem,
            static_cast<std::uint32_t>(cur_nd - nd_arr_ + left_branch_size + 1), dim, median);

        if (left_branch_size == 1) {
            MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, this_begin);
            if (RAI right_child_data_it = median + 1;
                right_child_data_it != this_end) {
                MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, right_child_data_it);
            } else if (nd_ptr_this_iter->right_idx = 0;
                (ar_it - 1)->end_idx_or_zero_as_traced
                - static_cast<std::uint32_t>(this_end - data_begin) == 2) {
                // single bottom right child up one level
                MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, ++this_end);
                ++this_end;
                // traceback: directly revert one left instead of set one right - revert one right
                --ar_it;
                RevertOneNode(highs, hl_spread, highs, lows);
            }
            // termination
            if (this_end == data_end) [[unlikely]]
                break;
            // traceback: revert right path
            while ((--ar_it)->end_idx_or_zero_as_traced == 0)
                RevertOneNode(lows, hl_spread, highs, lows);
            // traceback: revert one left set one right on ancestor node
            auto& [dim_on_depth, end_idx_or_traced, prev_high_low_on_dim] = *ar_it++;
            utility::CycleSwap(prev_high_low_on_dim, lows[dim_on_depth], highs[dim_on_depth]);
            hl_spread[dim_on_depth] = highs[dim_on_depth] - lows[dim_on_depth];
            this_begin = this_end + 1;
            this_end = data_begin + end_idx_or_traced;
            end_idx_or_traced = 0;
        } else {
            *ar_it++ = {dim, static_cast<std::uint32_t>(this_end - data_begin), highs[dim]};
            highs[dim] = nd_ptr_this_iter->key[dim];
            hl_spread[dim] = highs[dim] - lows[dim];
            this_end = median;
        }
    } */

    /*
    struct ActRecord {
        std::uint8_t dim;
        std::uint32_t end_idx_or_zero_as_traced; // use 0 to mark as traced
        FPType prev_high_low_on_dim;
    };
    std::array<ActRecord, kMaxBalancedTreeHeight> ar_stack;
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator ar_it = ar_stack.begin();
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator leaf_prt_ar_it = ar_stack.begin() + (height_ > 2 ? height_ - 2 : 0);
    RAI this_begin = data_begin, this_end = data_end;

    auto RevertOneNode = [&ar_it](auto& highs_or_lows, auto& hl_spread, auto& highs, auto& lows) {
        auto &[dim_on_depth, end_idx_or_traced, prev_high_low_on_dim] = *ar_it;
        highs_or_lows[dim_on_depth] = prev_high_low_on_dim;
        hl_spread[dim_on_depth] = highs[dim_on_depth] - lows[dim_on_depth];
    };
    
    std::uint32_t leaf_prt_idx = 0;
    while (true) {
        std::ptrdiff_t left_subtree_size = (this_end - this_begin)/2;
        RAI median = this_begin + left_subtree_size;
        for (; ar_it != leaf_prt_ar_it; ++ar_it) {
            std::uint8_t dim = static_cast<std::uint8_t>(std::max_element(hl_spread.cbegin(), hl_spread.cend())
                - hl_spread.cbegin());
            std::nth_element(this_begin, median, this_end,
                [dim](const auto& nh1, const auto& nh2) {return nh1.key[dim] < nh2.key[dim]; });
            MvConstructOneNdIncIter(cur_nd, cur_elem,
                static_cast<std::uint32_t>(cur_nd - nd_arr_ + left_subtree_size + 1), dim, median);
            *ar_it = {dim, static_cast<std::uint32_t>(this_end - data_begin), std::exchange(highs[dim], (cur_nd-1)->key[dim])};
            hl_spread[dim] = highs[dim] - lows[dim];
            left_subtree_size /= 2;
            this_end = std::exchange(median, this_begin + left_subtree_size);
        }
        
        if (left_subtree_size == 0) {
            MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, median);
        } else {
            std::uint8_t dim = static_cast<std::uint8_t>(std::max_element(hl_spread.cbegin(), hl_spread.cend())
                - hl_spread.cbegin());
            if (median + 1 == this_end) {
                if (this_begin->key[dim] >= median->key[dim])
                    std::swap(this_begin, median);
                MvConstructOneNdIncIter(cur_nd, cur_elem, 0, dim, median);
                MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, this_begin);
            } else {
                utility::Sort3(this_begin, median, median+1, [dim](const auto& nh1, const auto& nh2) {return nh1.key[dim] < nh2.key[dim]; });
                MvConstructOneNdIncIter(cur_nd, cur_elem, static_cast<std::uint32_t>(cur_nd - nd_arr_ + 2), dim, median);
                MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, this_begin);
                MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, median+1);
            }
        }
        
        if (this_end == data_end) [[unlikely]]
            return;

        for (std::uint32_t i = std::countr_zero(++leaf_prt_idx); i != 0; --i) {
            --ar_it;
            RevertOneNode(lows, hl_spread, highs, lows);
        }

        auto& [dim_on_depth, end_idx_or_traced, prev_high_low_on_dim] = *(--ar_it);
        utility::CycleSwap(prev_high_low_on_dim, lows[dim_on_depth], highs[dim_on_depth]);
        hl_spread[dim_on_depth] = highs[dim_on_depth] - lows[dim_on_depth];
        this_begin = this_end + 1;
        this_end = data_begin + end_idx_or_traced;
        ++ar_it; 
    } */
/*
DEBUG_PRINT:
    std::for_each_n(nd_arr_, size_, [](const MetaNode& nd){
        std::cout << "right_idx: " << nd.right_idx << '\t' << "dim_to_expand: " << static_cast<std::uint32_t>(nd.dim_to_expand)  << '\t';
        std::cout << "pt coord: ";
        std::copy(nd.key.cbegin(), nd.key.cend(), std::ostream_iterator<FPType>(std::cout,","));
        std::cout << '\n';
    }); */
}



template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <class ConstRAI>
std::tuple<std::array<FPType, N>, std::array<FPType, N>, std::array<FPType, N>> KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
ComputeInitBBoxHlSpread(ConstRAI cbegin, ConstRAI cend) {
    std::array<FPType, N> lows, highs, hl_spread;
    lows.fill(std::numeric_limits<FPType>::max());
    highs.fill(std::numeric_limits<FPType>::min());
    std::for_each(cbegin, cend, [&lows, &highs](const auto &nh) mutable {
        for (std::uint8_t i = 0; i < N; ++i) {
            FPType pt_val_on_dim = nh.key[i];
            auto &bbox_low = lows[i], &bbox_high = highs[i];
            bbox_low = std::min(bbox_low, pt_val_on_dim);
            bbox_high = std::max(bbox_high, pt_val_on_dim);
        }
    });
    std::transform(highs.cbegin(), highs.cend(), lows.cbegin(), hl_spread.begin(), std::minus<FPType>());
    return {lows, highs, hl_spread};
}



template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <class RAI>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
RangeCtorRecursion(MetaNode *&cur_nd, ElemType *&cur_elem, RAI begin, RAI end,
                   std::array<FPType, N> &lows, std::array<FPType, N> &highs, std::array<FPType, N> &hl_spread) {
    std::uint8_t dim = static_cast<std::uint8_t>(std::distance(hl_spread.cbegin(), std::max_element(hl_spread.cbegin(), hl_spread.cend())));
    MetaNode* nd_ptr_this_iter = cur_nd;
    std::ptrdiff_t left_branch_size = (end - begin)/2;
    RAI median = begin + left_branch_size;
    std::nth_element(begin, median, end, [dim](const auto& nh1, const auto& nh2) {return nh1.key[dim] < nh2.key[dim];});
    MvConstructOneNdIncIter(cur_nd, cur_elem, static_cast<std::uint32_t>(cur_nd - nd_arr_ + left_branch_size + 1), dim, median);
    
    if (left_branch_size == 1) {
        MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, begin);
        if (median + 1 != end) {
            MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, median + 1);
        } else {
            nd_ptr_this_iter->right_idx = 0;
        }
        return;
    } else {
        FPType prev_high_on_dim = std::exchange(highs[dim], nd_ptr_this_iter->key[dim]);
        hl_spread[dim] = highs[dim] - lows[dim];
        RangeCtorRecursion(cur_nd, cur_elem, begin, median, lows, highs, hl_spread);
        highs[dim] = prev_high_on_dim;
        hl_spread[dim] = highs[dim] - lows[dim];
    }
    
    if (median + 2 == end) {
        MvConstructOneNdIncIter(cur_nd, cur_elem, 0, N, median + 1);
    } else [[likely]] {
        FPType prev_low_on_dim = std::exchange(lows[dim], nd_ptr_this_iter->key[dim]);
        hl_spread[dim] = highs[dim] - lows[dim];
        RangeCtorRecursion(cur_nd, cur_elem, median+1, end, lows, highs, hl_spread);
        lows[dim] = prev_low_on_dim;
        hl_spread[dim] = highs[dim] - prev_low_on_dim;
    }
    /*
    std::uint8_t dim = static_cast<std::uint8_t>(std::distance(hl_spread.cbegin(), std::max_element(hl_spread.cbegin(), hl_spread.cend())));
    MetaNode* nd_ptr_this_iter = cur_nd;
    std::ptrdiff_t left_branch_size = (end - begin)/2;
    RAI median = begin + left_branch_size;
    std::nth_element(begin, median, end, [dim](const auto& p1, const auto& p2) {return p1.key[dim] < p2.key[dim];});
    //oneapi::dpl::nth_element(oneapi::dpl::execution::par_unseq, begin, median, end, [dim](const auto& p1, const auto& p2) {return p1.key[dim] < p2.key[dim]; });

    *cur_nd++ = {static_cast<std::uint32_t>(cur_nd - nd_arr_ + left_branch_size + 1), dim, std::move(median->key)};
    new (cur_elem++) ElemType (std::move(median->value));
    
    if (left_branch_size == 1) {
        *cur_nd++ = {0, N, std::move(begin->key)};
        new (cur_elem++) ElemType (std::move(begin->value));
        if (median + 1 == end) {
            nd_ptr_this_iter->right_idx = 0;
            return;
        }
    } else {
        FPType prev_high_on_dim = std::exchange(highs[dim], nd_ptr_this_iter->key[dim]);
        hl_spread[dim] = highs[dim] - lows[dim];
        RangeCtorRecursion(cur_nd, cur_elem, begin, median, lows, highs, hl_spread);
        highs[dim] = prev_high_on_dim;
        hl_spread[dim] = prev_high_on_dim - lows[dim];
    }
    
    if (std::distance(median, end) == 2) {
        *cur_nd++ = {0, N, std::move((median+1)->key)};
        new (cur_elem++) ElemType (std::move((median+1)->value));
    } else [[likely]] {
        FPType prev_low_on_dim = std::exchange(lows[dim], nd_ptr_this_iter->key[dim]);
        hl_spread[dim] = highs[dim] - lows[dim];
        RangeCtorRecursion(cur_nd, cur_elem, median+1, end, lows, highs, hl_spread);
        lows[dim] = prev_low_on_dim;
        hl_spread[dim] = highs[dim] - prev_low_on_dim;
    } */
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
constexpr std::uint8_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::dimension() const {
    return N;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
typename def::DistType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::distType() const {
    return DT;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
std::uint32_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::size() const {
    return size_;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
std::uint32_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::cap() const {
    return cap_;
}


template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
std::uint32_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::height() const {
    return height_;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::Empty() const {
    return size_ == 0;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator==(const KDTreeExpandLongestVec& rhs) const {
    if (size_ != rhs.size_)
        return false;
    return std::equal(elem_arr_, elem_arr_ + size_, rhs.elem_arr_) && std::equal(nd_arr_, nd_arr_ + size_, rhs.nd_arr_);
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator!=(const KDTreeExpandLongestVec& rhs) const {
    return !(*this == rhs);
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::PrintTreeInfo() const {
    std::cout << "Tree height is " << height() << "\n"
              << "Tree size is " << size() << "\n"
              << "Tree capacity is " << cap_ << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::Clear() {
    std::destroy_n(nd_arr_, size_);
    std::destroy_n(elem_arr_, size_);
    size_ = 0;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::ShrinkToFit() {
    if (cap_ > size_) {
        // TO DO
    }
    std::destroy_n(nd_arr_, size_);
    std::destroy_n(elem_arr_, size_);
    size_ = 0;
}


template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
insert(const PointND<FPType, N>& pt, const ElemType& value) {
    MetaNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        //(*ndPtr)->object = value;
    } else {
        //*ndPtr = new MetaNode(0, pt, value);
        //treeSize++;
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::contains(const PointND<FPType, N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator[] (const PointND<FPType, N>& pt) {
    MetaNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new MetaNode(pt, ElemType());
        //treeSize++;
    }
    return (*ndPtr)->object;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeExpandLongestVec&>(*this).at(pt));
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
const ElemType& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) const {
    MetaNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::MetaNode**
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) {
    return const_cast<MetaNode**>(static_cast<const KDTreeExpandLongestVec*>(this)->findNodePtr(pt));
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::MetaNode*const*
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) const {
    //MetaNode *const*n = &root;
    MetaNode *const*n;
    for (std::uint8_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1){}
    //n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
kNNValue(const PointND<FPType, N>& pt, std::size_t k) const {
    //if (empty())
    //    return ElemType();
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType, FPType> bpq(k);
    //kNNValueHelper(ndVec[0], 0, pt, bpq);
    
    std::multimap<std::uint32_t, ElemType, std::greater<std::size_t>> freqMap;
    while (!bpq.empty()) {
        ElemType elem = bpq.dequeueMin();
        for (typename std::multimap<std::uint32_t, ElemType>::iterator
             it = freqMap.begin();it != freqMap.end(); ++it) {
            if (it->second == elem) {
                freqMap.emplace(it->first+1, it->second);
                freqMap.erase(it);
                goto outer;
            }
        }
        freqMap.emplace(1, elem);
    outer:
        ;
    }
    return freqMap.begin()->second;
    
    /*
     -- This version requires ElemType to have its own hash function for best
     performance. But for the k size to be usually small, The difference
     between O(1) and O(lgk) and O(k) are hard to tell. So we use the
     uncommented version which does not require hash or comparison operators.
     
     std::unordered_map<ElemType, int> freqMap;
     while (!bpq.empty()) {
     freqMap[bpq.dequeueMin()]++;
     }
     int max = 0;
     ElemType frequent = ElemType();
     for (const std::pair<ElemType, int> &p : freqMap) {
     if (p.second > max) {
     max = p.second;
     frequent = p.first;
     }
     }
     return frequent; */
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::kNNValueHelper(MetaNode *cur, std::uint8_t dim,
                                                             const PointND<FPType, N>& pt, BoundedPQueue<ElemType, FPType> &bpq) const {
    bpq.enqueue(cur->object, PointND<FPType, N>::template dist<DT>(cur->key, pt));
    std::uint8_t next_dim = dim + 1 < N ? dim + 1 : 0;
    MetaNode *next = pt[dim] < cur->key[dim] ? cur->left : cur->right;
    if (next)
        kNNValueHelper(next, next_dim, pt, bpq);
    if (bpq.size() < bpq.maxSize()
        || branchMin(cur->key, pt, dim) < bpq.worst()) {
        MetaNode *other = next == cur->left ? cur->right : cur->left;
        if (other)
            kNNValueHelper(other, next_dim, pt, bpq);
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <typename NdTypeOutIt> 
    //requires std::output_iterator<NdTypeOutIt, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
NdTypeOutIt KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
NNsWithFence(const PointND<FPType, N>& pt, FPType fence_sq, NdTypeOutIt pe_out_it) const {
    struct ActRecord {
        FPType diff_on_dim_sq;
        std::uint32_t on_stack_idx;
    } ar_stack[kMaxBalancedTreeHeight], *ar_it = ar_stack;
    const MetaNode *cur_nd = nd_arr_;
    FPType cur_dist_sq, best_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(cur_nd->key, pt);
    FPType best_dist_plus_fence_sq = best_dist_sq + fence_sq + FPType(2.0)*sqrt(fence_sq*best_dist_sq);
    
    constexpr std::uint_fast16_t kMaxDpeSizeOnStack = 1152;
    struct DistPtElem {
        FPType dist_sq;
        const PointND<FPType, N>* pt;
        const ElemType* elem;
    };
    std::array<DistPtElem, kMaxDpeSizeOnStack> result_dpe_arr;
    typename std::array<DistPtElem, kMaxDpeSizeOnStack>::iterator result_dpe_arr_it = result_dpe_arr.begin();
    *result_dpe_arr_it++ = {best_dist_sq, &nd_arr_->key, elem_arr_};
    std::vector<DistPtElem> result_dpe_vec;
   
    auto CheckOneNd = [&](const auto& nd) {
        cur_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(nd->key, pt);
        if (cur_dist_sq < best_dist_plus_fence_sq) {
            if (cur_dist_sq < best_dist_sq) {
                best_dist_sq = cur_dist_sq;
                best_dist_plus_fence_sq = best_dist_sq + fence_sq + FPType(2.0)*sqrt(fence_sq*best_dist_sq);
            }
            *result_dpe_arr_it++ = {cur_dist_sq, &nd->key, &elem_arr_[nd-nd_arr_]};
            if (result_dpe_arr_it == result_dpe_arr.end()) [[unlikely]] {
                result_dpe_arr_it = result_dpe_arr.begin();
                result_dpe_vec.insert(result_dpe_vec.end(), result_dpe_arr.begin(), result_dpe_arr.end());
            }
        }
    };

    std::uint32_t cur_nd_idx = 0, right_idx;
    while (true) {
        if ((right_idx = cur_nd->right_idx)) {
            ++cur_nd_idx;
            FPType diff = pt[cur_nd->dim_to_expand] - cur_nd->key[cur_nd->dim_to_expand];
            if (diff < 0.0) {
                *ar_it = {diff*diff, right_idx};
                ++cur_nd;
            } else {
                *ar_it = {diff*diff, cur_nd_idx};
                cur_nd_idx = right_idx;
                cur_nd = nd_arr_ + right_idx;
            }
            ++ar_it;
            CheckOneNd(cur_nd); 
        } else {
            if (cur_nd->dim_to_expand != N)
                CheckOneNd(++cur_nd);
            do {
                if (ar_it == ar_stack) [[unlikely]] {
                    auto filter_transform_dpe_arr_to_pe_arr = [&pe_out_it, best_dist_plus_fence_sq](const auto &dpe) mutable {
                        if (dpe.dist_sq < best_dist_plus_fence_sq)
                            *pe_out_it++ = { *(dpe.pt), *(dpe.elem) };
                    };
                    std::for_each(std::make_reverse_iterator(result_dpe_arr_it), result_dpe_arr.rend(),
                                  filter_transform_dpe_arr_to_pe_arr);
                    std::for_each(result_dpe_vec.rbegin(), result_dpe_vec.rend(), filter_transform_dpe_arr_to_pe_arr); 
                    return pe_out_it;
                }
            } while ((--ar_it)->diff_on_dim_sq >= best_dist_plus_fence_sq);
            cur_nd = nd_arr_ + (cur_nd_idx = ar_it->on_stack_idx);
            CheckOneNd(cur_nd);
        }
    }
    return pe_out_it;
}

/*
using PtNodePtr = llama::Record <
        llama::Field<ns::Point, llama::Array<FPType, N>>,
        llama::Field<ns::NodePtr, MetaNode*>
    >;
llama::ArrayDims<N> pt_ndPtr_arr_size;
std::fill(pt_ndPtr_arr_size.begin(), pt_ndPtr_arr_size.end(), kMaxBalancedTreeHeight);
auto pt_ndPtr_soa_view = llama::allocView(llama::mapping::SoA<llama::ArrayDims<N>, PtNodePtr>(pt_ndPtr_arr_size));
pt_ndPtr_soa_view(0ull, 0ull, 0ull)(ns::NodePtr{}) = nd_arr_;
*/

//static std::size_t search_times = 0, check_count = 0;
/*
template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N> &search_pt) const {
    struct ActRecord {
        FPType diff_on_dim_sq;
        std::uint32_t on_stack_idx;
    };
    std::array<ActRecord, kMaxBalancedTreeHeight> ar_stack;
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator
        ar_it = ar_stack.begin(), leaf_prt_ar_it = ar_stack.begin() + (height_ > 2 ? height_ - 2 : 0);
    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    const MetaNode *cur_nd = nd_arr_, *best_node = nd_arr_;
    FPType best_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(cur_nd->key, search_pt);
    
    auto CheckOneNd = [&](const auto& nd) {
        if (FPType cur_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(nd->key, search_pt);
            cur_dist_sq < best_dist_sq) {
            best_dist_sq = cur_dist_sq;
            best_node = nd;
        }
        //check_count++;
    }; */
    /*
    search_times++;
    if (search_times > 22000) {
        std::cout << check_count << '\n';
        assert(false);
    }*/

/*
    while (true) {
        if (std::uint32_t right_idx = cur_nd->right_idx) {
            std::uint8_t dim = cur_nd->dim_to_expand;
            FPType diff = search_pt[dim] - cur_nd->key[dim];
            if (diff < 0.0) {
                *ar_it++ = {diff*diff, right_idx};
                ++cur_nd;
            } else {
                *ar_it++ = {diff*diff, static_cast<std::uint32_t>(cur_nd-nd_arr_+1)};
                cur_nd = nd_arr_ + right_idx;
            }
            CheckOneNd(cur_nd);
        } else {
            if (cur_nd->dim_to_expand != N)
                CheckOneNd(++cur_nd);
            do {
                if (ar_it == ar_stack.begin())
                    return elem_arr_[best_node - nd_arr_];
            } while ((--ar_it)->diff_on_dim_sq >= best_dist_sq);
            cur_nd = nd_arr_ + ar_it->on_stack_idx;
            CheckOneNd(cur_nd);
        }
    }
    return elem_arr_[best_node - nd_arr_];
*/
    /*
    while (true) {
        for (; ar_it != leaf_prt_ar_it; ++ar_it) {
            auto& [right_idx, dim, nd_pt] = *cur_nd;
            FPType diff = search_pt[dim] - nd_pt[dim];
            if (diff < 0.0) {
                *ar_it = {diff*diff, right_idx};
                ++cur_nd;
            } else {
                *ar_it = {diff*diff, static_cast<std::uint32_t>(cur_nd-nd_arr_+1)};
                cur_nd = nd_arr_ + right_idx;
            }
            CheckOneNd(cur_nd);
        }
        if (cur_nd->right_idx) {
            CheckOneNd(++cur_nd);
            CheckOneNd(++cur_nd);
        } else if (cur_nd->dim_to_expand != N) {
            CheckOneNd(++cur_nd);
        }
        do {
            if (ar_it == ar_stack.begin()) [[unlikely]]
                return elem_arr_[best_node - nd_arr_];
            auto& [diff_on_dim_sq, on_stack_idx] = *--ar_it;
            if (diff_on_dim_sq < best_dist_sq) {
                cur_nd = nd_arr_ + on_stack_idx;
                CheckOneNd(cur_nd);
                if (cur_nd->right_idx != 0) {
                    if ((cur_nd + 1)->dim_to_expand != N) {
                        break;
                    } else {
                        CheckOneNd(++cur_nd);
                        CheckOneNd(++cur_nd);
                        return elem_arr_[best_node - nd_arr_];
                    }
                } else if (cur_nd->dim_to_expand != N) {
                    CheckOneNd(++cur_nd);
                }
            }
        } while (true);
    }
    return elem_arr_[best_node - nd_arr_];
*/ /*
} */

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N>& search_pt) const {
    
    const MetaNode* cur_nd = nd_arr_;
    FPType best_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(cur_nd->key, search_pt);
    std::uint32_t best_nd_idx = 0, cur_nd_idx = 0, right_idx;
    struct ActRecord {
        FPType diff_on_dim_sq;
        std::uint32_t on_stack_idx;
    };
    std::array<ActRecord, kMaxBalancedTreeHeight> ar_stack;
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator ar_it = ar_stack.begin();
    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    
    auto CheckOneNd = [&](const auto& nd, const std::uint32_t& cur_nd_idx) {
        if (FPType cur_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(nd->key, search_pt);
            cur_dist_sq < best_dist_sq) {
            best_dist_sq = cur_dist_sq;
            best_nd_idx = cur_nd_idx;
        }
    };
    while (true) {
        if ((right_idx = cur_nd->right_idx)) {
            ++cur_nd_idx;
            FPType diff = search_pt[cur_nd->dim_to_expand] - cur_nd->key[cur_nd->dim_to_expand];
            if (diff < 0.0) {
                *ar_it = {diff*diff, right_idx};
                ++cur_nd;
            } else {
                *ar_it = {diff*diff, cur_nd_idx};
                cur_nd_idx = right_idx;
                cur_nd = nd_arr_ + right_idx;
            }
            ++ar_it;
            CheckOneNd(cur_nd, cur_nd_idx);
        } else {
            if (cur_nd->dim_to_expand != N)
                CheckOneNd(++cur_nd, ++cur_nd_idx);
            do {
                if (ar_it == ar_stack.begin())
                    return elem_arr_[best_nd_idx];
            } while ((--ar_it)->diff_on_dim_sq >= best_dist_sq);
            cur_nd = nd_arr_ + (cur_nd_idx = ar_it->on_stack_idx);
            CheckOneNd(cur_nd, cur_nd_idx);
        }
    }
    return elem_arr_[best_nd_idx];
}

/*
template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N>& search_pt) const {
    struct ActRecord {
        FPType diff_on_dim_sq;
        std::uint32_t on_stack_idx;
        std::uint32_t cand_nd_idx;
    };
    std::array<ActRecord, kMaxBalancedTreeHeight> ar_stack;
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator
        ar_it = ar_stack.begin() + 1, cand_nd_ar_it, prev_ar_it = ar_stack.begin();

    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    const MetaNode* cur_nd = nd_arr_;
    FPType best_dist_sq = std::numeric_limits<FPType>::max();
    prev_ar_it->cand_nd_idx = 0;
    std::uint32_t best_nd_idx = 0;
    while (true) {
        if (std::uint32_t right_idx = cur_nd->right_idx) {
            std::uint8_t dim = cur_nd->dim_to_expand;
            FPType diff = search_pt[dim] - cur_nd->key[dim];
            if (diff < 0.0) {
                *ar_it++ = {diff*diff, right_idx, static_cast<std::uint32_t>(++cur_nd-nd_arr_)};
            } else {
                *ar_it++ = {diff*diff, static_cast<std::uint32_t>(cur_nd-nd_arr_+1), right_idx};
                cur_nd = nd_arr_ + right_idx;
            }
        } else { 
            cand_nd_ar_it = ar_it;
            if (cur_nd->dim_to_expand != N)
                cand_nd_ar_it++->cand_nd_idx = static_cast<std::uint32_t>(++cur_nd - nd_arr_);
            std::for_each(std::make_reverse_iterator(cand_nd_ar_it), std::make_reverse_iterator(prev_ar_it), [&](const auto& ar) {
                if (FPType cur_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>((nd_arr_ + ar.cand_nd_idx)->key, search_pt);
                    cur_dist_sq < best_dist_sq) {
                    best_dist_sq = cur_dist_sq;
                    best_nd_idx = ar.cand_nd_idx;
                }
            });
            do {
                if (ar_it == ar_stack.begin() + 1)
                    return elem_arr_[best_nd_idx];
            } while ((--ar_it)->diff_on_dim_sq >= best_dist_sq);
            prev_ar_it = ar_it - 1;
            prev_ar_it->cand_nd_idx = ar_it->on_stack_idx;
            cur_nd = nd_arr_ + ar_it->on_stack_idx;
        }
    }
    return elem_arr_[best_nd_idx];
}
*/

/*
template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N>& search_pt) const {
    struct ActRecord {
        FPType diff_on_dim_sq;
        std::uint32_t on_stack_nd_idx;
        std::uint32_t cand_nd_idx;
    };
    std::array<ActRecord, kMaxBalancedTreeHeight> ar_stack;
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator ar_it = ar_stack.begin();
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator cand_nd_ar_it = ar_stack.begin();
    constexpr std::uint32_t num_lanes = 32/sizeof(FPType);

    Eigen::Matrix<FPType, num_lanes, N> cand_pts_matrix;
    std::uint32_t lane_idx = 0;
    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    const MetaNode* cur_nd = nd_arr_;
    std::uint32_t best_idx = 0;
    FPType best_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(cur_nd->key, search_pt);
    Eigen::Matrix<FPType, 1, N> search_pt_v = Eigen::Map<Eigen::Matrix<FPType, 1, N>>(const_cast<FPType*>(search_pt.DataArray().data()));
    while (true) {
        if (std::uint32_t right_idx = cur_nd->right_idx) {
            std::uint8_t dim = cur_nd->dim_to_expand;
            FPType diff = search_pt[dim] - cur_nd->key[dim];
            ar_it->diff_on_dim_sq = diff * diff;
            if (diff < 0.0) {
                ar_it++->on_stack_nd_idx = right_idx;
                ++cur_nd;
            } else {
                ar_it++->on_stack_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_ + 1);
                cur_nd = nd_arr_ + right_idx;
            }
        } else if (cur_nd++->dim_to_expand == N) {
            if (lane_idx != 0) {
                if (FPType best_dist_sq_this_iter =
                    (cand_pts_matrix(Eigen::seq(0, lane_idx - 1), Eigen::all).rowwise() - search_pt_v)
                    .rowwise().squaredNorm().minCoeff(&lane_idx);
                    best_dist_sq_this_iter < best_dist_sq) {
                    best_dist_sq = best_dist_sq_this_iter;
                    best_idx = ar_stack[lane_idx].cand_nd_idx;
                }
                lane_idx = 0;
                cand_nd_ar_it = ar_stack.begin();
            }
            do {
                if (ar_it == ar_stack.begin()) [[unlikely]]
                    return *(elem_arr_ + best_idx);
            } while ((--ar_it)->diff_on_dim_sq >= best_dist_sq);
            cur_nd = nd_arr_ + ar_it->on_stack_nd_idx;
        }
        cand_pts_matrix.row(lane_idx++) = Eigen::Map<Eigen::Matrix<FPType, 1, N>>(const_cast<FPType*>(cur_nd->key.DataArray().data()));
        cand_nd_ar_it++->cand_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_);
        if (lane_idx == num_lanes) {
            if (FPType best_dist_sq_this_iter = (cand_pts_matrix.rowwise() - search_pt_v)
                .rowwise().squaredNorm().minCoeff(&lane_idx);
                best_dist_sq_this_iter < best_dist_sq) {
                best_dist_sq = best_dist_sq_this_iter;
                best_idx = ar_stack[lane_idx].cand_nd_idx;
            }
            lane_idx = 0;
            cand_nd_ar_it = ar_stack.begin();
        }
    }
    return *(elem_arr_ + best_idx);
} */


/*
template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N>& search_pt) const {
    struct ActRecord {
        FPType diff_on_dim_sq;
        std::uint32_t on_stack_idx;
        std::uint32_t cand_nd_idx;
    };
    std::array<ActRecord, kMaxBalancedTreeHeight> ar_stack;
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator
        ar_it = ar_stack.begin() + 1, cand_nd_ar_it, prev_ar_it = ar_stack.begin();

    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    const MetaNode* cur_nd = nd_arr_;
    FPType best_dist_sq = std::numeric_limits<FPType>::max();
    prev_ar_it->cand_nd_idx = 0;
    std::uint32_t best_nd_idx = 0;
    
    Eigen::Matrix<FPType, 1, N> search_pt_v = Eigen::Map<Eigen::Matrix<FPType, 1, N>>(const_cast<FPType*>(search_pt.DataArray().data()));
    Eigen::Matrix<FPType, kMaxBalancedTreeHeight, N> cand_pts_matrix;
    std::uint32_t cand_pts_matrix_idx = 0;
    
    while (true) {
        if (std::uint32_t right_idx = cur_nd->right_idx) {
            std::uint8_t dim = cur_nd->dim_to_expand;
            FPType diff = search_pt[dim] - cur_nd->key[dim];
            if (diff < 0.0) {
                *ar_it++ = {diff*diff, right_idx, static_cast<std::uint32_t>(++cur_nd-nd_arr_)};
            } else {
                *ar_it++ = {diff*diff, static_cast<std::uint32_t>(cur_nd-nd_arr_+1), right_idx};
                cur_nd = nd_arr_ + right_idx;
            }
        } else {
            cand_nd_ar_it = ar_it;
            if (cur_nd->dim_to_expand != N)
                cand_nd_ar_it++->cand_nd_idx = static_cast<std::uint32_t>(++cur_nd - nd_arr_);
            std::for_each(std::make_reverse_iterator(cand_nd_ar_it), std::make_reverse_iterator(prev_ar_it), [&](const auto& ar) {
                cand_pts_matrix.row(cand_pts_matrix_idx++) = Eigen::Map<Eigen::Matrix<FPType, 1, N>>(const_cast<FPType*>((nd_arr_ + ar.cand_nd_idx)->key.DataArray().data()));
            });
            if (FPType best_dist_sq_this_iter =
                (cand_pts_matrix(Eigen::seq(0, cand_pts_matrix_idx - 1), Eigen::all).rowwise() - search_pt_v)
                .rowwise().squaredNorm().minCoeff(&cand_pts_matrix_idx);
                best_dist_sq_this_iter < best_dist_sq) {
                best_dist_sq = best_dist_sq_this_iter;
                best_nd_idx = (std::make_reverse_iterator(cand_nd_ar_it) + cand_pts_matrix_idx)->cand_nd_idx;
            }
            cand_pts_matrix_idx = 0;
            do {
                if (ar_it == ar_stack.begin() + 1)
                    return elem_arr_[best_nd_idx];
            } while ((--ar_it)->diff_on_dim_sq >= best_dist_sq);
            prev_ar_it = ar_it - 1;
            prev_ar_it->cand_nd_idx = ar_it->on_stack_idx;
            cur_nd = nd_arr_ + ar_it->on_stack_idx;
        }
    }
    return elem_arr_[best_nd_idx];
}
*/

/*
template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N> &search_pt) const {
    struct ActRecord {
        FPType diff_on_dim_sq;
        std::uint32_t on_stack_nd_idx;
        std::uint32_t cand_nd_idx;
    };
    std::array<ActRecord, kMaxBalancedTreeHeight> ar_stack;
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator ar_it = ar_stack.begin();
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator cand_nd_ar_it = ar_stack.begin();
    
    Eigen::Matrix<FPType, kMaxBalancedTreeHeight, N> cand_pts_matrix;
    std::uint32_t cand_pts_matrix_idx = 0;

    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    const MetaNode *cur_nd = nd_arr_;
    std::uint32_t best_idx = 0;
    FPType best_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(cur_nd->key, search_pt);
    Eigen::Matrix<FPType, 1, N> search_pt_v = Eigen::Map<Eigen::Matrix<FPType, 1, N>>(const_cast<FPType*>(search_pt.DataArray().data()));
    while (true) {
        if (std::uint32_t right_idx = cur_nd->right_idx) {
            std::uint8_t dim = cur_nd->dim_to_expand;
            FPType diff = search_pt[dim] - cur_nd->key[dim];
            ar_it->diff_on_dim_sq = diff * diff;
            if (diff < 0.0) {
                ar_it++->on_stack_nd_idx = right_idx;
                ++cur_nd;
            } else {
                ar_it++->on_stack_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_ + 1);
                cur_nd = nd_arr_ + right_idx;
            }
            cand_pts_matrix.row(cand_pts_matrix_idx++) = Eigen::Map<Eigen::Matrix<FPType, 1, N>>(const_cast<FPType*>(cur_nd->key.DataArray().data()));
            cand_nd_ar_it++->cand_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_);
        } else {
            if (cur_nd++->dim_to_expand != N) {
                cand_pts_matrix.row(cand_pts_matrix_idx++) = Eigen::Map<Eigen::Matrix<FPType, 1, N>>(const_cast<FPType*>(cur_nd->key.DataArray().data()));
                cand_nd_ar_it++->cand_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_);
            }
            if (FPType best_dist_sq_this_iter = 
                (cand_pts_matrix(Eigen::seq(0, cand_pts_matrix_idx - 1), Eigen::all).rowwise() - search_pt_v)
                .rowwise().squaredNorm().minCoeff(&cand_pts_matrix_idx);
                best_dist_sq_this_iter < best_dist_sq) {
                best_dist_sq = best_dist_sq_this_iter;
                best_idx = ar_stack[cand_pts_matrix_idx].cand_nd_idx;
            }
            cand_pts_matrix_idx = 0;
            cand_nd_ar_it = ar_stack.begin();
            do {
                if (ar_it == ar_stack.begin()) [[unlikely]]
                    return *(elem_arr_ + best_idx);
            } while ((--ar_it)->diff_on_dim_sq >= best_dist_sq);
            cur_nd = nd_arr_ + ar_it->on_stack_nd_idx;
            cand_pts_matrix.row(cand_pts_matrix_idx++) = Eigen::Map<Eigen::Matrix<FPType, 1, N>>(const_cast<FPType*>(cur_nd->key.DataArray().data()));
            cand_nd_ar_it++->cand_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_);
        }
    }
    return *(elem_arr_ + best_idx);
}
*/


/*
template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N>& search_pt) const {
    struct ActRecord {
        FPType diff_on_dim_sq;
        std::uint32_t on_stack_nd_idx;
        std::uint32_t cand_nd_idx;
    };
    std::array<ActRecord, kMaxBalancedTreeHeight> ar_stack;
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator ar_it = ar_stack.begin();
    typename std::array<ActRecord, kMaxBalancedTreeHeight>::iterator cand_nd_ar_it = ar_stack.begin();

    std::array<FPType, N*kMaxBalancedTreeHeight> cand_pts_matrix;
    std::int32_t cand_pts_matrix_idx = 0;

    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    const MetaNode* cur_nd = nd_arr_;
    std::uint32_t best_idx = 0;
    FPType best_dist_sq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(cur_nd->key, search_pt);
    Eigen::Matrix<FPType, 1, N> search_pt_v = Eigen::Map<Eigen::Matrix<FPType, 1, N>>(const_cast<FPType*>(search_pt.DataArray().data()));
    while (true) {
        if (std::uint32_t right_idx = cur_nd->right_idx) {
            std::uint8_t dim = cur_nd->dim_to_expand;
            FPType diff = search_pt[dim] - cur_nd->key[dim];
            ar_it->diff_on_dim_sq = diff * diff;
            if (diff < 0.0) {
                ar_it++->on_stack_nd_idx = right_idx;
                ++cur_nd;
            } else {
                ar_it++->on_stack_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_ + 1);
                cur_nd = nd_arr_ + right_idx;
            }
            for (std::int32_t cur_dim = N - 1; cur_dim >= 0; --cur_dim) 
                cand_pts_matrix[cur_dim*kMaxBalancedTreeHeight + cand_pts_matrix_idx] = cur_nd->key[cur_dim];
            ++cand_pts_matrix_idx;
            cand_nd_ar_it++->cand_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_);
        } else {
            if (cur_nd++->dim_to_expand != N) {
                for (std::int32_t cur_dim = N - 1; cur_dim >= 0; --cur_dim)
                    cand_pts_matrix[cur_dim*kMaxBalancedTreeHeight + cand_pts_matrix_idx] = cur_nd->key[cur_dim];
                ++cand_pts_matrix_idx;
                cand_nd_ar_it++->cand_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_);
            }

            for (; cand_pts_matrix_idx >= 0; --cand_pts_matrix_idx) {
                FPType cand_pt_dist_sq = 0;
                for (std::int32_t cur_dim = N - 1; cur_dim >= 0; --cur_dim) {
                    auto& cand_pt_val_on_dim = cand_pts_matrix[cur_dim*kMaxBalancedTreeHeight + cand_pts_matrix_idx];
                    cand_pt_val_on_dim -= search_pt[cur_dim];
                    cand_pt_dist_sq += cand_pt_val_on_dim * cand_pt_val_on_dim;
                }
                if (cand_pt_dist_sq < best_dist_sq) {
                    best_dist_sq = cand_pt_dist_sq;
                    best_idx = ar_stack[cand_pts_matrix_idx].cand_nd_idx;
                }
            }
            cand_pts_matrix_idx = 0;
            cand_nd_ar_it = ar_stack.begin();
            do {
                if (ar_it == ar_stack.begin()) [[unlikely]]
                    return *(elem_arr_ + best_idx);
            } while ((--ar_it)->diff_on_dim_sq >= best_dist_sq);
            cur_nd = nd_arr_ + ar_it->on_stack_nd_idx;
            for (std::int32_t cur_dim = N - 1; cur_dim >= 0; --cur_dim)
                cand_pts_matrix[cur_dim*kMaxBalancedTreeHeight + cand_pts_matrix_idx] = cur_nd->key[cur_dim];
            ++cand_pts_matrix_idx;
            cand_nd_ar_it++->cand_nd_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_);
        }
    }
    return *(elem_arr_ + best_idx);
}
*/


template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <typename def::DistType thisDt,
typename std::enable_if<thisDt == def::DistType::kEuc, int>::type>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
NNValueHelper(MetaNode *cur, std::uint8_t dim, const PointND<FPType, N> &pt,
              const ElemType *&bestValue, FPType &bestDist) const {
    FPType curDist = PointND<FPType, N>::template
    dist<def::DistType::kEucSq>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    std::uint8_t next_dim = dim == N - 1 ? 0 : dim + 1;
    FPType diff = pt[dim] - cur->key[dim];
    MetaNode *next = diff < 0.0 ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    if (diff*diff < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <typename def::DistType thisDt,
typename std::enable_if<thisDt != def::DistType::kEuc, int>::type>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
NNValueHelper(MetaNode *cur, std::uint8_t dim, const PointND<FPType, N> &pt,
              const ElemType *&bestValue, FPType &bestDist) const {
    FPType curDist = PointND<FPType, N>::template dist<DT>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    std::uint8_t next_dim = dim == N - 1 ? 0 : dim + 1;
    MetaNode *next = pt[dim] < cur->key[dim] ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    if (branchMin(cur->key, pt, dim) < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
FPType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::branchMin(const PointND<FPType, N> &trPt,
                                                          const PointND<FPType, N> &searchPt, std::uint8_t dim) const {
    switch (DT) {
        case def::DistType::kEuc:
        case def::DistType::kMan:
            return std::fabs(trPt[dim] - searchPt[dim]);
            /*
             case DistType::kHav:
             PointND<FPType, N> pt;
             pt[idx] = searchPt[idx];
             pt[1-idx] = trPt[1-idx];
             return PointND<FPType, N>::havDist(trPt, pt);
             */
    }
}




#endif /* KDTreeExpandLongestVec_hpp */
