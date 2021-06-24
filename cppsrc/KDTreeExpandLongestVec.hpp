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
//#include "PoolAllocator.hpp"



#include "BoundedPQueue.hpp"

#include <llama/llama.hpp>
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
#include <array>

//#include <ranges>

namespace ns {
    struct NodePtr {};
    struct PtNodePtrArray {};
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
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
    
    // std::uint_fast8_t Dim() const;
    // Usage: std::uint_fast8_t dim = kd.Dim();
    // ----------------------------------------------------
    // Returns the Dim of the tree.
    constexpr std::uint_fast8_t dimension() const;
    typename PointND<FPType, N>::DistType distType() const;
    
    // std::uint_fast32_t size() const;
    // std::uint_fast32_t height() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree, the max height
    // and whether the tree is empty.
    std::uint_fast32_t size() const;
    std::uint_fast32_t cap() const;
    std::uint_fast32_t height() const;
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
        std::uint_fast8_t dim_to_expand;
        PointND<FPType, N> key;
        bool operator==(const MetaNode&) const = default;
        bool operator!=(const MetaNode&) const = default;
    };
    
    MetaNode *nd_arr_;
    ElemType *elem_arr_;
    uint32_t size_;
    uint32_t cap_;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    std::uint_fast32_t heightHelper(const MetaNode *n) const;
    
    void DeAlloc();
    
    // ----------------------------------------------------
    // Helper method for range constructor
    //template <std::random_access_iterator RAI>
    template <typename RAI>
    void RangeCtorHelper(RAI, RAI);

    template <class RAI>
    void RangeCtorRecursion(MetaNode*&, ElemType*&, RAI, RAI, std::array<FPType, N>&,
                                  std::array<FPType, N>&, std::array<FPType, N>&);
    
    template <class ConstRAI>
    static std::tuple<std::array<FPType, N>, std::array<FPType,N>, std::array<FPType, N>>
           ComputeInitBBoxHlSpread(ConstRAI, ConstRAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(MetaNode *cur, std::uint_fast8_t dim, const PointND<FPType, N> &pt,
                        BoundedPQueue<ElemType, FPType> &bpq) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const PointND<FPType, N>& key) const;
    
    template <typename PointND<FPType, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt == PointND<FPType, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(MetaNode*, std::uint_fast8_t, const PointND<FPType, N>&,
                       const ElemType *&, FPType&) const;
    
    template <typename PointND<FPType, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt != PointND<FPType, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(MetaNode*, std::uint_fast8_t, const PointND<FPType, N>&,
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
    
    FPType branchMin(const PointND<FPType, N>&, const PointND<FPType, N>&, std::uint_fast8_t) const;
    
};

/** KDTreeExpandLongestVec class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::KDTreeExpandLongestVec() {
    size_ = cap_ = 0;
    nd_arr_ = nullptr;
    elem_arr_ = nullptr;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs) : size_(rhs.size_), cap_(rhs.cap_) {
    nd_arr_ = new MetaNode[size_];
    std::copy_n(rhs.nd_arr_, size_, nd_arr_);
    elem_arr_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType), std::nothrow));
    std::uninitialized_copy_n(rhs.elem_arr_, size_, elem_arr_);
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>&
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator=(const KDTreeExpandLongestVec& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        DeAlloc();
        size_ = cap_ = rhs.size_;
        nd_arr_ = new MetaNode[size_];
        std::copy_n(rhs.nd_arr_, size_, nd_arr_);
        elem_arr_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType)));
        std::uninitialized_copy_n(rhs.elem_arr_, size_, elem_arr_);
    }
    return *this;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
KDTreeExpandLongestVec(KDTreeExpandLongestVec&& rhs) noexcept
: nd_arr_(rhs.nd_arr_), elem_arr_(rhs.elem_arr_), size_(rhs.size_), cap_(rhs.cap_) {
    rhs.nd_arr_ = nullptr;
    rhs.elem_arr_ = nullptr;
    rhs.size_ = rhs.cap_ = 0;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
operator=(KDTreeExpandLongestVec&& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        DeAlloc();
        nd_arr_ = rhs.nd_arr_;
        rhs.nd_arr_ = nullptr;
        elem_arr_ = rhs.elem_arr_;
        rhs.elem_arr_ = nullptr;
        size_ = rhs.size_;
        cap_ = rhs.cap_;
        rhs.size_ = rhs.cap_ = 0;
    }
    return *this;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::~KDTreeExpandLongestVec() {
    DeAlloc();
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::DeAlloc() {
    delete[] nd_arr_;
    std::destroy_n(elem_arr_, size_);
    ::operator delete(static_cast<void*>(elem_arr_));
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
//template <std::random_access_iterator ConstRAI> requires def::const_iterator<ConstRAI>
//&& std::same_as<typename std::iter_value_t<ConstRAI>, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
template <typename ConstRAI, std::enable_if_t<std::is_const_v<typename std::remove_pointer_t<typename std::iterator_traits<ConstRAI>::pointer>>, int>>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::KDTreeExpandLongestVec(ConstRAI cbegin, ConstRAI cend) : size_(static_cast<std::uint32_t>(cend - cbegin)), cap_(size_) {
    std::vector<node_type> constructData(cbegin, cend);
    RangeCtorHelper(constructData.begin(), constructData.end());
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
//template <std::random_access_iterator RAI> requires def::non_const_iterator<RAI>
//requires std::same_as<typename std::iterator_traits<RAI>::value_type, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
template <typename RAI, std::enable_if_t<!std::is_const_v<typename std::remove_pointer_t<typename std::iterator_traits<RAI>::pointer>>, int>>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::KDTreeExpandLongestVec(RAI begin, RAI end) :size_(static_cast<std::uint32_t>(end - begin)), cap_(size_) {
    RangeCtorHelper(begin, end);
}
/*
#include <tuple>
thread_local static constinit std::array<std::tuple<std::size_t, std::size_t, std::size_t, std::size_t>, 4> count_score_arr{ std::make_tuple(0, 0, 0, 0), std::make_tuple(0, 0, 0, 1), std::make_tuple(0, 0, 0, 2), std::make_tuple(0, 0, 0, 3) };
*/

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
//template <std::random_access_iterator RAI>
template <typename RAI>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::RangeCtorHelper(RAI begin, RAI end) {
    if (size_ == 0) [[unlikely]]
        return;
    nd_arr_ = new MetaNode[size_];
    elem_arr_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType), std::nothrow));
    if (size_ == 1) [[unlikely]] {
        *nd_arr_ = {0, N, std::move(begin->key)};
        new (elem_arr_) ElemType(std::move(begin->value));
        return;
    }
    auto [lows, highs, hl_spread] = ComputeInitBBoxHlSpread(std::as_const(begin), std::as_const(end));
    MetaNode* cur_nd = nd_arr_;
    ElemType* cur_elem = elem_arr_;
    RangeCtorRecursion(cur_nd, cur_elem, begin, end, lows, highs, hl_spread);
    /*
    //std::cout << twos << '\t' << threes << '\t' << fours << '\t' << others << std::endl;
    std::sort(count_score_arr.begin(), count_score_arr.end(), [](auto t1, auto t2) {return std::get<0>(t1) < std::get<0>(t2); });
    std::for_each(count_score_arr.begin(), count_score_arr.end(), [score = 0u](auto& t) mutable {std::get<1>(t) += std::get<0>(t); std::get<0>(t) = 0; std::get<2>(t) += score++; });
    std::sort(count_score_arr.begin(), count_score_arr.end(), [](auto t1, auto t2) {return std::get<3>(t1) < std::get<3>(t2); });
    std::cout << "total num of executions this branch in the order of 2s 3s 4s others : \n";
    std::for_each(count_score_arr.begin(), count_score_arr.end(), [](auto t) {std::cout  << std::get<1>(t) << '\t'; });
    std::cout << '\n';
    std::cout << "Frequency scores in the order of 2s 3s 4s others : \n";
    std::for_each(count_score_arr.begin(), count_score_arr.end(), [](auto t) {std::cout << std::get<2>(t) << '\t'; });
    std::cout << '\n';
    std::cout << '\n';
    std::cout << '\n';
    */





    /*
    struct ActRecord {
        std::uint32_t nd_idx_to_construct_right_child;
        std::uint_fast8_t dim;
        RAI this_begin, this_end;
    };
    
    constexpr std::uint_fast8_t kARStackSize = 32;
    //actRecord actSt[stackSize], *actStIt = actSt;
    std::tuple<std::uint32_t, std::uint32_t, RAI, RAI> actSt[kARStackSize], *actStIt = actSt;
    std::pair<FPType*, FPType> bboxChangeSt[kARStackSize], *curBboxChangeStIt = bboxChangeSt-1;
    RAI this_begin_it = begin, this_end_it = end;
    std::uint32_t depth = 0;
    
    
    while (true) {
        std::uint_fast8_t dim = static_cast<std::uint_fast8_t>(std::distance(hl_spread.cbegin(), std::max_element(hl_spread.cbegin(), hl_spread.cend())));
        RAI median = this_begin_it + (this_end_it-this_begin_it)/2;
        std::nth_element(this_begin_it, median, this_end_it, [dim](const auto& p1, const auto& p2) {return p1.key[dim] < p2.key[dim];});
        new (cur_nd++) MetaNode {0, dim, std::move(median->key)};
        new (cur_elem++) ElemType (std::move(median->value));
     
        switch (std::ptrdiff_t num_items_this_iteration = end - begin);
                num_items_this_iteration) {
        case 2:
            *cur_nd++ = { 0, N, std::move(begin->key) };
            new (cur_elem++) ElemType(std::move(begin->value));
            break;
        case 3:
                        *cur_nd++ = {0, N, std::move(begin->key)};
            *cur_nd = {0, N, std::move((median+1)->key)};
            nd_ptr_this_iter->right_idx = static_cast<std::uint32_t>(cur_nd++ - nd_arr_);
            new (cur_elem++) ElemType (std::move(begin->value));
            new (cur_elem++) ElemType (std::move((median+1)->value));
            break;
        case 4:
            FPType prev_high_on_dim = std::exchange(highs[dim], cur_nd->key[dim]);
            hl_spread[dim] = highs[dim] - lows[dim];
            new (actStIt++) std::tuple<std::uint32_t, std::uint32_t, RAI, RAI> {static_cast<std::uint32_t>(cur_nd-nd_arr_), dim, median + 1, this_end_it};
            this_end_it = median;
        } else {
            if (num_nds_this_iteration == 2) {
                new (cur_nd++) MetaNode {0, N, std::move(this_begin_it->key)};
                new (cur_elem++) ElemType (std::move(this_begin_it->value));
            }
        TRACEBACK:
            if (cur_nd - nd_arr_ == size_)
                return;
            const auto &[ndIdxToAssignRightIdx, stDepth, stThisBeginIt, stThisEndIt] = *--actStIt;
            (nd_arr_ + ndIdxToAssignRightIdx)->right_idx = static_cast<std::uint32_t>(cur_nd - nd_arr_);
            if ((this_begin_it = stThisBeginIt) + 1 == (this_end_it = stThisEndIt)) {
                new (cur_nd++) MetaNode {0, N, this_begin_it->key};
                new (cur_elem++) ElemType (this_begin_it->value);
                goto TRACEBACK;
            }
            
            
            while (depth != stDepth + 1 && --depth ) {
                *curBboxChangeStIt->first = curBboxChangeStIt->second;
                --curBboxChangeStIt;
            }
           
            FPType rightVal = *(curBboxChangeStIt->first-1);
            *(curBboxChangeStIt->first-1) = *curBboxChangeStIt->first;
            *curBboxChangeStIt->first = curBboxChangeStIt->second;
            *curBboxChangeStIt = {curBboxChangeStIt->first-1, rightVal};
        }
    }*/
    
     
     
    /*
    std::for_each_n(nd_arr_, size_, [](const MetaNode& nd){
        std::cout << nd.dim_to_expand << '\t' << nd.right_idx << '\t';
        std::copy(nd.key.begin(), nd.key.end(), std::ostream_iterator<FPType>(std::cout,","));
        std::cout << '\n';
    });
    */
}



template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <class ConstRAI>
std::tuple<std::array<FPType, N>, std::array<FPType, N>, std::array<FPType, N>> KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
ComputeInitBBoxHlSpread(ConstRAI cbegin, ConstRAI cend) {
    std::array<FPType, N> lows, highs, hl_spread;
    lows.fill(std::numeric_limits<FPType>::max());
    highs.fill(std::numeric_limits<FPType>::min());
    std::for_each(cbegin, cend, [&lows, &highs](const auto &nh) mutable {
        for (std::uint_fast8_t i = 0; i < N; ++i) {
            FPType pt_val_on_dim = nh.key[i];
            auto &bbox_low = lows[i], &bbox_high = highs[i];
            bbox_low = std::min(bbox_low, pt_val_on_dim);
            bbox_high = std::max(bbox_high, pt_val_on_dim);
        }
    });
    std::transform(highs.cbegin(), highs.cend(), lows.cbegin(), hl_spread.begin(), std::minus<FPType>());
    return {lows, highs, hl_spread};
}



template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <class RAI>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
RangeCtorRecursion(MetaNode *&cur_nd, ElemType *&cur_elem, RAI begin, RAI end,
                   std::array<FPType, N> &lows, std::array<FPType, N> &highs, std::array<FPType, N> &hl_spread) {
    std::uint_fast8_t dim = static_cast<std::uint_fast8_t>(std::distance(hl_spread.cbegin(), std::max_element(hl_spread.cbegin(), hl_spread.cend())));
    RAI median = begin + (end - begin)/2;
    std::nth_element(begin, median, end, [dim](const auto& p1, const auto& p2) {return p1.key[dim] < p2.key[dim];});
    //oneapi::dpl::nth_element(oneapi::dpl::execution::par_unseq, begin, median, end, [dim](const auto& p1, const auto& p2) {return p1.key[dim] < p2.key[dim]; });

    new (cur_elem++) ElemType (std::move(median->value));
    

    switch (std::ptrdiff_t num_items_this_iteration = end - begin;
    num_items_this_iteration) {
    [[likely]] case 5 ... std::numeric_limits<std::uint32_t>::max(): {// ++std::get<0>(count_score_arr[3]);
        *cur_nd = { static_cast<std::uint32_t>((cur_nd - nd_arr_) + (median - begin)) + 1, dim, std::move(median->key) };
        FPType prev_high_low_on_dim = std::exchange(highs[dim], cur_nd++->key[dim]);
        hl_spread[dim] = highs[dim] - lows[dim];
        RangeCtorRecursion(cur_nd, cur_elem, begin, median, lows, highs, hl_spread);
        hl_spread[dim] = prev_high_low_on_dim - highs[dim];
        std::swap(lows[dim], highs[dim]);
        std::swap(highs[dim], prev_high_low_on_dim);
        RangeCtorRecursion(cur_nd, cur_elem, median + 1, end, lows, highs, hl_spread);
        lows[dim] = prev_high_low_on_dim;
        hl_spread[dim] = highs[dim] - prev_high_low_on_dim;
        break;
    }
    case 2: {//++std::get<0>(count_score_arr[0]);
        *cur_nd++ = { 0, dim, std::move(median->key) };
        *cur_nd++ = { 0, N, std::move(begin->key) };
        new (cur_elem++) ElemType(std::move(begin->value));
        break;
    }
    case 4: {//++std::get<0>(count_score_arr[2]);
        *cur_nd = { static_cast<std::uint32_t>(cur_nd - nd_arr_) + 3, dim, std::move(median->key) };
        hl_spread[dim] = cur_nd++->key[dim] - lows[dim];
        dim = static_cast<std::uint_fast8_t>(std::distance(hl_spread.cbegin(), std::max_element(hl_spread.cbegin(), hl_spread.cend())));
        RAI child_or_parent = (begin->key[dim] < (begin + 1)->key[dim]) ? begin++ : begin + 1;
        *cur_nd++ = { 0, dim, std::move(begin->key) };
        new (cur_elem++) ElemType(std::move(begin->value));
        *cur_nd++ = { 0, N, std::move(child_or_parent->key) };
        new (cur_elem++) ElemType(std::move(child_or_parent->value));
        *cur_nd++ = { 0, N, std::move((median + 1)->key) };
        new (cur_elem++) ElemType(std::move((median + 1)->value));
        hl_spread[dim] = highs[dim] - lows[dim];
        break;
    }
    case 3: {// ++std::get<0>(count_score_arr[1]);
        *cur_nd++ = { static_cast<std::uint32_t>(cur_nd - nd_arr_) + 2, dim, std::move(median->key) };
        *cur_nd++ = { 0, N, std::move(begin->key) };
        *cur_nd++ = { 0, N, std::move((median + 1)->key) };
        new (cur_elem++) ElemType(std::move(begin->value));
        new (cur_elem++) ElemType(std::move((median + 1)->value));
        break;
    }
    }

    

    /*
    std::uint_fast8_t dim = static_cast<std::uint_fast8_t>(std::distance(hl_spread.cbegin(), std::max_element(hl_spread.cbegin(), hl_spread.cend())));
    RAI median = begin + (end - begin) / 2;
    std::nth_element(begin, median, end, [dim](const auto& p1, const auto& p2) {return p1.key[dim] < p2.key[dim]; });
    //oneapi::dpl::nth_element(oneapi::dpl::execution::par_unseq, begin, median, end, [dim](const auto& p1, const auto& p2) {return p1.key[dim] < p2.key[dim]; });

    *cur_nd = { median + 1 != end ? static_cast<std::uint32_t>((cur_nd - nd_arr_) + (median - begin)) + 1 : 0, dim, std::move(median->key) };
    new (cur_elem++) ElemType(std::move(median->value));
    FPType this_pt_val_on_dim = cur_nd->key[dim];

    if (begin == median - 1) {
        *++cur_nd = {0, N, std::move(begin->key)};
        ++cur_nd;
        new (cur_elem++) ElemType (std::move(begin->value));
    } else {
        FPType prev_high_on_dim = std::exchange(highs[dim], this_pt_val_on_dim);
        hl_spread[dim] = highs[dim] - lows[dim];
        RangeCtorRecursion(++cur_nd, cur_elem, begin, median, lows, highs, hl_spread);
        highs[dim] = prev_high_on_dim;
        hl_spread[dim] = prev_high_on_dim - lows[dim];
    }
    
    if (median + 1 != end) {
        if (median + 2 == end) {
            *cur_nd++ = {0, N, std::move((median+1)->key)};
            new (cur_elem++) ElemType (std::move((median+1)->value));
        } else {
            FPType prev_low_on_dim = std::exchange(lows[dim], this_pt_val_on_dim);
            hl_spread[dim] = highs[dim] - lows[dim];
            RangeCtorRecursion(cur_nd, cur_elem, median+1, end, lows, highs, hl_spread);
            lows[dim] = prev_low_on_dim;
            hl_spread[dim] = highs[dim] - prev_low_on_dim;
        }
    } */
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
constexpr std::uint_fast8_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::dimension() const {
    return N;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename PointND<FPType, N>::DistType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::distType() const {
    return DT;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
std::uint_fast32_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::size() const {
    return size_;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
std::uint_fast32_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::cap() const {
    return cap_;
}


template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
std::uint_fast32_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::height() const {
    const MetaNode *root = &nd_arr_[0];
    return heightHelper(root);
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
std::uint_fast32_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::heightHelper(const MetaNode *n) const {
    //return n ? 1 + std::max(heightHelper(n->left), heightHelper(n->right)) : -1;
    return 1;
}


template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::Empty() const {
    return size_ == 0;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator==(const KDTreeExpandLongestVec& rhs) const {
    if (size_ != rhs.size_)
        return false;
    return std::equal(elem_arr_, elem_arr_ + size_, rhs.elem_arr_) && std::equal(nd_arr_, nd_arr_ + size_, rhs.nd_arr_);
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator!=(const KDTreeExpandLongestVec& rhs) const {
    return !(*this == rhs);
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::PrintTreeInfo() const {
    std::cout << "Tree height is " << height()
    << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::Clear() {
    std::destroy_n(nd_arr_, size_);
    std::destroy_n(elem_arr_, size_);
    size_ = 0;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::ShrinkToFit() {
    if (cap_ > size_) {
        // TO DO
    }
    std::destroy_n(nd_arr_, size_);
    std::destroy_n(elem_arr_, size_);
    size_ = 0;
}


template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
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

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::contains(const PointND<FPType, N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator[] (const PointND<FPType, N>& pt) {
    MetaNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new MetaNode(pt, ElemType());
        //treeSize++;
    }
    return (*ndPtr)->object;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeExpandLongestVec&>(*this).at(pt));
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
const ElemType& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) const {
    MetaNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::MetaNode**
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) {
    return const_cast<MetaNode**>(static_cast<const KDTreeExpandLongestVec*>(this)->findNodePtr(pt));
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::MetaNode*const*
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) const {
    //MetaNode *const*n = &root;
    MetaNode *const*n;
    for (std::uint_fast8_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1){}
    //n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
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

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::kNNValueHelper(MetaNode *cur, std::uint_fast8_t dim,
                                                             const PointND<FPType, N>& pt, BoundedPQueue<ElemType, FPType> &bpq) const {
    bpq.enqueue(cur->object, PointND<FPType, N>::template dist<DT>(cur->key, pt));
    std::uint_fast8_t next_dim = dim + 1 < N ? dim + 1 : 0;
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

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename NdTypeOutIt> 
    //requires std::output_iterator<NdTypeOutIt, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
NdTypeOutIt KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
NNsWithFence(const PointND<FPType, N>& pt, FPType fence_sq, NdTypeOutIt pe_out_it) const {
    constexpr std::uint_fast8_t kMaxBalancedTreeHeight = 32;
    struct ActRecord {
        FPType diff_on_dim_sq;
        const MetaNode* nd;
    } act_record_stack[kMaxBalancedTreeHeight], *act_record_it = act_record_stack;
    const MetaNode *cur_nd = nd_arr_;
    FPType cur_dist_sq, best_dist_sq = PointND<FPType, N>::template dist<PointND<FPType, N>::DistType::EUCSQ>(cur_nd->key, pt);
    FPType best_dist_plus_fence_sq = best_dist_sq + fence_sq + FPType(2.0)*sqrt(fence_sq*best_dist_sq);
    
    constexpr std::uint_fast16_t kMaxDpeSizeOnStack = 1152;
    struct DistPtElem {
        FPType dist_sq;
        const PointND<FPType, N>* pt;
        const ElemType* elem;
    } result_dpe_arr[kMaxDpeSizeOnStack], *result_dpe_arr_it = result_dpe_arr,
      *result_dist_pe_arr_end = result_dpe_arr + kMaxDpeSizeOnStack;
    *result_dpe_arr_it++ = {best_dist_sq, &nd_arr_->key, elem_arr_};
    std::vector<DistPtElem> result_dpe_vec;
   
    while (true) {
        if (std::uint32_t right_idx = cur_nd->right_idx) {
            std::uint_fast8_t dim = cur_nd->dim_to_expand;
            FPType diff = pt[dim] - cur_nd->key[dim];
            if (diff < 0.0) {
                new(act_record_it++) ActRecord {diff*diff, nd_arr_ + right_idx};
                ++cur_nd;
            } else {
                new(act_record_it++) ActRecord {diff*diff, cur_nd+1};
                cur_nd = nd_arr_ + right_idx;
            }
        } else if (cur_nd++->dim_to_expand == N) {
            do {
                if (act_record_it == act_record_stack) [[unlikely]] {
                    auto filter_transform_dpe_arr_to_pe_arr = [&pe_out_it, best_dist_plus_fence_sq](const auto &dpe) mutable {
                        if (dpe.dist_sq < best_dist_plus_fence_sq)
                            *pe_out_it++ = { *(dpe.pt), *(dpe.elem) };
                    };
                    std::for_each(std::make_reverse_iterator(result_dpe_arr_it), std::make_reverse_iterator(result_dpe_arr), filter_transform_dpe_arr_to_pe_arr);
                    std::for_each(result_dpe_vec.rbegin(), result_dpe_vec.rend(), filter_transform_dpe_arr_to_pe_arr); 
                    return pe_out_it;
                }
            } while ((--act_record_it)->diff_on_dim_sq > best_dist_plus_fence_sq);
            cur_nd = act_record_it->nd;
        }
        
        cur_dist_sq = PointND<FPType, N>::template dist<PointND<FPType, N>::DistType::EUCSQ>(cur_nd->key, pt);
        if (cur_dist_sq < best_dist_plus_fence_sq) {
            if (cur_dist_sq < best_dist_sq) {
                best_dist_sq = cur_dist_sq;
                best_dist_plus_fence_sq = best_dist_sq + fence_sq + FPType(2.0)*sqrt(fence_sq*best_dist_sq);
            }
            *result_dpe_arr_it++ = {cur_dist_sq, &cur_nd->key, &elem_arr_[cur_nd-nd_arr_]};
            if (result_dpe_arr_it == result_dist_pe_arr_end) [[unlikely]] {
                result_dpe_arr_it = result_dpe_arr;
                result_dpe_vec.insert(result_dpe_vec.end(), result_dpe_arr, result_dist_pe_arr_end);
            }
        }
    }
    return pe_out_it;
}

/*
constexpr std::uint_fast8_t kMaxBalancedTreeHeight = 32;
using PtNodePtr = llama::Record <
        llama::Field<ns::Point, llama::Array<FPType, N>>,
        llama::Field<ns::NodePtr, MetaNode*>
    >;
llama::ArrayDims<N> pt_ndPtr_arr_size;
std::fill(pt_ndPtr_arr_size.begin(), pt_ndPtr_arr_size.end(), kMaxBalancedTreeHeight);
auto pt_ndPtr_soa_view = llama::allocView(llama::mapping::SoA<llama::ArrayDims<N>, PtNodePtr>(pt_ndPtr_arr_size));
pt_ndPtr_soa_view(0ull, 0ull, 0ull)(ns::NodePtr{}) = nd_arr_;
*/

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N> &search_pt) const {

    constexpr std::uint_fast8_t kMaxBalancedTreeHeight = 32;
    struct ActRecord {
        FPType diff_on_dim_sq;
        const MetaNode* nd;
    } act_record_stack[kMaxBalancedTreeHeight], *act_record_it = act_record_stack;
    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    const MetaNode *cur = nd_arr_, *best_node = nd_arr_;
    FPType cur_dist_sq, best_dist_sq = PointND<FPType, N>::template dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, search_pt);
    
    while (true) {
        if (std::uint32_t right_idx = cur->right_idx) {
            std::uint_fast8_t dim = cur->dim_to_expand;
            FPType diff = search_pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(act_record_it++) ActRecord {diff*diff, nd_arr_ + right_idx};
                ++cur;
            } else {
                new(act_record_it++) ActRecord {diff*diff, cur+1};
                cur = nd_arr_ + right_idx;
            }
        } else if (cur++->dim_to_expand == N) {



            do {
                if (act_record_it == act_record_stack)
                    return elem_arr_[best_node - nd_arr_];
            } while ((--act_record_it)->diff_on_dim_sq >= best_dist_sq);
            cur = act_record_it->nd;
        }
        cur_dist_sq = PointND<FPType, N>::template dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, search_pt);
        if (cur_dist_sq < best_dist_sq) {
            best_dist_sq = cur_dist_sq;
            best_node = cur;
        }
    }
    return elem_arr_[best_node - nd_arr_];
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename PointND<FPType, N>::DistType thisDt,
typename std::enable_if<thisDt == PointND<FPType, N>::DistType::EUC, int>::type>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
NNValueHelper(MetaNode *cur, std::uint_fast8_t dim, const PointND<FPType, N> &pt,
              const ElemType *&bestValue, FPType &bestDist) const {
    FPType curDist = PointND<FPType, N>::template
    dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    std::uint_fast8_t next_dim = dim == N - 1 ? 0 : dim + 1;
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

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename PointND<FPType, N>::DistType thisDt,
typename std::enable_if<thisDt != PointND<FPType, N>::DistType::EUC, int>::type>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
NNValueHelper(MetaNode *cur, std::uint_fast8_t dim, const PointND<FPType, N> &pt,
              const ElemType *&bestValue, FPType &bestDist) const {
    FPType curDist = PointND<FPType, N>::template dist<DT>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    std::uint_fast8_t next_dim = dim == N - 1 ? 0 : dim + 1;
    MetaNode *next = pt[dim] < cur->key[dim] ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    if (branchMin(cur->key, pt, dim) < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    }
}

template <typename FPType, std::uint_fast8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
FPType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::branchMin(const PointND<FPType, N> &trPt,
                                                          const PointND<FPType, N> &searchPt, std::uint_fast8_t dim) const {
    switch (DT) {
        case PointND<FPType, N>::DistType::EUC:
        case PointND<FPType, N>::DistType::MAN:
            return std::fabs(trPt[dim] - searchPt[dim]);
            /*
             case DistType::HAV:
             PointND<FPType, N> pt;
             pt[idx] = searchPt[idx];
             pt[1-idx] = trPt[1-idx];
             return PointND<FPType, N>::havDist(trPt, pt);
             */
    }
}




#endif /* KDTreeExpandLongestVec_hpp */
