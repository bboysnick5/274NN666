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
#include <ranges>

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
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
    template <std::random_access_iterator RAI> requires def::non_const_iterator<RAI>
    //requires std::same_as<typename std::iterator_traits<RAI>::value_type, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
    KDTreeExpandLongestVec(RAI, RAI);
    
    template <std::random_access_iterator ConstRAI> requires def::const_iterator<ConstRAI>
    //&& std::same_as<typename std::iter_value_t<ConstRAI>, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
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
    
    // std::size_t Dimension() const;
    // Usage: std::size_t dim = kd.Dimension();
    // ----------------------------------------------------
    // Returns the Dimension of the tree.
    constexpr std::size_t dimension() const;
    typename PointND<FPType, N>::DistType distType() const;
    
    // std::size_t size() const;
    // std::size_t height() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree, the max height
    // and whether the tree is empty.
    std::size_t size() const;
    std::size_t cap() const;
    int height() const;
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
    
    
    struct TreeNode {
        unsigned int right_idx;
        unsigned int dim_to_expand;
        PointND<FPType, N> key;
        bool operator==(const TreeNode&) const = default;
        bool operator!=(const TreeNode&) const = default;
    };
    
    TreeNode *nd_vec_;
    ElemType *elem_vec_;
    uint32_t size_;
    uint32_t cap_;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    int heightHelper(const TreeNode *n) const;
    
    void DeAlloc();
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <std::random_access_iterator RAI>
    void RangeCtorHelper(RAI, RAI);

    template <class RAI>
    void RangeCtorRecursionHelper(TreeNode*&, ElemType*&, RAI, RAI, std::array<FPType, N>&,
                                  std::array<FPType, N>&, std::array<FPType, N>&);
    
    template <class ConstRAI>
    static std::tuple<std::array<FPType, N>, std::array<FPType,N>, std::array<FPType, N>>
           ComputeInitBBoxSpans(ConstRAI, ConstRAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(TreeNode *cur, std::size_t dim, const PointND<FPType, N> &pt,
                        BoundedPQueue<ElemType, FPType> &bpq) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const PointND<FPType, N>& key) const;
    
    template <typename PointND<FPType, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt == PointND<FPType, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, std::size_t, const PointND<FPType, N>&,
                       const ElemType *&, FPType&) const;
    
    template <typename PointND<FPType, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt != PointND<FPType, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, std::size_t, const PointND<FPType, N>&,
                       const ElemType*&, FPType&) const;
    
    
    // TreeNode** findNodePtr(const PointND<FPType, N>& pt);
    // TreeNode*const* findNodePtr(const PointND<FPType, N>& pt) const;
    // Usage: TreeNode **nodePtr = findNodePtr(pt);
    // ----------------------------------------------------
    // Returns the pointer pointing to the node address
    // corresponding to the given PointND. In this FPType pointing
    // fashion, we can construct a node at that location.
    TreeNode** findNodePtr(const PointND<FPType, N>& pt);
    TreeNode*const* findNodePtr(const PointND<FPType, N>& pt) const;
    
    FPType branchMin(const PointND<FPType, N>&, const PointND<FPType, N>&, std::size_t) const;
    
};

/** KDTreeExpandLongestVec class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::KDTreeExpandLongestVec() {
    size_ = cap_ = 0;
    nd_vec_ = nullptr;
    elem_vec_ = nullptr;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs) : size_(rhs.size_), cap_(rhs.cap_) {
    nd_vec_ = new TreeNode[size_];
    std::copy_n(rhs.nd_vec_, size_, nd_vec_);
    elem_vec_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType), std::nothrow));
    std::uninitialized_copy_n(rhs.elem_vec_, size_, elem_vec_);
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>&
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator=(const KDTreeExpandLongestVec& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        DeAlloc();
        size_ = cap_ = rhs.size_;
        nd_vec_ = new TreeNode[size_];
        std::copy_n(rhs.nd_vec_, size_, nd_vec_);
        elem_vec_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType)));
        std::uninitialized_copy_n(rhs.elem_vec_, size_, elem_vec_);
    }
    return *this;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
KDTreeExpandLongestVec(KDTreeExpandLongestVec&& rhs) noexcept
: nd_vec_(rhs.nd_vec_), elem_vec_(rhs.elem_vec_), size_(rhs.size_), cap_(rhs.cap_) {
    rhs.nd_vec_ = nullptr;
    rhs.elem_vec_ = nullptr;
    rhs.size_ = rhs.cap_ = 0;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
operator=(KDTreeExpandLongestVec&& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        DeAlloc();
        nd_vec_ = rhs.nd_vec_;
        rhs.nd_vec_ = nullptr;
        elem_vec_ = rhs.elem_vec_;
        rhs.elem_vec_ = nullptr;
        size_ = rhs.size_;
        cap_ = rhs.cap_;
        rhs.size_ = rhs.cap_ = 0;
    }
    return *this;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::~KDTreeExpandLongestVec() {
    DeAlloc();
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::DeAlloc() {
    delete[] nd_vec_;
    std::destroy_n(elem_vec_, size_);
    ::operator delete(static_cast<void*>(elem_vec_));
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <std::random_access_iterator ConstRAI> requires def::const_iterator<ConstRAI>
//&& std::same_as<typename std::iter_value_t<ConstRAI>, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::KDTreeExpandLongestVec(ConstRAI cbegin, ConstRAI cend) : size_(static_cast<uint32_t>(cend - cbegin)), cap_(size_) {
    std::vector<node_type> constructData(cbegin, cend);
    RangeCtorHelper(constructData.begin(), constructData.end());
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <std::random_access_iterator RAI> requires def::non_const_iterator<RAI>
//requires std::same_as<typename std::iterator_traits<RAI>::value_type, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::KDTreeExpandLongestVec(RAI begin, RAI end) :size_(static_cast<uint32_t>(end - begin)), cap_(size_) {
    RangeCtorHelper(begin, end);
    
    
    /*
    struct actRecord {
        unsigned int ndIdxToAssignRightIdx;
        unsigned int depth;
        RAI thisBeginIt, thisEndIt;
    }; */
    
    /*
    
    std::size_t stackSize = static_cast<std::size_t>(log2(size_+1));
    //actRecord actSt[stackSize], *actStIt = actSt;
    std::tuple<unsigned int, unsigned int, RAI, RAI> actSt[stackSize], *actStIt = actSt;
    std::pair<FPType*, FPType> bboxChangeSt[stackSize], *curBboxChangeStIt = bboxChangeSt-1;
    RAI thisBeginIt = begin, thisEndIt = end;
    unsigned int depth = 0;
    
    
    while (true) {
        auto dim = static_cast<unsigned int>(std::distance(
                           spans.cbegin(),std::max_element(spans.cbegin(), spans.cend())));
        
        TreeNode* curNdPtr = curNd;
        RAI median = thisBeginIt + (thisEndIt-thisBeginIt)/2;
        std::nth_element(thisBeginIt, median, thisEndIt,
                         [=](const auto& p1, const auto& p2) {
                             return p1.key[dim] < p2.key[dim];});
        new (curNd++) TreeNode {0, dim, median->key};
        new (curObj++) ElemType (median->value);
     
        auto numElems = thisEndIt - thisBeginIt;
        if (numElems > 2) {
            *++curBboxChangeStIt = {highs.data() + dim, highs[dim]};
            highs[dim] = curNdPtr->key[dim];
            //new (actStIt++) actRecord {static_cast<unsigned int>(curNdPtr-nd_vec_), depth++, median + 1, thisEndIt};
            new (actStIt++) std::tuple<unsigned int, unsigned int, RAI, RAI> {static_cast<unsigned int>(curNdPtr-nd_vec_), depth++, median + 1, thisEndIt};
            thisEndIt = median;
        } else {
            if (numElems == 2) {
                new (curNd++) TreeNode {0, N, thisBeginIt->key};
                new (curObj++) ElemType (thisBeginIt->value);
            }
        TRACEBACK:
            if (curNd - nd_vec_ == size_)
                return;
            const auto &[ndIdxToAssignRightIdx, stDepth, stThisBeginIt, stThisEndIt] = *--actStIt;
            (nd_vec_ + ndIdxToAssignRightIdx)->right_idx = static_cast<unsigned int>(curNd - nd_vec_);
            if ((thisBeginIt = stThisBeginIt) + 1 == (thisEndIt = stThisEndIt)) {
                new (curNd++) TreeNode {0, N, thisBeginIt->key};
                new (curObj++) ElemType (thisBeginIt->value);
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
    }
    
     */
     
    /*
    std::for_each_n(nd_vec_, size_, [](const TreeNode& nd){
        std::cout << nd.dim_to_expand << '\t' << nd.right_idx << '\t';
        std::copy(nd.key.begin(), nd.key.end(), std::ostream_iterator<FPType>(std::cout,","));
        std::cout << '\n';
    });
    */
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <std::random_access_iterator RAI> 
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::RangeCtorHelper(RAI begin, RAI end) {
    if (size_ == 0) [[unlikely]]
        return;
    nd_vec_ = new TreeNode[size_];
    elem_vec_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType), std::nothrow));
    if (size_ == 1) [[unlikely]] {
        *nd_vec_ = {0, N, std::move(begin->key)};
        new (elem_vec_) ElemType(std::move(begin->value));
        return;
    }
    auto [lows, highs, spans] = ComputeInitBBoxSpans(std::as_const(begin), std::as_const(end));
    auto curNd = nd_vec_;
    auto curObj = elem_vec_;
    RangeCtorRecursionHelper(curNd, curObj, begin, end, lows, highs, spans);
}


template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <class ConstRAI>
std::tuple<std::array<FPType, N>, std::array<FPType, N>, std::array<FPType, N>> KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
ComputeInitBBoxSpans(ConstRAI cbegin, ConstRAI cend) {
    std::array<FPType, N> lows, highs, spans;
    lows.fill(std::numeric_limits<FPType>::max());
    highs.fill(std::numeric_limits<FPType>::min());
    std::for_each(cbegin, cend, [&](const auto &nh) mutable {
        for (std::size_t i = 0; i < N; ++i) {
            FPType ptValOnithDim = nh.key[i];
            auto &bboxLow = lows[i], &bboxHigh = highs[i];
            bboxLow = std::min(bboxLow, ptValOnithDim);
            bboxHigh = std::max(bboxHigh, ptValOnithDim);
        }
    });
    std::transform(highs.crbegin(), highs.crend(), lows.crbegin(), spans.rbegin(), std::minus<FPType>());
    return {lows, highs, spans};
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <class RAI>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
RangeCtorRecursionHelper(TreeNode *&curNd, ElemType *&curObj, RAI begin, RAI end,
                         std::array<FPType, N> &lows, std::array<FPType, N> &highs, std::array<FPType, N> &spans) {
    auto dim = static_cast<unsigned int>(std::distance(spans.cbegin(), std::max_element(spans.cbegin(), spans.cend())));
    TreeNode* ndPtrThisIter = curNd;
    RAI median = begin + (end - begin)/2;
    std::nth_element(begin, median, end, [=](const auto& p1, const auto& p2) {return p1.key[dim] < p2.key[dim];});
    //oneapi::dpl::nth_element(oneapi::dpl::execution::par_unseq, begin, median, end, [=](const auto& p1, const auto& p2) {return p1.key[dim] < p2.key[dim]; });

    *curNd++ = {0, dim, std::move(median->key)};
    new (curObj++) ElemType (std::move(median->value));
        
    if (begin == median - 1) {
        *curNd++ = {0, N, std::move(begin->key)};
        new (curObj++) ElemType (std::move(begin->value));
    } else {
        FPType curValOnDim = ndPtrThisIter->key[dim];
        auto prevDimHigh = highs[dim];
        highs[dim] = curValOnDim;
        spans[dim] = curValOnDim - lows[dim];
        RangeCtorRecursionHelper(curNd, curObj, begin, median, lows, highs, spans);
        highs[dim] = prevDimHigh;
        spans[dim] = prevDimHigh - lows[dim];
    }
    
    if (median + 1 != end) {
        if (median + 2 == end) {
            ndPtrThisIter->right_idx = static_cast<unsigned int>(curNd - nd_vec_);
            *curNd++ = {0, N, std::move((median+1)->key)};
            new (curObj++) ElemType (std::move((median+1)->value));
        } else {
            ndPtrThisIter->right_idx = static_cast<unsigned int>(curNd - nd_vec_);
            FPType curValOnDim = ndPtrThisIter->key[dim];
            auto prevDimLow = lows[dim];
            lows[dim] = curValOnDim;
            spans[dim] = highs[dim] - curValOnDim;
            RangeCtorRecursionHelper(curNd, curObj, median+1, end, lows, highs, spans);
            lows[dim] = prevDimLow;
            spans[dim] = highs[dim] - prevDimLow;
        }
    }
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
constexpr std::size_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::dimension() const {
    return N;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename PointND<FPType, N>::DistType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::distType() const {
    return DT;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
std::size_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::size() const {
    return size_;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
std::size_t KDTreeExpandLongestVec<FPType, N, ElemType, DT>::cap() const {
    return cap_;
}


template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
int KDTreeExpandLongestVec<FPType, N, ElemType, DT>::height() const {
    const TreeNode *root = &nd_vec_[0];
    return heightHelper(root);
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
int KDTreeExpandLongestVec<FPType, N, ElemType, DT>::heightHelper(const TreeNode *n) const {
    //return n ? 1 + std::max(heightHelper(n->left), heightHelper(n->right)) : -1;
    return 1;
}


template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::Empty() const {
    return size_ == 0;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator==(const KDTreeExpandLongestVec& rhs) const {
    if (size_ != rhs.size_)
        return false;
    return std::equal(elem_vec_, elem_vec_ + size_, rhs.elem_vec_) && std::equal(nd_vec_, nd_vec_ + size_, rhs.nd_vec_);
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator!=(const KDTreeExpandLongestVec& rhs) const {
    return !(*this == rhs);
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::PrintTreeInfo() const {
    std::cout << "Tree height is " << height()
    << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::Clear() {
    std::destroy_n(nd_vec_, size_);
    std::destroy_n(elem_vec_, size_);
    size_ = 0;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::ShrinkToFit() {
    if (cap_ > size_) {
        // TO DO
    }
    std::destroy_n(nd_vec_, size_);
    std::destroy_n(elem_vec_, size_);
    size_ = 0;
}


template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
insert(const PointND<FPType, N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        //(*ndPtr)->object = value;
    } else {
        //*ndPtr = new TreeNode(0, pt, value);
        //treeSize++;
    }
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeExpandLongestVec<FPType, N, ElemType, DT>::contains(const PointND<FPType, N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::operator[] (const PointND<FPType, N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        //treeSize++;
    }
    return (*ndPtr)->object;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeExpandLongestVec&>(*this).at(pt));
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
const ElemType& KDTreeExpandLongestVec<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::TreeNode**
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTreeExpandLongestVec*>(this)->findNodePtr(pt));
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::TreeNode*const*
KDTreeExpandLongestVec<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) const {
    //TreeNode *const*n = &root;
    TreeNode *const*n;
    for (std::size_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1){}
    //n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
kNNValue(const PointND<FPType, N>& pt, std::size_t k) const {
    //if (empty())
    //    return ElemType();
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType, FPType> bpq(k);
    //kNNValueHelper(ndVec[0], 0, pt, bpq);
    
    std::multimap<std::size_t, ElemType, std::greater<std::size_t>> freqMap;
    while (!bpq.empty()) {
        ElemType elem = bpq.dequeueMin();
        for (typename std::multimap<std::size_t, ElemType>::iterator
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

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::kNNValueHelper(TreeNode *cur, std::size_t dim,
                                                             const PointND<FPType, N>& pt, BoundedPQueue<ElemType, FPType> &bpq) const {
    bpq.enqueue(cur->object, PointND<FPType, N>::template dist<DT>(cur->key, pt));
    std::size_t nextDim = dim + 1 < N ? dim + 1 : 0;
    TreeNode *next = pt[dim] < cur->key[dim] ? cur->left : cur->right;
    if (next)
        kNNValueHelper(next, nextDim, pt, bpq);
    if (bpq.size() < bpq.maxSize()
        || branchMin(cur->key, pt, dim) < bpq.worst()) {
        TreeNode *other = next == cur->left ? cur->right : cur->left;
        if (other)
            kNNValueHelper(other, nextDim, pt, bpq);
    }
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename NdTypeOutIt> 
    //requires std::output_iterator<NdTypeOutIt, typename KDTreeExpandLongestVec<FPType, N, ElemType, DT>::node_type>
NdTypeOutIt KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
NNsWithFence(const PointND<FPType, N>& pt, FPType fence_sq, NdTypeOutIt pe_out_it) const {
    
    struct ActRecord {
        FPType dist;
        const TreeNode* nd;
    } act_record_stack[32], *act_record_it = act_record_stack;
    FPType best_dist_sq = std::numeric_limits<FPType>::max(),
        best_dist_plus_fence_sq = best_dist_sq;
    const TreeNode *cur = nd_vec_;
    
    const std::size_t MAX_DISTPTELEMS_ON_STACK = 1152;
    struct DistPtElem {
        FPType dist;
        const PointND<FPType, N>* pt;
        const ElemType* elem;
    } result_dpe_arr[MAX_DISTPTELEMS_ON_STACK], *result_dpe_arr_it = result_dpe_arr,
      *result_dist_pe_arr_end = result_dpe_arr + MAX_DISTPTELEMS_ON_STACK;
    std::vector<DistPtElem> result_dpe_vec;
   
    while (true) {
        FPType curDistSq = PointND<FPType, N>::template dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDistSq < best_dist_plus_fence_sq) {
            if (curDistSq < best_dist_sq) {
                best_dist_sq = curDistSq;
                best_dist_plus_fence_sq = best_dist_sq + fence_sq + FPType(2.0)*sqrt(fence_sq*best_dist_sq);
            }
            *result_dpe_arr_it++ = {curDistSq, &cur->key, &elem_vec_[cur-nd_vec_]};
            if (result_dpe_arr_it == result_dist_pe_arr_end) [[unlikely]] {
                result_dpe_arr_it = result_dpe_arr;
                result_dpe_vec.insert(result_dpe_vec.end(), result_dpe_arr, result_dist_pe_arr_end);
            }
        }
        
        if (unsigned int right_idx = cur->right_idx) {
            unsigned int dim = cur->dim_to_expand;
            FPType diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(act_record_it++) ActRecord {diff*diff, nd_vec_ + right_idx};
                ++cur;
            } else {
                new(act_record_it++) ActRecord {diff*diff, cur+1};
                cur = nd_vec_ + right_idx;
            }
        } else if (cur++->dim_to_expand == N) {
            do {
                if (act_record_it == act_record_stack) [[unlikely]] {
                    auto within_dist_plus_fence_sq_filter = std::views::filter([=](const DistPtElem& dpe) { return dpe.dist < best_dist_plus_fence_sq; });
                    auto dpe_to_pe = std::views::transform([](const DistPtElem& dpe)->node_type {return {*(dpe.pt), *(dpe.elem)}; });
                    std::ranges::copy(std::ranges::subrange(result_dpe_arr, result_dpe_arr_it)
                                        | std::views::reverse
                                        | within_dist_plus_fence_sq_filter
                                        | dpe_to_pe,
                                      pe_out_it);
                    std::ranges::copy(std::views::reverse(result_dpe_vec)
                                        | within_dist_plus_fence_sq_filter
                                        | dpe_to_pe,
                                      pe_out_it);
                    /*
                    std::for_each(std::make_reverse_iterator(result_dist_pe_arr_it), std::make_reverse_iterator(result_dist_pe_arr),
                                  [&pe_out_it, best_dist_plus_fence_sq](const auto& dpe) mutable {
                        if (dpe.dist < best_dist_plus_fence_sq)
                            *pe_out_it++ = {*dpe.pt, *dpe.elem};
                    }); 
                    std::for_each(result_dist_pe_vec.rbegin(), result_dist_pe_vec.rend(), [&pe_out_it, best_dist_plus_fence_sq](const auto& dpe) mutable {
                        if (dpe.dist < best_dist_plus_fence_sq)
                            *pe_out_it++ = {*dpe.pt, *dpe.elem};
                    }); */
                    return pe_out_it;
                }
            } while ((--act_record_it)->dist > best_dist_plus_fence_sq);
            cur = act_record_it->nd;
        }
    }
    return pe_out_it;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N> &pt) const {
    
    struct ActRecord {
        FPType dist;
        const TreeNode* nd;
    } act_record_stack[32], *act_record_it = act_record_stack;
    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    const TreeNode *cur = nd_vec_, *best_node = nd_vec_;
    FPType cur_dist, best_dist = PointND<FPType, N>::template dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
    
    while (true) {
        if (unsigned int right_idx = cur->right_idx) {
            unsigned int dim = cur->dim_to_expand;
            FPType diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(act_record_it++) ActRecord {diff*diff, nd_vec_ + right_idx};
                ++cur;
            } else {
                new(act_record_it++) ActRecord {diff*diff, cur+1};
                cur = nd_vec_ + right_idx;
            }
        } else if (cur++->dim_to_expand == N) {
            do {
                if (act_record_it == act_record_stack)
                    return elem_vec_[best_node - nd_vec_];
            } while ((--act_record_it)->dist >= best_dist);
            cur = act_record_it->nd;
        }
        cur_dist = PointND<FPType, N>::template dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
        if (cur_dist < best_dist) {
            best_dist = cur_dist;
            best_node = cur;
        }
    }
    return elem_vec_[best_node - nd_vec_];
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename PointND<FPType, N>::DistType thisDt,
typename std::enable_if<thisDt == PointND<FPType, N>::DistType::EUC, int>::type>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
NNValueHelper(TreeNode *cur, std::size_t dim, const PointND<FPType, N> &pt,
              const ElemType *&bestValue, FPType &bestDist) const {
    FPType curDist = PointND<FPType, N>::template
    dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    std::size_t nextDim = dim == N - 1 ? 0 : dim + 1;
    FPType diff = pt[dim] - cur->key[dim];
    TreeNode *next = diff < 0.0 ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    if (diff*diff < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    }
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename PointND<FPType, N>::DistType thisDt,
typename std::enable_if<thisDt != PointND<FPType, N>::DistType::EUC, int>::type>
void KDTreeExpandLongestVec<FPType, N, ElemType, DT>::
NNValueHelper(TreeNode *cur, std::size_t dim, const PointND<FPType, N> &pt,
              const ElemType *&bestValue, FPType &bestDist) const {
    FPType curDist = PointND<FPType, N>::template dist<DT>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    std::size_t nextDim = dim == N - 1 ? 0 : dim + 1;
    TreeNode *next = pt[dim] < cur->key[dim] ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    if (branchMin(cur->key, pt, dim) < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    }
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
FPType KDTreeExpandLongestVec<FPType, N, ElemType, DT>::branchMin(const PointND<FPType, N> &trPt,
                                                          const PointND<FPType, N> &searchPt, std::size_t idx) const {
    switch (DT) {
        case PointND<FPType, N>::DistType::EUC:
        case PointND<FPType, N>::DistType::MAN:
            return std::fabs(trPt[idx] - searchPt[idx]);
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
