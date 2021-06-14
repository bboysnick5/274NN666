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


template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
class KDTreeExpandLongestVec {
public:
    
    typedef _Tp                                   value_type;
    typedef ElemType*                             value_iterator;
    typedef const ElemType*                       const_value_iterator;

    
    struct node_type {
        PointND<value_type, N> key;
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
    template <std::random_access_iterator RAI>
    requires std::same_as<typename std::iterator_traits<RAI>::value_type, typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::node_type>
    KDTreeExpandLongestVec(RAI, RAI);
    
    template <std::random_access_iterator ConstRAI> requires def::const_iterator<ConstRAI>
    && std::same_as<typename std::iterator_traits<ConstRAI>::value_type, typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::node_type const>
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
    
    // std::size_t dimension() const;
    // Usage: std::size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTreeExpandLongestVec.
    constexpr std::size_t dimension() const;
    typename PointND<value_type, N>::DistType distType() const;
    
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
    bool empty() const;
    
    bool operator==(const KDTreeExpandLongestVec& rhs) const;
    bool operator!=(const KDTreeExpandLongestVec& rhs) const;

    
    void clear();
    void shrink_to_fit();
    
    void printTreeInfo() const;
    
    // bool contains(const PointND<_Tp, N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTreeExpandLongestVec.
    bool contains(const PointND<value_type, N>&) const;
    
    // void insert(const PointND<_Tp, N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTreeExpandLongestVec, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const PointND<value_type, N>&, const ElemType&);
    
    // ElemType& operator[](const PointND<_Tp, N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTreeExpandLongestVec.
    // If the point does not exist, then it is added to the KDTreeExpandLongestVec using the
    // default value of ElemType as its key.
    ElemType& operator[](const PointND<value_type, N>& pt);
    
    // ElemType& at(const PointND<_Tp, N>& pt);
    // const ElemType& at(const PointND<_Tp, N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function //throws an out_of_range exception.
    ElemType& at(const PointND<value_type, N>& pt);
    const ElemType& at(const PointND<value_type, N>& pt) const;
    
    // ElemType kNNValue(const PointND<_Tp, N>& key, std::size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTreeExpandLonge*stVec
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const PointND<value_type, N>& key, std::size_t k) const;
    
    // Iter rangeDiffKNNPairs(const PointND<_Tp, N>&, _Tp, Iter) const
    // Usage: Iter end = kd.rangeDiffKNNPairs(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a value_type fence IN SQUARE, fill the container iterator result
    // with a list of all points such that every point in the list is at least
    // fence distance closer to p than the rest of the points in the tree,
    // and return the iterator with final value +1 position.
    // The forward iterator is passed in and filled and the end will be returned.
    template <class Forward_Iter>
    Forward_Iter rangeDiffKNNPairs(const PointND<value_type, N>& p, value_type fenseSq, Forward_Iter result) const;
    
private:
    
    
    struct TreeNode {
        unsigned int rightIdx;
        unsigned int dimToExpand;
        PointND<value_type, N> key;
        bool operator==(const TreeNode&) const = default;
        bool operator!=(const TreeNode&) const = default;
    };
    
    TreeNode *ndVec_;
    ElemType *objVec_;
    uint32_t size_;
    uint32_t cap_;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    int heightHelper(const TreeNode *n) const;
    
    void deAlloc();
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void rangeCtorHelper(TreeNode*&, ElemType*&, RAI, RAI, std::array<value_type, N>&,
                         std::array<value_type, N>&, std::array<value_type, N>&);
    
    
    template <class Const_RAI>
    static std::tuple<std::array<value_type, N>, std::array<value_type,N>, std::array<value_type, N>>
           computeInitBBoxSpans(Const_RAI, Const_RAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(TreeNode *cur, std::size_t dim, const PointND<value_type, N> &pt,
                        BoundedPQueue<ElemType, value_type> &bpq) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const PointND<value_type, N>& key) const;
    
    template <typename PointND<value_type, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt == PointND<value_type, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, std::size_t, const PointND<value_type, N>&,
                       const ElemType *&, value_type&) const;
    
    template <typename PointND<value_type, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt != PointND<value_type, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, std::size_t, const PointND<value_type, N>&,
                       const ElemType*&, value_type&) const;
    
    
    // TreeNode** findNodePtr(const PointND<_Tp, N>& pt);
    // TreeNode*const* findNodePtr(const PointND<_Tp, N>& pt) const;
    // Usage: TreeNode **nodePtr = findNodePtr(pt);
    // ----------------------------------------------------
    // Returns the pointer pointing to the node address
    // corresponding to the given PointND. In this _Tp pointing
    // fashion, we can construct a node at that location.
    TreeNode** findNodePtr(const PointND<value_type, N>& pt);
    TreeNode*const* findNodePtr(const PointND<value_type, N>& pt) const;
    
    value_type branchMin(const PointND<value_type, N>&, const PointND<value_type, N>&, std::size_t) const;
    
};

/** KDTreeExpandLongestVec class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::KDTreeExpandLongestVec() {
    size_ = cap_ = 0;
    ndVec_ = nullptr;
    objVec_ = nullptr;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs) : size_(rhs.size_), cap_(rhs.cap_) {
    ndVec_ = new TreeNode[size_];
    std::copy_n(rhs.ndVec_, size_, ndVec_);
    objVec_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType), std::nothrow));
    std::uninitialized_copy_n(rhs.objVec_, size_, objVec_);
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>&
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::operator=(const KDTreeExpandLongestVec& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        deAlloc();
        size_ = cap_ = rhs.size_;
        ndVec_ = new TreeNode[size_];
        std::copy_n(rhs.ndVec_, size_, ndVec_);
        objVec_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType)));
        std::uninitialized_copy_n(rhs.objVec_, size_, objVec_);
    }
    return *this;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
KDTreeExpandLongestVec(KDTreeExpandLongestVec&& rhs) noexcept
: ndVec_(rhs.ndVec_), objVec_(rhs.objVec_), size_(rhs.size_), cap_(rhs.cap_) {
    rhs.ndVec_ = nullptr;
    rhs.objVec_ = nullptr;
    rhs.size_ = rhs.cap_ = 0;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>& KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
operator=(KDTreeExpandLongestVec&& rhs) noexcept {
    if (this != &rhs) [[likely]] {
        deAlloc();
        ndVec_ = rhs.ndVec_;
        rhs.ndVec_ = nullptr;
        objVec_ = rhs.objVec_;
        rhs.objVec_ = nullptr;
        size_ = rhs.size_;
        cap_ = rhs.cap_;
        rhs.size_ = rhs.cap_ = 0;
    }
    return *this;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::~KDTreeExpandLongestVec() {
    deAlloc();
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::deAlloc() {
    delete[] ndVec_;
    std::destroy_n(objVec_, size_);
    ::operator delete(static_cast<void*>(objVec_));
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
template <std::random_access_iterator ConstRAI> requires def::const_iterator<ConstRAI>
&& std::same_as<typename std::iterator_traits<ConstRAI>::value_type, typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::node_type const>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::KDTreeExpandLongestVec(ConstRAI cbegin, ConstRAI cend)
: size_(static_cast<uint32_t>(cend - cbegin)) {
    if (size_ == 0) [[unlikely]]
        return;
    ndVec_ = new TreeNode[size_];
    objVec_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType), std::nothrow));
    if (size_ == 1) [[unlikely]] {
        *ndVec_ = {0, N, cbegin->key};
        new (objVec_) ElemType (cbegin->value);
        return;
    }
    std::vector<node_type> constructData(cbegin, cend);
    auto[lows, highs, spans] = computeInitBBoxSpans(constructData.begin(), constructData.end());
    auto curNd = ndVec_;
    auto curObj = objVec_;
    rangeCtorHelper(curNd, curObj, constructData.begin(), constructData.end(), lows, highs, spans);
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
template <std::random_access_iterator RAI> 
requires std::same_as<typename std::iterator_traits<RAI>::value_type, typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::node_type>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::KDTreeExpandLongestVec(RAI begin, RAI end) :size_(static_cast<uint32_t>(end - begin)) {
    if (size_ == 0) [[unlikely]]
        return;
    ndVec_ = new TreeNode[size_];
    objVec_ = static_cast<ElemType*>(::operator new(size_ * sizeof(ElemType), std::nothrow));
    if (size_ == 1) [[unlikely]] {
        *ndVec_ = {0, N, std::move(begin->key)};
        new (objVec_) ElemType (std::move(begin->value));
        return;
    }
    auto[lows, highs, spans] = computeInitBBoxSpans(begin, end);
    auto curNd = ndVec_;
    auto curObj = objVec_;
    rangeCtorHelper(curNd, curObj, begin, end, lows, highs, spans);
    
    
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
    std::pair<_Tp*, _Tp> bboxChangeSt[stackSize], *curBboxChangeStIt = bboxChangeSt-1;
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
            //new (actStIt++) actRecord {static_cast<unsigned int>(curNdPtr-ndVec_), depth++, median + 1, thisEndIt};
            new (actStIt++) std::tuple<unsigned int, unsigned int, RAI, RAI> {static_cast<unsigned int>(curNdPtr-ndVec_), depth++, median + 1, thisEndIt};
            thisEndIt = median;
        } else {
            if (numElems == 2) {
                new (curNd++) TreeNode {0, N, thisBeginIt->key};
                new (curObj++) ElemType (thisBeginIt->value);
            }
        TRACEBACK:
            if (curNd - ndVec_ == size_)
                return;
            const auto &[ndIdxToAssignRightIdx, stDepth, stThisBeginIt, stThisEndIt] = *--actStIt;
            (ndVec_ + ndIdxToAssignRightIdx)->rightIdx = static_cast<unsigned int>(curNd - ndVec_);
            if ((thisBeginIt = stThisBeginIt) + 1 == (thisEndIt = stThisEndIt)) {
                new (curNd++) TreeNode {0, N, thisBeginIt->key};
                new (curObj++) ElemType (thisBeginIt->value);
                goto TRACEBACK;
            }
            
            
            while (depth != stDepth + 1 && --depth ) {
                *curBboxChangeStIt->first = curBboxChangeStIt->second;
                --curBboxChangeStIt;
            }
           
            _Tp rightVal = *(curBboxChangeStIt->first-1);
            *(curBboxChangeStIt->first-1) = *curBboxChangeStIt->first;
            *curBboxChangeStIt->first = curBboxChangeStIt->second;
            *curBboxChangeStIt = {curBboxChangeStIt->first-1, rightVal};
        }
    }
    
     */
     
    /*
    std::for_each_n(ndVec_, size_, [](const TreeNode& nd){
        std::cout << nd.dimToExpand << '\t' << nd.rightIdx << '\t';
        std::copy(nd.key.begin(), nd.key.end(), std::ostream_iterator<_Tp>(std::cout,","));
        std::cout << '\n';
    });
    */
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
template <class Const_RAI>
std::tuple<std::array<_Tp, N>, std::array<_Tp, N>, std::array<_Tp, N>> KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
computeInitBBoxSpans(Const_RAI begin, Const_RAI end) {
    std::array<_Tp, N> lows, highs, spans;
    lows.fill(std::numeric_limits<_Tp>::max());
    highs.fill(std::numeric_limits<_Tp>::min());
    std::for_each(begin, end, [&](const auto &nh) mutable {
        for (std::size_t i = 0; i < N; ++i) {
            _Tp ptValOnithDim = nh.key[i];
            auto &bboxLow = lows[i], &bboxHigh = highs[i];
            bboxLow = std::min(bboxLow, ptValOnithDim);
            bboxHigh = std::max(bboxHigh, ptValOnithDim);
        }
    });
    std::transform(highs.cbegin(), highs.cend(), lows.cbegin(), spans.begin(), std::minus<_Tp>());
    return {lows, highs, spans};
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
template <class RAI>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
rangeCtorHelper(TreeNode *&curNd, ElemType *&curObj, RAI begin, RAI end,
                std::array<_Tp, N> &lows, std::array<_Tp, N> &highs, std::array<_Tp, N> &spans) {
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
        _Tp curValOnDim = ndPtrThisIter->key[dim];
        auto prevDimHigh = highs[dim];
        highs[dim] = curValOnDim;
        spans[dim] = curValOnDim - lows[dim];
        rangeCtorHelper(curNd, curObj, begin, median, lows, highs, spans);
        highs[dim] = prevDimHigh;
        spans[dim] = prevDimHigh - lows[dim];
    }
    
    if (median + 1 != end) {
        if (median + 2 == end) {
            ndPtrThisIter->rightIdx = static_cast<unsigned int>(curNd - ndVec_);
            *curNd++ = {0, N, std::move((median+1)->key)};
            new (curObj++) ElemType (std::move((median+1)->value));
        } else {
            ndPtrThisIter->rightIdx = static_cast<unsigned int>(curNd - ndVec_);
            _Tp curValOnDim = ndPtrThisIter->key[dim];
            auto prevDimLow = lows[dim];
            lows[dim] = curValOnDim;
            spans[dim] = highs[dim] - curValOnDim;
            rangeCtorHelper(curNd, curObj, median+1, end, lows, highs, spans);
            lows[dim] = prevDimLow;
            spans[dim] = highs[dim] - prevDimLow;
        }
    }
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
constexpr std::size_t KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::dimension() const {
    return N;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
typename PointND<_Tp, N>::DistType KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::distType() const {
    return DT;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
std::size_t KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::size() const {
    return size_;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
std::size_t KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::cap() const {
    return cap_;
}


template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
int KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::height() const {
    const TreeNode *root = &ndVec_[0];
    return heightHelper(root);
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
int KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::heightHelper(const TreeNode *n) const {
    //return n ? 1 + std::max(heightHelper(n->left), heightHelper(n->right)) : -1;
    return 1;
}


template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
bool KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::empty() const {
    return size_ == 0;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
bool KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::operator==(const KDTreeExpandLongestVec& rhs) const {
    if (size_ != rhs.size_)
        return false;
    return std::equal(objVec_, objVec_ + size_, rhs.objVec_) && std::equal(ndVec_, ndVec_ + size_, rhs.ndVec_);
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
bool KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::operator!=(const KDTreeExpandLongestVec& rhs) const {
    return !(*this == rhs);
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::printTreeInfo() const {
    std::cout << "Tree height is " << height()
    << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::clear() {
    std::destroy_n(ndVec_, size_);
    std::destroy_n(objVec_, size_);
    size_ = 0;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::shrink_to_fit() {
    if (cap_ > size_) {
        // TO DO
    }
    std::destroy_n(ndVec_, size_);
    std::destroy_n(objVec_, size_);
    size_ = 0;
}


template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
insert(const PointND<_Tp, N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        //(*ndPtr)->object = value;
    } else {
        //*ndPtr = new TreeNode(0, pt, value);
        //treeSize++;
    }
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
bool KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::contains(const PointND<_Tp, N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
ElemType& KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::operator[] (const PointND<_Tp, N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        //treeSize++;
    }
    return (*ndPtr)->object;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
ElemType& KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::at(const PointND<_Tp, N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeExpandLongestVec&>(*this).at(pt));
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
const ElemType& KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::at(const PointND<_Tp, N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::TreeNode**
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::findNodePtr(const PointND<_Tp, N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTreeExpandLongestVec*>(this)->findNodePtr(pt));
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::TreeNode*const*
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::findNodePtr(const PointND<_Tp, N>& pt) const {
    //TreeNode *const*n = &root;
    TreeNode *const*n;
    for (std::size_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1){}
    //n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
ElemType KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
kNNValue(const PointND<_Tp, N>& pt, std::size_t k) const {
    //if (empty())
    //    return ElemType();
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType, _Tp> bpq(k);
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

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::kNNValueHelper(TreeNode *cur, std::size_t dim,
                                                             const PointND<_Tp, N>& pt, BoundedPQueue<ElemType, _Tp> &bpq) const {
    bpq.enqueue(cur->object, PointND<_Tp, N>::template dist<DT>(cur->key, pt));
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

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
template <class OutputIter>
OutputIter KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
rangeDiffKNNPairs(const PointND<_Tp, N>& pt, _Tp fence_sq, OutputIter out_nd_type_it) const {
    
    struct ActRecord {
        _Tp dist;
        const TreeNode* nd;
    } act_record_stack[32], *act_record_it = act_record_stack;
    _Tp best_dist_sq = std::numeric_limits<_Tp>::max(),
        best_dist_plus_fence_sq = best_dist_sq;
    const TreeNode *cur = ndVec_;
    
    const std::size_t MAX_DISTPTELEMS_ON_STACK = 1152;
    struct DistKV {
        _Tp dist;
        const PointND<_Tp, N>* pt;
        const ElemType* elem;
    } result_dist_kv_arr[MAX_DISTPTELEMS_ON_STACK], *result_dist_kv_arr_it = result_dist_kv_arr,
      *result_dist_kv_arr_end = result_dist_kv_arr + MAX_DISTPTELEMS_ON_STACK;
    std::vector<DistKV> result_dist_kv_vec;
   
    while (true) {
        _Tp curDistSq = PointND<_Tp, N>::template dist<PointND<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDistSq < best_dist_plus_fence_sq) {
            if (curDistSq < best_dist_sq) {
                best_dist_sq = curDistSq;
                best_dist_plus_fence_sq = best_dist_sq + fence_sq + _Tp(2.0)*sqrt(fence_sq*best_dist_sq);
            }
            *result_dist_kv_arr_it++ = {curDistSq, &cur->key, &objVec_[cur-ndVec_]};
            if (result_dist_kv_arr_it == result_dist_kv_arr_end) [[unlikely]] {
                result_dist_kv_arr_it = result_dist_kv_arr;
                result_dist_kv_vec.insert(result_dist_kv_vec.end(), result_dist_kv_arr, result_dist_kv_arr_end);
            }
        }
        
        if (unsigned int rightIdx = cur->rightIdx) {
            unsigned int dim = cur->dimToExpand;
            _Tp diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(act_record_it++) ActRecord {diff*diff, ndVec_ + rightIdx};
                ++cur;
            } else {
                new(act_record_it++) ActRecord {diff*diff, cur+1};
                cur = ndVec_ + rightIdx;
            }
        } else if (cur++->dimToExpand == N) {
            do {
                if (act_record_it == act_record_stack) [[unlikely]] {
                    std::for_each(result_dist_kv_vec.begin(), result_dist_kv_vec.end(), [&out_nd_type_it, best_dist_plus_fence_sq](const auto& dpe) mutable {
                        if (dpe.dist < best_dist_plus_fence_sq)
                            *out_nd_type_it++ = {*dpe.pt, *dpe.elem};
                    });
                    std::for_each(result_dist_kv_arr, result_dist_kv_arr_it, [&out_nd_type_it, best_dist_plus_fence_sq](const auto& dpe) mutable {
                        if (dpe.dist < best_dist_plus_fence_sq)
                            *out_nd_type_it++ = {*dpe.pt, *dpe.elem};
                    });
                    return out_nd_type_it;
                }
            } while ((--act_record_it)->dist > best_dist_plus_fence_sq);
            cur = act_record_it->nd;
        }
    }
    return out_nd_type_it;
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
ElemType KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::NNValue(const PointND<_Tp, N> &pt) const {
    
    struct act_record {
        _Tp dist;
        const TreeNode* nd;
    } act_record_stack[32], *act_record_it = act_record_stack;
    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    // used cast not std::floor for speed
    const TreeNode *cur = ndVec_, *best_node = ndVec_;
    _Tp cur_dist, best_dist = PointND<_Tp, N>::template dist<PointND<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
    
    while (true) {
        if (unsigned int right_idx = cur->rightIdx) {
            unsigned int dim = cur->dimToExpand;
            _Tp diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(act_record_it++) act_record {diff*diff, ndVec_ + right_idx};
                ++cur;
            } else {
                new(act_record_it++) act_record {diff*diff, cur+1};
                cur = ndVec_ + right_idx;
            }
        } else if (cur++->dimToExpand == N) {
            do {
                if (act_record_it == act_record_stack)
                    return objVec_[best_node - ndVec_];
            } while ((--act_record_it)->dist >= best_dist);
            cur = act_record_it->nd;
        }
        cur_dist = PointND<_Tp, N>::template dist<PointND<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
        if (cur_dist < best_dist) {
            best_dist = cur_dist;
            best_node = cur;
        }
    }
    return objVec_[best_node - ndVec_];
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
template <typename PointND<_Tp, N>::DistType thisDt,
typename std::enable_if<thisDt == PointND<_Tp, N>::DistType::EUC, int>::type>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
NNValueHelper(TreeNode *cur, std::size_t dim, const PointND<_Tp, N> &pt,
              const ElemType *&bestValue, _Tp &bestDist) const {
    _Tp curDist = PointND<_Tp, N>::template
    dist<PointND<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    std::size_t nextDim = dim == N - 1 ? 0 : dim + 1;
    _Tp diff = pt[dim] - cur->key[dim];
    TreeNode *next = diff < 0.0 ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    if (diff*diff < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    }
}

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
template <typename PointND<_Tp, N>::DistType thisDt,
typename std::enable_if<thisDt != PointND<_Tp, N>::DistType::EUC, int>::type>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
NNValueHelper(TreeNode *cur, std::size_t dim, const PointND<_Tp, N> &pt,
              const ElemType *&bestValue, _Tp &bestDist) const {
    _Tp curDist = PointND<_Tp, N>::template dist<DT>(cur->key, pt);
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

template <typename _Tp, std::size_t N, typename ElemType, typename PointND<_Tp, N>::DistType DT>
_Tp KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::branchMin(const PointND<_Tp, N> &trPt,
                                                          const PointND<_Tp, N> &searchPt, std::size_t idx) const {
    switch (DT) {
        case PointND<_Tp, N>::DistType::EUC:
        case PointND<_Tp, N>::DistType::MAN:
            return std::fabs(trPt[idx] - searchPt[idx]);
            /*
             case DistType::HAV:
             PointND<_Tp, N> pt;
             pt[idx] = searchPt[idx];
             pt[1-idx] = trPt[1-idx];
             return PointND<_Tp, N>::havDist(trPt, pt);
             */
    }
}




#endif /* KDTreeExpandLongestVec_hpp */
