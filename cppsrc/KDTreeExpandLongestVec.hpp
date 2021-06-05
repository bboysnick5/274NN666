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


template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT
= Point<_Tp, N>::DistType::EUC>
class KDTreeExpandLongestVec {
public:
    
    typedef _Tp                                   value_type;
    typedef ElemType*                             value_iterator;
    typedef const ElemType*                       const_value_iterator;

    
    struct node_type {
        Point<value_type, N> key;
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
    template <typename RAI,
    typename std::enable_if<
    std::is_same<typename std::iterator_traits<RAI>::iterator_category,
    std::random_access_iterator_tag>::value
    && !std::is_const<typename std::remove_pointer<
    typename std::iterator_traits<RAI>::pointer>::type>::value
    && std::is_same<typename std::iterator_traits<RAI>::value_type, node_type>::value,
    int>::type = 0>
    KDTreeExpandLongestVec(RAI, RAI);
    
    template <typename Const_RAI,
    typename std::enable_if<
    std::is_same<typename std::iterator_traits<typename std::remove_const_t<Const_RAI>>::iterator_category,
    std::random_access_iterator_tag>::value
    && std::is_const<typename std::remove_pointer<typename std::iterator_traits<Const_RAI>::pointer>
    ::type>::value
    && std::is_same<typename std::iterator_traits<Const_RAI>::value_type, node_type>::value,
    int>::type = 0>
    KDTreeExpandLongestVec(Const_RAI, Const_RAI);
    
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
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTreeExpandLongestVec.
    constexpr size_t dimension() const;
    typename Point<value_type, N>::DistType distType() const;
    
    // size_t size() const;
    // size_t height() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree, the max height
    // and whether the tree is empty.
    size_t size() const;
    size_t cap() const;
    int height() const;
    bool empty() const;
    
    bool operator==(const KDTreeExpandLongestVec& rhs) const;
    bool operator!=(const KDTreeExpandLongestVec& rhs) const;

    
    void clear();
    void shrink_to_fit();
    
    void printTreeInfo() const;
    
    // bool contains(const Point<_Tp, N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTreeExpandLongestVec.
    bool contains(const Point<value_type, N>&) const;
    
    // void insert(const Point<_Tp, N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTreeExpandLongestVec, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<value_type, N>&, const ElemType&);
    
    // ElemType& operator[](const Point<_Tp, N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTreeExpandLongestVec.
    // If the point does not exist, then it is added to the KDTreeExpandLongestVec using the
    // default value of ElemType as its key.
    ElemType& operator[](const Point<value_type, N>& pt);
    
    // ElemType& at(const Point<_Tp, N>& pt);
    // const ElemType& at(const Point<_Tp, N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function //throws an out_of_range exception.
    ElemType& at(const Point<value_type, N>& pt);
    const ElemType& at(const Point<value_type, N>& pt) const;
    
    // ElemType kNNValue(const Point<_Tp, N>& key, size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTreeExpandLonge*stVec
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<value_type, N>& key, size_t k) const;
    
    // Iter rangeDiffKNNPairs(const Point<_Tp, N>&, _Tp, Iter) const
    // Usage: Iter end = kd.rangeDiffKNNPairs(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a value_type fence IN SQUARE, fill the container iterator result
    // with a list of all points such that every point in the list is at least
    // fence distance closer to p than the rest of the points in the tree,
    // and return the iterator with final value +1 position.
    // The forward iterator is passed in and filled and the end will be returned.
    template <class Forward_Iter>
    Forward_Iter rangeDiffKNNPairs(const Point<value_type, N>& p, value_type fenseSq, Forward_Iter result) const;
    
private:
    
    
    struct TreeNode {
        unsigned int rightIdx;
        unsigned int dimToExpand;
        Point<value_type, N> key;
        bool operator==(const TreeNode&) const = default;
        bool operator!=(const TreeNode&) const = default;
    };
    
    TreeNode *_ndVec;
    ElemType *_objVec;
    uint32_t _size;
    uint32_t _cap;
    
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
    void kNNValueHelper(TreeNode *cur, size_t dim, const Point<value_type, N> &pt,
                        BoundedPQueue<ElemType, value_type> &bpq) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const Point<value_type, N>& key) const;
    
    template <typename Point<value_type, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt == Point<value_type, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, size_t, const Point<value_type, N>&,
                       const ElemType *&, value_type&) const;
    
    template <typename Point<value_type, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt != Point<value_type, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, size_t, const Point<value_type, N>&,
                       const ElemType*&, value_type&) const;
    
    
    // TreeNode** findNodePtr(const Point<_Tp, N>& pt);
    // TreeNode*const* findNodePtr(const Point<_Tp, N>& pt) const;
    // Usage: TreeNode **nodePtr = findNodePtr(pt);
    // ----------------------------------------------------
    // Returns the pointer pointing to the node address
    // corresponding to the given Point. In this _Tp pointing
    // fashion, we can construct a node at that location.
    TreeNode** findNodePtr(const Point<value_type, N>& pt);
    TreeNode*const* findNodePtr(const Point<value_type, N>& pt) const;
    
    value_type branchMin(const Point<value_type, N>&, const Point<value_type, N>&, size_t) const;
    
};

/** KDTreeExpandLongestVec class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::KDTreeExpandLongestVec() {
    _size = _cap = 0;
    _ndVec = nullptr;
    _objVec = nullptr;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs) : _size(rhs._size), _cap(rhs._cap) {
    _ndVec = new TreeNode[_size];
    std::copy_n(rhs._ndVec, _size, _ndVec);
    _objVec = static_cast<ElemType*>(::operator new(_size * sizeof(ElemType), std::nothrow));
    std::uninitialized_copy_n(rhs._objVec, _size, _objVec);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>&
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::operator=(const KDTreeExpandLongestVec& rhs) noexcept {
    if (this != &rhs) {
        deAlloc();
        _size = _cap = rhs._size;
        _ndVec = new TreeNode[_size];
        std::copy_n(rhs._ndVec, _size, _ndVec);
        _objVec = static_cast<ElemType*>(::operator new(_size * sizeof(ElemType)));
        std::uninitialized_copy_n(rhs._objVec, _size, _objVec);
    }
    return *this;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
KDTreeExpandLongestVec(KDTreeExpandLongestVec&& rhs) noexcept
: _ndVec(rhs._ndVec), _objVec(rhs._objVec), _size(rhs._size), _cap(rhs._cap) {
    rhs._ndVec = nullptr;
    rhs._objVec = nullptr;
    rhs._size = rhs._cap = 0;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>& KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
operator=(KDTreeExpandLongestVec&& rhs) noexcept {
    if (this != &rhs) {
        deAlloc();
        _ndVec = rhs._ndVec;
        rhs._ndVec = nullptr;
        _objVec = rhs._objVec;
        rhs._objVec = nullptr;
        _size = rhs._size;
        _cap = rhs._cap;
        rhs._size = rhs._cap = 0;
    }
    return *this;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::~KDTreeExpandLongestVec() {
    deAlloc();
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::deAlloc() {
    delete[] _ndVec;
    std::destroy_n(_objVec, _size);
    ::operator delete(static_cast<void*>(_objVec));
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <typename Const_RAI,
typename std::enable_if<std::is_same<typename std::iterator_traits<typename
std::remove_const_t<Const_RAI>>::iterator_category,
std::random_access_iterator_tag>::value
&& std::is_const<typename std::remove_pointer<typename
std::iterator_traits<Const_RAI>::pointer>::type>::value
&& std::is_same<typename std::iterator_traits<Const_RAI>::value_type,
                typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::node_type>::value,
int>::type>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::KDTreeExpandLongestVec(Const_RAI cbegin, Const_RAI cend)
: _size(static_cast<uint32_t>(cend - cbegin)) {
    if (_size == 0)
        return;
    _ndVec = new TreeNode[_size];
    _objVec = static_cast<ElemType*>(::operator new(_size * sizeof(ElemType), std::nothrow));
    if (_size == 1) {
        *_ndVec = {0, N, cbegin->key};
        new (_objVec) ElemType (cbegin->value);
        return;
    }
    std::vector<node_type> constructData(cbegin, cend);
    auto[lows, highs, spans] = computeInitBBoxSpans(constructData.begin(), constructData.end());
    auto curNd = _ndVec;
    auto curObj = _objVec;
    rangeCtorHelper(curNd, curObj, constructData.begin(), constructData.end(), lows, highs, spans);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <typename RAI, typename std::enable_if<std::is_same<typename
std::iterator_traits<RAI>::iterator_category,
std::random_access_iterator_tag>::value
&& !std::is_const<typename std::remove_pointer<typename std::iterator_traits<RAI>::pointer>::type>::value
&& std::is_same<typename std::iterator_traits<RAI>::value_type, typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::node_type>::value, int>::type>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::KDTreeExpandLongestVec(RAI begin, RAI end) :_size(static_cast<uint32_t>(end - begin)) {
    if (_size == 0)
        return;
    _ndVec = new TreeNode[_size];
    _objVec = static_cast<ElemType*>(::operator new(_size * sizeof(ElemType), std::nothrow));
    if (_size == 1) {
        *_ndVec = {0, N, std::move(begin->key)};
        new (_objVec) ElemType (std::move(begin->value));
        return;
    }
    auto[lows, highs, spans] = computeInitBBoxSpans(begin, end);
    auto curNd = _ndVec;
    auto curObj = _objVec;
    rangeCtorHelper(curNd, curObj, begin, end, lows, highs, spans);
    
    
    /*
    struct actRecord {
        unsigned int ndIdxToAssignRightIdx;
        unsigned int depth;
        RAI thisBeginIt, thisEndIt;
    }; */
    
    /*
    
    size_t stackSize = static_cast<size_t>(log2(_size+1));
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
            //new (actStIt++) actRecord {static_cast<unsigned int>(curNdPtr-_ndVec), depth++, median + 1, thisEndIt};
            new (actStIt++) std::tuple<unsigned int, unsigned int, RAI, RAI> {static_cast<unsigned int>(curNdPtr-_ndVec), depth++, median + 1, thisEndIt};
            thisEndIt = median;
        } else {
            if (numElems == 2) {
                new (curNd++) TreeNode {0, N, thisBeginIt->key};
                new (curObj++) ElemType (thisBeginIt->value);
            }
        TRACEBACK:
            if (curNd - _ndVec == _size)
                return;
            const auto &[ndIdxToAssignRightIdx, stDepth, stThisBeginIt, stThisEndIt] = *--actStIt;
            (_ndVec + ndIdxToAssignRightIdx)->rightIdx = static_cast<unsigned int>(curNd - _ndVec);
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
    std::for_each_n(_ndVec, _size, [](const TreeNode& nd){
        std::cout << nd.dimToExpand << '\t' << nd.rightIdx << '\t';
        std::copy(nd.key.begin(), nd.key.end(), std::ostream_iterator<_Tp>(std::cout,","));
        std::cout << '\n';
    });
    */
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <class Const_RAI>
std::tuple<std::array<_Tp, N>, std::array<_Tp, N>, std::array<_Tp, N>> KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
computeInitBBoxSpans(Const_RAI begin, Const_RAI end) {
    std::array<_Tp, N> lows, highs, spans;
    lows.fill(std::numeric_limits<_Tp>::max());
    highs.fill(std::numeric_limits<_Tp>::min());
    std::for_each(begin, end, [&](const auto &nh) mutable {
        for (size_t i = 0; i < N; ++i) {
            _Tp ptValOnithDim = nh.key[i];
            auto &bboxLow = lows[i], &bboxHigh = highs[i];
            bboxLow = std::min(bboxLow, ptValOnithDim);
            bboxHigh = std::max(bboxHigh, ptValOnithDim);
        }
    });
    std::transform(highs.cbegin(), highs.cend(), lows.cbegin(), spans.begin(), std::minus<_Tp>());
    return {lows, highs, spans};
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
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
    
    if (median + 2 == end) {
        ndPtrThisIter->rightIdx = static_cast<unsigned int>(curNd - _ndVec);
        *curNd++ = {0, N, std::move((median+1)->key)};
        new (curObj++) ElemType (std::move((median+1)->value));
    } else if (median + 1 != end) {
        ndPtrThisIter->rightIdx = static_cast<unsigned int>(curNd - _ndVec);
        _Tp curValOnDim = ndPtrThisIter->key[dim];
        auto prevDimLow = lows[dim];
        lows[dim] = curValOnDim;
        spans[dim] = highs[dim] - curValOnDim;
        rangeCtorHelper(curNd, curObj, median+1, end, lows, highs, spans);
        lows[dim] = prevDimLow;
        spans[dim] = highs[dim] - prevDimLow;
    }
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
constexpr size_t KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::dimension() const {
    return N;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
typename Point<_Tp, N>::DistType KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::distType() const {
    return DT;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
size_t KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::size() const {
    return _size;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
size_t KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::cap() const {
    return _cap;
}


template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
int KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::height() const {
    const TreeNode *root = &_ndVec[0];
    return heightHelper(root);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
int KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::heightHelper(const TreeNode *n) const {
    //return n ? 1 + std::max(heightHelper(n->left), heightHelper(n->right)) : -1;
    return 1;
}


template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
bool KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::empty() const {
    return _size == 0;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
bool KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::operator==(const KDTreeExpandLongestVec& rhs) const {
    if (_size != rhs._size)
        return false;
    return std::equal(_objVec, _objVec + _size, rhs._objVec) && std::equal(_ndVec, _ndVec + _size, rhs._ndVec);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
bool KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::operator!=(const KDTreeExpandLongestVec& rhs) const {
    return !(*this == rhs);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::printTreeInfo() const {
    std::cout << "Tree height is " << height()
    << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::clear() {
    std::destroy_n(_ndVec, _size);
    std::destroy_n(_objVec, _size);
    _size = 0;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::shrink_to_fit() {
    if (_cap > _size) {
        // TO DO
    }
    std::destroy_n(_ndVec, _size);
    std::destroy_n(_objVec, _size);
    _size = 0;
}


template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
insert(const Point<_Tp, N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        //(*ndPtr)->object = value;
    } else {
        //*ndPtr = new TreeNode(0, pt, value);
        //treeSize++;
    }
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
bool KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::contains(const Point<_Tp, N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
ElemType& KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::operator[] (const Point<_Tp, N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        //treeSize++;
    }
    return (*ndPtr)->object;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
ElemType& KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::at(const Point<_Tp, N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeExpandLongestVec&>(*this).at(pt));
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
const ElemType& KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::at(const Point<_Tp, N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::TreeNode**
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::findNodePtr(const Point<_Tp, N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTreeExpandLongestVec*>(this)->findNodePtr(pt));
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
typename KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::TreeNode*const*
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::findNodePtr(const Point<_Tp, N>& pt) const {
    //TreeNode *const*n = &root;
    TreeNode *const*n;
    for (size_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1){}
    //n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
ElemType KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
kNNValue(const Point<_Tp, N>& pt, size_t k) const {
    //if (empty())
    //    return ElemType();
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType, _Tp> bpq(k);
    //kNNValueHelper(ndVec[0], 0, pt, bpq);
    
    std::multimap<size_t, ElemType, std::greater<size_t>> freqMap;
    while (!bpq.empty()) {
        ElemType elem = bpq.dequeueMin();
        for (typename std::multimap<size_t, ElemType>::iterator
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

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::kNNValueHelper(TreeNode *cur, size_t dim,
                                                             const Point<_Tp, N>& pt, BoundedPQueue<ElemType, _Tp> &bpq) const {
    bpq.enqueue(cur->object, Point<_Tp, N>::template dist<DT>(cur->key, pt));
    size_t nextDim = dim + 1 < N ? dim + 1 : 0;
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

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <class Iter>
Iter KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
rangeDiffKNNPairs(const Point<_Tp, N>& pt, _Tp fenceSq, Iter returnIt) const {
    
    struct actRecord {
        _Tp dist;
        const TreeNode* nd;
    } st[static_cast<size_t>(log2(_size+1))], *it = st;
    _Tp bestDistSq = std::numeric_limits<_Tp>::max(),
        bestDistDiffSq = bestDistSq;
    const TreeNode *cur = _ndVec;
    
    const size_t MAX_DISTPTELEMS_ON_STACK = 8192;
    struct distPtElem {
        _Tp dist;
        const Point<_Tp, N>* pt;
        const ElemType* elem;
    } distPtElems[MAX_DISTPTELEMS_ON_STACK], *distPtElemsIt = distPtElems,
      *distPtElemsEnd = distPtElems + MAX_DISTPTELEMS_ON_STACK;
    std::vector<distPtElem> distPtElemVec;
   
    while (true) {
        _Tp curDistSq = Point<_Tp, N>::template dist<Point<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDistSq < bestDistDiffSq) {
            if (curDistSq < bestDistSq) {
                bestDistSq = curDistSq;
                bestDistDiffSq = bestDistSq + fenceSq + 2.0*sqrt(fenceSq*bestDistSq);
            }
            *distPtElemsIt++ = {curDistSq, &cur->key, &_objVec[cur-_ndVec]};
            if (distPtElemsIt == distPtElemsEnd) {
                distPtElemsIt = distPtElems;
                distPtElemVec.insert(distPtElemVec.end(), distPtElems, distPtElemsEnd);
            }
        }
        
        if (unsigned int rightIdx = cur->rightIdx) {
            unsigned int dim = cur->dimToExpand;
            _Tp diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(it++) actRecord {diff*diff, _ndVec + rightIdx};
                ++cur;
            } else {
                new(it++) actRecord {diff*diff, cur+1};
                cur = _ndVec + rightIdx;
            }
        } else if (cur++->dimToExpand == N) {
            do {
                if (it == st)
                    goto FINAL;
            } while ((--it)->dist > bestDistDiffSq);
            cur = it->nd;
        }
    }
    
FINAL:
    std::for_each(distPtElemVec.begin(), distPtElemVec.end(), [&returnIt, bestDistDiffSq](const auto& dpe) mutable {
        if (dpe.dist < bestDistDiffSq)
            *returnIt++ = {*dpe.pt, *dpe.elem};
    });
    std::for_each(distPtElems, distPtElemsIt, [&returnIt, bestDistDiffSq](const auto& dpe) mutable {
        if (dpe.dist < bestDistDiffSq)
            *returnIt++ = {*dpe.pt, *dpe.elem};
    });
    return returnIt;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
ElemType KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::NNValue(const Point<_Tp, N> &pt) const {
    
    struct actRecord {
        _Tp dist;
        const TreeNode* nd;
    } st[static_cast<size_t>(log2(_size+1))], *it = st;
    // BIG ASSUMPTION TREE IS BALANCED, otherwise stackoverflow
    // used cast not std::floor for speed
    const TreeNode *cur = _ndVec, *best = _ndVec;
    _Tp curDist, bestDist = Point<_Tp, N>::template dist<Point<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
    
    while (true) {
        if (unsigned int rightIdx = cur->rightIdx) {
            unsigned int dim = cur->dimToExpand;
            _Tp diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(it++) actRecord {diff*diff, _ndVec + rightIdx};
                ++cur;
            } else {
                new(it++) actRecord {diff*diff, cur+1};
                cur = _ndVec + rightIdx;
            }
        } else if (cur++->dimToExpand == N) {
            do {
                if (it == st)
                    return _objVec[best - _ndVec];
            } while ((--it)->dist >= bestDist);
            cur = it->nd;
        }
        curDist = Point<_Tp, N>::template dist<Point<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDist < bestDist) {
            bestDist = curDist;
            best = cur;
        }
    }
    return _objVec[best - _ndVec];
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <typename Point<_Tp, N>::DistType thisDt,
typename std::enable_if<thisDt == Point<_Tp, N>::DistType::EUC, int>::type>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
NNValueHelper(TreeNode *cur, size_t dim, const Point<_Tp, N> &pt,
              const ElemType *&bestValue, _Tp &bestDist) const {
    _Tp curDist = Point<_Tp, N>::template
    dist<Point<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    size_t nextDim = dim == N - 1 ? 0 : dim + 1;
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

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <typename Point<_Tp, N>::DistType thisDt,
typename std::enable_if<thisDt != Point<_Tp, N>::DistType::EUC, int>::type>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
NNValueHelper(TreeNode *cur, size_t dim, const Point<_Tp, N> &pt,
              const ElemType *&bestValue, _Tp &bestDist) const {
    _Tp curDist = Point<_Tp, N>::template dist<DT>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    size_t nextDim = dim == N - 1 ? 0 : dim + 1;
    TreeNode *next = pt[dim] < cur->key[dim] ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    if (branchMin(cur->key, pt, dim) < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    }
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
_Tp KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::branchMin(const Point<_Tp, N> &trPt,
                                                          const Point<_Tp, N> &searchPt, size_t idx) const {
    switch (DT) {
        case Point<_Tp, N>::DistType::EUC:
        case Point<_Tp, N>::DistType::MAN:
            return std::fabs(trPt[idx] - searchPt[idx]);
            /*
             case DistType::HAV:
             Point<_Tp, N> pt;
             pt[idx] = searchPt[idx];
             pt[1-idx] = trPt[1-idx];
             return Point<_Tp, N>::havDist(trPt, pt);
             */
    }
}




#endif /* KDTreeExpandLongestVec_hpp */
