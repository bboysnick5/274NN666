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
#include <stdexcept>
#include <unordered_map>
#include <map>
#include <utility>
#include <set>
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

    
    // Constructor: KDTreeExpandLongestVec();
    // Usage: KDTreeExpandLongestVec<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTreeExpandLongestVec.
    KDTreeExpandLongestVec() = default;
    
    // Constructor: KDTreeExpandLongestVec(FwdItType begin, FwdItType end);
    // Usage: KDTreeExpandLongestVec<3, int> myTree(vec.begin(), bec.end());
    // ----------------------------------------------------
    // Constructs a KDTreeExpandLongestVec from a collection. The tree will
    // be balanced using median constructing method
    // NOTE: The tree will not eliminate duplicates and the
    //       intended behavior will not be comprimised, tho
    //       less efficient with extra wasteful space.
    template <typename RAI,
    typename std::enable_if<std::is_same<
    typename std::iterator_traits<RAI>::iterator_category,
    std::random_access_iterator_tag>::value &&
    !std::is_const<typename std::remove_pointer<
    typename std::iterator_traits<RAI>::pointer>::type>::value, int>::type = 0>
    KDTreeExpandLongestVec(RAI, RAI);
    
    template <typename Const_RAI,
    typename std::enable_if<std::is_same<
    typename std::iterator_traits<typename
    std::remove_const_t<Const_RAI>>::iterator_category,
    std::random_access_iterator_tag>::value && std::is_const<typename
    std::remove_pointer<typename std::iterator_traits<Const_RAI>::pointer>
    ::type>::value, int>::type = 0>
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
    KDTreeExpandLongestVec& operator=(const KDTreeExpandLongestVec&)&;
    
    // KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs);
    // KDTreeExpandLongestVec& operator=(const KDTreeExpandLongestVec& rhs);
    // Usage: KDTreeExpandLongestVec<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Move constructor and move assignment operator.
    KDTreeExpandLongestVec(KDTreeExpandLongestVec&&) noexcept;
    KDTreeExpandLongestVec& operator=(KDTreeExpandLongestVec&&)& noexcept;
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTreeExpandLongestVec.
    size_t dimension() const;
    typename Point<value_type, N>::DistType distType() const;
    
    // size_t size() const;
    // size_t height() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree, the max height
    // and whether the tree is empty.
    size_t size() const;
    int height() const;
    bool empty() const;
    
    void clear();
    
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
    // Given a point v and an integer k, finds the k points in the KDTreeExpandLongestVec
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<value_type, N>& key, size_t k) const;
    
    // Iter rangeDiffKNNPairs(const Point<_Tp, N>&, _Tp, Iter) const
    // Usage: Iter end = kd.rangeDiffKNNPairs(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a _Tp offset, return a set of points in the KDTreeExpandLongestVec
    // nearest to p such that the farthest one in the set is at least offset
    // distance close to p than the rest of the points in the tree.
    // The forward iterator is passed in and filled and the end will be returned.
    template <class Iter>
    Iter rangeDiffKNNPairs(const Point<value_type, N>&, value_type, Iter) const;
    
private:
    
    
    struct TreeNode {
        unsigned int rightIdx;
        unsigned int dimToExpand;
        Point<value_type, N> key;
    };
    
    TreeNode *_ndVec;
    ElemType *_objVec;
    size_t _size;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    int heightHelper(const TreeNode *n) const;
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void rangeCtorHelper(TreeNode*&, ElemType*&, RAI, RAI, std::array<value_type, N*2>&);
    
    
    template <class RAI>
    static std::array<value_type, N*2> computeInitBBox(RAI, RAI);
    
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
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs) : _size(rhs._size) {
    _ndVec = static_cast<TreeNode*>(::operator new[](_size, std::nothrow));
    _objVec = static_cast<ElemType*>(::operator new[](_size, std::nothrow));
    std::uninitialized_copy_n(rhs._ndVec, _size, _ndVec);
    std::uninitialized_copy_n(rhs._objVec, _size, _objVec);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>&
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::operator=(const KDTreeExpandLongestVec& rhs) & {
    if (this != &rhs) {
        std::destroy_n(_ndVec, _size);
        ::operator delete(_ndVec);
        std::destroy_n(_objVec, _size);
        ::operator delete(_objVec);
        _size = rhs._size;
        _ndVec = static_cast<TreeNode*>(::operator new[](_size, std::nothrow));
        _objVec = static_cast<ElemType*>(::operator new[](_size, std::nothrow));
        std::uninitialized_copy_n(rhs._ndVec, _size, _ndVec);
        std::uninitialized_copy_n(rhs._objVec, _size, _objVec);
    }
    return *this;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
KDTreeExpandLongestVec(KDTreeExpandLongestVec&& rhs) noexcept
: _ndVec(rhs._ndVec), _objVec(rhs._objVec), _size(rhs._size) {
    rhs._ndVec = nullptr;
    rhs._objVec = nullptr;
    rhs._size = 0;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>& KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
operator=(KDTreeExpandLongestVec&& rhs) & noexcept {
    if (this != &rhs) {
        _ndVec = rhs._ndVec;
        rhs._ndVec = nullptr;
        _objVec = rhs._objVec;
        rhs._objVec = nullptr;
        _size = rhs._size;
        rhs._size = 0;
    }
    return *this;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::~KDTreeExpandLongestVec<_Tp, N, ElemType, DT>() {
    std::destroy_n(_ndVec, _size);
    ::operator delete(_ndVec);
    std::destroy_n(_objVec, _size);
    ::operator delete(_objVec);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <typename Const_RAI,
typename std::enable_if<std::is_same<typename std::iterator_traits<typename
std::remove_const_t<Const_RAI>>::iterator_category,
std::random_access_iterator_tag>::value &&
std::is_const<typename std::remove_pointer<typename
std::iterator_traits<Const_RAI>::pointer>::type>::value, int>::type>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::KDTreeExpandLongestVec(Const_RAI cbegin, Const_RAI cend) :
_size(cend - cbegin) {
    _ndVec = static_cast<TreeNode*>(::operator new(_size * sizeof(TreeNode), std::nothrow));
    _objVec = static_cast<ElemType*>(::operator new(_size * sizeof(ElemType), std::nothrow));
    
    std::vector<std::pair<Point<_Tp, N>, ElemType>> constructData(cbegin, cend);
    auto bbox = computeInitBBox(cbegin, cend);
    auto curNd = _ndVec;
    auto curObj = _objVec;
    rangeCtorHelper(curNd, curObj, constructData.begin(), constructData.end(), bbox);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <typename RAI, typename std::enable_if<std::is_same<typename
std::iterator_traits<RAI>::iterator_category,
std::random_access_iterator_tag>::value && !std::is_const<typename
std::remove_pointer< typename std::iterator_traits<RAI>::pointer>::type>::value, int>::type>
KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::KDTreeExpandLongestVec(RAI begin, RAI end) :_size(end - begin) {
    _ndVec = static_cast<TreeNode*>(::operator new(_size * sizeof(TreeNode), std::nothrow));
    _objVec = static_cast<ElemType*>(::operator new(_size * sizeof(ElemType), std::nothrow));
    auto bbox = computeInitBBox(begin, end);
    auto curNd = _ndVec;
    auto curObj = _objVec;
    rangeCtorHelper(curNd, curObj, begin, end, bbox);
    
    
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
        int dim = 0;
        _Tp maxSpan = bbox[1] - bbox[0];
        for (int i = 1; i < N; ++i) {
            int bboxLowIdx = i * 2, bboxHighIdx = bboxLowIdx + 1;
            _Tp span = bbox[bboxHighIdx] - bbox[bboxLowIdx];
            if (span > maxSpan) {
                maxSpan = span;
                dim = i;
            }
        }
        
        TreeNode* curNdPtr = curNd;
        RAI median = thisBeginIt + (thisEndIt-thisBeginIt)/2;
        std::nth_element(thisBeginIt, median, thisEndIt,
                         [=](const auto& p1, const auto& p2) {
                             return p1.first[dim] < p2.first[dim];});
        new (curNd++) TreeNode {0, dim, std::move(median->first)};
        new (curObj++) ElemType (std::move(median->second));
     
        auto span = thisEndIt - thisBeginIt;
        if (span > 2) {
            _Tp *bboxValToChangeLeftPtr = bbox.data() + dim*2+1;
            *++curBboxChangeStIt = {bboxValToChangeLeftPtr, *bboxValToChangeLeftPtr};
            *bboxValToChangeLeftPtr = curNdPtr->key[dim];
            //new (actStIt++) actRecord {static_cast<unsigned int>(curNdPtr-_ndVec), depth++, median + 1, thisEndIt};
            new (actStIt++) std::tuple<unsigned int, unsigned int, RAI, RAI> {static_cast<unsigned int>(curNdPtr-_ndVec), depth++, median + 1, thisEndIt};
            thisEndIt = median;
        } else {
            if (span == 2) {
                new (curNd++) TreeNode {0, -1, std::move(thisBeginIt->first)};
                new (curObj++) ElemType (std::move(thisBeginIt->second));
            }
        TRACEBACK:
            if (curNd - _ndVec == _size)
                return;
            const auto &[ndIdxToAssignRightIdx, stDepth, stThisBeginIt, stThisEndIt] = *--actStIt;
            (_ndVec+ ndIdxToAssignRightIdx)->rightIdx = static_cast<unsigned int>(curNd - _ndVec);
            if ((thisBeginIt = stThisBeginIt) + 1 == (thisEndIt = stThisEndIt)) {
                new (curNd++) TreeNode {0, -1, std::move(thisBeginIt->first)};
                new (curObj++) ElemType (std::move(thisBeginIt->second));
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
template <class RAI>
std::array<_Tp, N*2> KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
computeInitBBox(RAI begin, RAI end) {
    std::array<_Tp, N*2> bbox;
    std::generate(bbox.begin(), bbox.end(),
                  [lowHighToggle = 0lu, lowHigh = std::array<_Tp, 2>{
                   std::numeric_limits<_Tp>::max(),
                   std::numeric_limits<_Tp>::min()}]() mutable {
                       return lowHigh[lowHighToggle++%2];});
    std::for_each(begin, end, [&](const auto &p) mutable {
        for (size_t i = 0; i < N; ++i) {
            _Tp ptValOnithDim = p.first[i];
            auto &bboxLow = bbox[i*2], &bboxHigh = bbox[i*2+1];
            bboxLow = std::min(bboxLow, ptValOnithDim);
            bboxHigh = std::max(bboxHigh, ptValOnithDim);
        }
    });
    return bbox;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <class RAI>
void KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::
rangeCtorHelper(TreeNode *&curNd, ElemType *&curObj, RAI begin, RAI end,
                std::array<_Tp, N*2> &bbox) {
    unsigned int dim = 0;
    _Tp maxSpan = bbox[1] - bbox[0];
    for (unsigned int i = 1; i != N; ++i) {
        unsigned int bboxLowIdx = i * 2, bboxHighIdx = bboxLowIdx + 1;
        _Tp span = bbox[bboxHighIdx] - bbox[bboxLowIdx];
        if (span > maxSpan) {
            maxSpan = span;
            dim = i;
        }
    }
    
    TreeNode* curNdPtr = curNd;
    RAI median = begin + (end - begin)/2;
    std::nth_element(begin, median, end, [=](const auto& p1, const auto& p2) {
        return p1.first[dim] < p2.first[dim];});
    new (curNd++) TreeNode {0, dim, median->first};
    new (curObj++) ElemType (median->second);
    _Tp curValOnDim = curNdPtr->key[dim];
    _Tp *bboxChangePtr = bbox.data() + dim*2+1;
    
    if (begin == median - 1) {
        new (curNd++) TreeNode {0, N, begin->first};
        new (curObj++) ElemType (begin->second);
    } else if (begin != median) {
        auto prevDimHigh = *bboxChangePtr;
        *bboxChangePtr = curValOnDim;
        rangeCtorHelper(curNd, curObj, begin, median, bbox);
        *bboxChangePtr = prevDimHigh;
    }
    
    if (median + 2 == end) {
        curNdPtr->rightIdx = static_cast<unsigned int>(curNd - _ndVec);
        new (curNd++) TreeNode {0, N, (median+1)->first};
        new (curObj++) ElemType ((median+1)->second);
    } else if (median + 1 != end) {
        curNdPtr->rightIdx = static_cast<unsigned int>(curNd - _ndVec);
        auto prevDimLow = *--bboxChangePtr;
        *bboxChangePtr = curValOnDim;
        rangeCtorHelper(curNd, curObj, median+1, end, bbox);
        *bboxChangePtr = prevDimLow;
    }
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
size_t KDTreeExpandLongestVec<_Tp, N, ElemType, DT>::dimension() const {
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
rangeDiffKNNPairs(const Point<_Tp, N>& pt, _Tp fence, Iter returnIt) const {
    
    
    struct actRecord {
        _Tp dist;
        const TreeNode* nd;
    } st[static_cast<size_t>(log2(_size+1))], *it = st;
    _Tp bestDistSq = std::numeric_limits<_Tp>::max(),
           bestDistDiffSq = bestDistSq, fenceSq = fence*fence;
    const TreeNode *cur = _ndVec;
    
    const size_t MAX_DISTPTELEMS_SIZE = 16384;
    struct distPtElem {
        _Tp dist;
        const Point<_Tp, N>* pt;
        const ElemType* elem;
    } distPtElems[MAX_DISTPTELEMS_SIZE], *distPtElemsIt = distPtElems,
      *distPtElemsNextToEnd = distPtElems + MAX_DISTPTELEMS_SIZE - 1;
    std::vector<distPtElem> distPtElemVec;
   
    while (true) {
        _Tp curDistSq = Point<_Tp, N>::template
        dist<Point<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDistSq < bestDistDiffSq) {
            if (curDistSq < bestDistSq) {
                bestDistSq = curDistSq;
                bestDistDiffSq = bestDistSq + fenceSq + 2*fence*sqrt(bestDistSq);
            }
            *distPtElemsIt++ = {curDistSq, &cur->key, &_objVec[cur-_ndVec]};
            if (distPtElemsIt == distPtElemsNextToEnd) {
                distPtElemsIt = distPtElems;
                distPtElemVec.insert(distPtElemVec.end(), distPtElems, distPtElemsNextToEnd);
                distPtElemVec.push_back({curDistSq, &cur->key, &_objVec[cur-_ndVec]});
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
            } while ((--it)->dist > bestDistDiffSq || !(cur = it->nd));
        }
    }
    
FINAL:
    std::for_each(distPtElemVec.begin(), distPtElemVec.end(), [&returnIt, bestDistDiffSq](const auto& dpe){
        if (dpe.dist < bestDistDiffSq)
            *returnIt++ = {*dpe.pt, *dpe.elem};
    });
    std::for_each(distPtElems, distPtElemsIt, [&returnIt, bestDistDiffSq](const auto& dpe){
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
    _Tp bestDist = std::numeric_limits<_Tp>::max();
    const TreeNode *cur = _ndVec, *best = cur;
    
    // LOGGGGGGGGGGGGGGG
    
    /*
     static size_t totalNumNodesSearches = 0, numNNSearches = 0, totalTreeSize = 0;
     static size_t numOfFullSearch = 0;
     size_t thisNumNodesSearches = 0;
     static bool logCondition;
     static constexpr size_t TREE_SIZE_LOWER_BOUND = 1200, TREE_SIZE_UPPER_BOUND = 60000000;
     logCondition = treeSize <= TREE_SIZE_UPPER_BOUND && treeSize >= TREE_SIZE_LOWER_BOUND;
     
     */
    
    while (true) {
        _Tp curDist = Point<_Tp, N>::template dist<Point<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDist < bestDist) {
            bestDist = curDist;
            best = cur;
        }
        
        // LOGGGGGGGGGGGGGGG
        //if (logCondition)
        // thisNumNodesSearches++;
        
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
            } while ((--it)->dist >= bestDist || !(cur = it->nd));
        }
    }
    
    // LOGGGGGGGGGGGGGGG
    /*
     if (logCondition) {
     numNNSearches++;
     totalTreeSize += treeSize;
     totalNumNodesSearches += thisNumNodesSearches;
     if (thisNumNodesSearches == treeSize)
     numOfFullSearch++;
     std::cout << "***** Log of treesize from " << TREE_SIZE_LOWER_BOUND
     << " to " << TREE_SIZE_UPPER_BOUND << " *****"
     << "\nTotal num of NN searches with this criteria: " << numNNSearches
     << "\n\nTreesize: " << treeSize
     << "\nNum of nodes searched in this NN search: " << thisNumNodesSearches
     << "\n\nAve treesize: " << totalTreeSize/numNNSearches
     << "\nAve num of nodes searched in each NN search: " << totalNumNodesSearches/numNNSearches
     << "\n\nnum of full searches percentage: " << numOfFullSearch*100.0/numNNSearches
     << "%\nSearched nodes over total num of nodes percentage: " << totalNumNodesSearches*100.0/totalTreeSize << "%\n\n\n\n";
     
     }
     */
    
    /*
     _Tp bestDist = std::numeric_limits<_Tp>::max();
     const ElemType *bestValue = nullptr;
     NNValueHelper(root, 0, pt, bestValue, bestDist);
     */
    
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
