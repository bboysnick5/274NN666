//
//  KDTreeExpandLongest.hpp
//  274F16NearestSB
//
//  Created by nick on 1/18/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef KDTreeExpandLongest_hpp
#define KDTreeExpandLongest_hpp


#include "Point.hpp"
#include "PoolAllocator.hpp"
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


template <size_t N, typename ElemType, typename Point<N>::DistType DT
= Point<N>::DistType::EUC>
class KDTreeExpandLongest {
public:
    // Constructor: KDTreeExpandLongest();
    // Usage: KDTreeExpandLongest<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTreeExpandLongest.
    KDTreeExpandLongest() = default;
    
    // Constructor: KDTreeExpandLongest(FwdItType begin, FwdItType end);
    // Usage: KDTreeExpandLongest<3, int> myTree(vec.begin(), bec.end());
    // ----------------------------------------------------
    // Constructs a KDTreeExpandLongest from a collection. The tree will
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
    KDTreeExpandLongest(RAI, RAI);
    
    template <typename Const_RAI,
    typename std::enable_if<std::is_same<
    typename std::iterator_traits<typename
    std::remove_const_t<Const_RAI>>::iterator_category,
    std::random_access_iterator_tag>::value && std::is_const<typename
    std::remove_pointer<typename std::iterator_traits<Const_RAI>::pointer>
    ::type>::value, int>::type = 0>
    KDTreeExpandLongest(Const_RAI, Const_RAI);
    
    // Destructor: ~KDTreeExpandLongest()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTreeExpandLongest.
    ~KDTreeExpandLongest();
    
    // KDTreeExpandLongest(const KDTreeExpandLongest& rhs);
    // KDTreeExpandLongest& operator=(const KDTreeExpandLongest& rhs);
    // Usage: KDTreeExpandLongest<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Copy constructor and copy assignment operator.
    KDTreeExpandLongest(const KDTreeExpandLongest&);
    KDTreeExpandLongest& operator=(const KDTreeExpandLongest&)&;
    
    // KDTreeExpandLongest(const KDTreeExpandLongest& rhs);
    // KDTreeExpandLongest& operator=(const KDTreeExpandLongest& rhs);
    // Usage: KDTreeExpandLongest<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Move constructor and move assignment operator.
    KDTreeExpandLongest(KDTreeExpandLongest&&) noexcept;
    KDTreeExpandLongest& operator=(KDTreeExpandLongest&&)& noexcept;
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTreeExpandLongest.
    size_t dimension() const;
    typename Point<N>::DistType distType() const;
    
    // size_t size() const;
    // size_t height() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree, the max height
    // and whether the tree is empty.
    size_t size() const;
    size_t height() const;
    bool empty() const;
    
    void clear();
    
    void printTreeInfo() const;
    
    // bool contains(const Point<N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTreeExpandLongest.
    bool contains(const Point<N>&) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTreeExpandLongest, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>&, const ElemType&);
    
    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTreeExpandLongest.
    // If the point does not exist, then it is added to the KDTreeExpandLongest using the
    // default value of ElemType as its key.
    ElemType& operator[](const Point<N>& pt);
    
    // ElemType& at(const Point<N>& pt);
    // const ElemType& at(const Point<N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function //throws an out_of_range exception.
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;
    
    // ElemType kNNValue(const Point<N>& key, size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTreeExpandLongest
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, size_t k) const;
    
    // Iter rangeDiffKNNPairs(const Point<N>&, double, Iter) const
    // Usage: Iter end = kd.rangeDiffKNNPairs(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a double offset, return a set of points in the KDTreeExpandLongest
    // nearest to p such that the farthest one in the set is at least offset
    // distance close to p than the rest of the points in the tree.
    // The forward iterator is passed in and filled and the end will be returned.
    template <class Iter>
    Iter rangeDiffKNNPairs(const Point<N>&, double, Iter) const;
    
private:
    
    struct TreeNode {
        size_t dimToExpand;
        Point<N> key;
        TreeNode *left;
        TreeNode *right;
        ElemType object;
        
        TreeNode() = default;
        // TreeNode(const TreeNode&) = default;
        // TreeNode& operator = (const TreeNode&) = default;
        // TreeNode(TreeNode&&) = default;
        // TreeNode& operator = (TreeNode&&) = default;
        
        TreeNode(size_t dimToExpand, const Point<N>& k, const ElemType& obj)
        : dimToExpand(dimToExpand), key(k), left(nullptr), right(nullptr),
          object(obj) {}
        
        TreeNode(size_t dimToExpand, const Point<N>&& k, const ElemType&& obj)
        : dimToExpand(dimToExpand), key(std::move(k)), left(nullptr), right(nullptr),
          object(std::move(obj)) {}
        
        ~TreeNode() = default;
    };
    
    
    TreeNode *root;
    size_t treeSize;
    std::unique_ptr<PooledAllocator> pool;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    size_t heightHelper(TreeNode *n) const;
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void rangeCtorHelper(TreeNode*&, RAI, RAI, RAI, std::array<double, N*2>&);
    
    
    template <class RAI>
    static std::array<double, N*2> computeInitBBox(RAI, RAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(TreeNode *cur, size_t dim, const Point<N> &pt,
                        BoundedPQueue<ElemType> &bpq) const;
    
    void rangeDiffKNNPairsHelper(TreeNode*, size_t, const Point<N>&, double,
                               std::vector<std::pair<double,
                            std::pair<Point<N>, ElemType>>>&, double&, double&) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const Point<N>& key) const;
    
    template <typename Point<N>::DistType thisDt = DT,
    typename std::enable_if<thisDt == Point<N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, size_t, const Point<N>&,
                       const ElemType *&, double&) const;
    
    template <typename Point<N>::DistType thisDt = DT,
    typename std::enable_if<thisDt != Point<N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, size_t, const Point<N>&,
                       const ElemType*&, double&) const;
    
    // ----------------------------------------------------
    // Helper meothod for deep copy
    void treeCopy(TreeNode*& thisNd, TreeNode *otherNd //TreeNode* ndPoolPtr
    );
    
    // TreeNode** findNodePtr(const Point<N>& pt);
    // TreeNode*const* findNodePtr(const Point<N>& pt) const;
    // Usage: TreeNode **nodePtr = findNodePtr(pt);
    // ----------------------------------------------------
    // Returns the pointer pointing to the node address
    // corresponding to the given Point. In this double pointing
    // fashion, we can construct a node at that location.
    TreeNode** findNodePtr(const Point<N>& pt);
    TreeNode*const* findNodePtr(const Point<N>& pt) const;
    
    double branchMin(const Point<N>&, const Point<N>&, size_t) const;
    
};

/** KDTreeExpandLongest class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------


template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <typename Const_RAI,
typename std::enable_if<std::is_same<typename std::iterator_traits<typename
std::remove_const_t<Const_RAI>>::iterator_category,
std::random_access_iterator_tag>::value &&
std::is_const<typename std::remove_pointer<typename
std::iterator_traits<Const_RAI>::pointer>::type>::value, int>::type>
KDTreeExpandLongest<N, ElemType, DT>::KDTreeExpandLongest(Const_RAI cbegin, Const_RAI cend)
: treeSize(cend-cbegin), pool(std::make_unique<PooledAllocator>()) {
    std::vector<std::pair<Point<N>, ElemType>> constructData(cbegin, cend);
    TreeNode* ndPoolPtr = root = pool->allocateExact<TreeNode>(treeSize);
    auto bbox = computeInitBBox(cbegin, cend);
    rangeCtorHelper(ndPoolPtr, constructData.begin(), constructData.end(),
                    constructData.begin() +
                    (constructData.end() - constructData.begin())/2, bbox);
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <typename RAI, typename std::enable_if<std::is_same<typename
std::iterator_traits<RAI>::iterator_category,
std::random_access_iterator_tag>::value && !std::is_const<typename
std::remove_pointer< typename std::iterator_traits<RAI>::pointer>::type>::value, int>::type>
KDTreeExpandLongest<N, ElemType, DT>::KDTreeExpandLongest(RAI begin, RAI end)
//, std::array<double, N> bboxHint)
: treeSize(end-begin), pool(std::make_unique<PooledAllocator>()) {
    auto bbox = computeInitBBox(begin, end);
    TreeNode* ndPoolPtr = root = pool->allocateExact<TreeNode>(treeSize);
    rangeCtorHelper(ndPoolPtr, begin, begin + (end - begin)/2, end, bbox);
    
    /*
     struct actRecord {
     TreeNode** curNdPtr;
     size_t dim;
     RAI thisBeginIt, median, thisEndIt;
     };
     
     actRecord st[static_cast<size_t>(log2(treeSize+1))], *it = st;
     bool hasChild = true;
     TreeNode* curNd;
     size_t dim = 0;
     RAI thisBeginIt = begin, thisEndIt = end, median = thisBeginIt + (thisEndIt-thisBeginIt)/2;
     while (it != st || hasChild) {
     if (!hasChild) {
     hasChild = true;
     *(--it)->curNdPtr = ndPoolPtr;
     dim = it->dim;
     thisBeginIt = it->thisBeginIt;
     thisEndIt = it->thisEndIt;
     median = it->median;
     }
     
     std::nth_element(thisBeginIt, median, thisEndIt,
     [=](const auto& p1, const auto& p2) {
     return p1.first[dim] < p2.first[dim];});
     pool->construct(ndPoolPtr, std::move(median->first),
     std::move(median->second));
     curNd = ndPoolPtr++;
     dim = dim == N - 1 ? 0 : dim + 1;
     
     if (thisBeginIt != median) {
     if (median + 1 != thisEndIt) {
     *it++ = {&curNd->right, dim, median + 1,
     median + (thisEndIt-median+1)/2, thisEndIt};
     }
     thisEndIt = median;
     median = thisBeginIt + (median-thisBeginIt)/2;
     curNd->left = ndPoolPtr;
     } else if (median + 1 != thisEndIt) {
     thisBeginIt = median+1;
     median = median + (thisEndIt-median+1)/2;
     curNd->right = ndPoolPtr;
     } else {
     hasChild = false;
     }
     }
     */
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTreeExpandLongest<N, ElemType, DT>::KDTreeExpandLongest(const KDTreeExpandLongest& rhs)
: root(new TreeNode()), treeSize(rhs.treeSize) {
    // wrong logic.
    // should be check whether this size is greater than other.
    // if not allocate additional space. if yes
    //std::destroy_n(root, treeSize);
    //pool->free_all();
    //root = pool->allocate<TreeNode>(rhs.treeSize);
    //auto *ndPoolIt = addressof(*root);
    //pool->construct(root, {rhs.root->key, rhs.root->obj});
    treeCopy(root, rhs.root);
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTreeExpandLongest<N, ElemType, DT>&
KDTreeExpandLongest<N, ElemType, DT>::operator=(const KDTreeExpandLongest& rhs) & {
    if (this != &rhs) {
        delete root;
        treeSize = rhs.treeSize;
        root = new TreeNode();
        treeCopy(root, rhs.root);
    }
    return *this;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTreeExpandLongest<N, ElemType, DT>::KDTreeExpandLongest(KDTreeExpandLongest&& rhs) noexcept
: root(rhs.root), treeSize(rhs.treeSize), pool(std::move(rhs.pool)) {
    rhs.root = nullptr;
    rhs.treeSize = 0;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTreeExpandLongest<N, ElemType, DT>& KDTreeExpandLongest<N, ElemType, DT>::
operator=(KDTreeExpandLongest&& rhs) & noexcept {
    if (this != &rhs) {
        root = rhs.root;
        treeSize = rhs.treeSize;
        pool = std::move(rhs.pool);
        rhs.root = nullptr;
        rhs.treeSize = 0;
    }
    return *this;
}

/*
 template <size_t N, typename ElemType, typename Point<N>::DistType DT>
 template <class RAI>
 void KDTreeExpandLongest<N, ElemType, DT>::
 rangeCtorHelper(TreeNode*& ndPoolPtr, RAI begin,
 RAI median, RAI end, std::array<double, N>& bbox) {
 size_t dim = std::max_element(bbox.begin(), bbox.end()) - bbox.begin();
 std::nth_element(begin, median, end, [=](const auto& p1, const auto& p2) {
 return p1.first[dim] < p2.first[dim];});
 pool->construct(ndPoolPtr, std::move(median->first),
 std::move(median->second));
 auto curNdPtr = ndPoolPtr;
 //size_t nextDim = dim == N - 1 ? 0 : dim + 1;
 double dec;
 
 if (begin == median - 1) {
 curNdPtr->left = ++ndPoolPtr;
 pool->construct(ndPoolPtr, std::move(begin->first),
 std::move(begin->second));
 } else if (begin != median) {
 dec = (end-1)->first[dim] - (median-1)->first[dim];
 bbox[dim] -= dec;
 //nextDim = std::max_element(bbox, bbox+N)-bbox;
 curNdPtr->left = ++ndPoolPtr;
 rangeCtorHelper(ndPoolPtr, begin,
 begin + (median-begin)/2, median, bbox);
 bbox[dim] += dec;
 }
 
 if (median + 2 == end) {
 curNdPtr->right = ++ndPoolPtr;
 pool->construct(ndPoolPtr, std::move((median + 1)->first),
 std::move((median+1)->second));
 } else if (median + 1 != end) {
 dec = (median+1)->first[dim] - begin->first[dim];
 bbox[dim] -= dec;
 curNdPtr->right = ++ndPoolPtr;
 rangeCtorHelper(ndPoolPtr, median+1,
 median + (end-median+1)/2, end, bbox);
 bbox[dim] += dec;
 }
 }
 */

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <class RAI>
std::array<double, N*2> KDTreeExpandLongest<N, ElemType, DT>::
computeInitBBox(RAI begin, RAI end) {
    std::array<double, N*2> bbox;
    std::generate(bbox.begin(), bbox.end(),
                  [lowHighToggle = 0lu, lowHigh = std::array<double, 2>{
                   std::numeric_limits<double>::max(),
                   std::numeric_limits<double>::min()}]() mutable {
                      return lowHigh[lowHighToggle++%2];});
    std::for_each(begin, end, [&](const auto &p) mutable {
        for (size_t i = 0; i < N; ++i) {
            size_t bboxLowIdx = i * 2, bboxHighIdx = bboxLowIdx + 1;
            double ptValOnithDim = p.first[i];
            bbox[bboxLowIdx] = std::min(bbox[bboxLowIdx], ptValOnithDim);
            bbox[bboxHighIdx] = std::max(bbox[bboxHighIdx], ptValOnithDim);
        }
    });
    return bbox;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <class RAI>
void KDTreeExpandLongest<N, ElemType, DT>::
rangeCtorHelper(TreeNode*& ndPoolPtr, RAI begin, RAI median, RAI end,
                std::array<double, N*2> &bbox) {
    size_t dim = 0;
    double maxSpan = bbox[1] - bbox[0];
    for (size_t i = 1; i < N; ++i) {
        size_t bboxLowIdx = i * 2, bboxHighIdx = bboxLowIdx + 1;
        double span = bbox[bboxHighIdx] - bbox[bboxLowIdx];
        if (span > maxSpan) {
            dim = i;
        }
    }
    
    std::nth_element(begin, median, end, [=](const auto& p1, const auto& p2) {
        return p1.first[dim] < p2.first[dim];});
    pool->construct(ndPoolPtr, dim, std::move(median->first),
                    std::move(median->second));
    auto curNdPtr = ndPoolPtr;
    
    if (begin == median - 1) {
        curNdPtr->left = ++ndPoolPtr;
        pool->construct(ndPoolPtr, dim, std::move(begin->first),
                        std::move(begin->second));
    } else if (begin != median) {
        curNdPtr->left = ++ndPoolPtr;
        auto prevDimHigh = bbox[dim*2+1];
        bbox[dim*2+1] = median->first[dim];
        rangeCtorHelper(ndPoolPtr, begin, begin + (median-begin)/2, median, bbox);
        bbox[dim*2+1] = prevDimHigh;
    }
    
    if (median + 2 == end) {
        curNdPtr->right = ++ndPoolPtr;
        pool->construct(ndPoolPtr, dim, std::move((median + 1)->first),
                        std::move((median+1)->second));
    } else if (median + 1 != end) {
        curNdPtr->right = ++ndPoolPtr;
        auto prevDimLow = bbox[dim*2];
        bbox[dim*2] = median->first[dim];
        rangeCtorHelper(ndPoolPtr, median+1, median + (end-median+1)/2, end, bbox);
        bbox[dim*2] = prevDimLow;
    }
}


template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTreeExpandLongest<N, ElemType, DT>::treeCopy(TreeNode*& thisNode,
                                             TreeNode* otherNode
//TreeNode* ndPoolIt
) {
    if (otherNode) {
        if (thisNode) {
            thisNode->object = otherNode->object;
            thisNode->key = otherNode->key;
        } else {
            //pool->construct(thisNode, {otherNode->key, otherNode->object});
            //++ndPoolIt;
            thisNode = new TreeNode(otherNode->key, otherNode->object);
        }
        treeCopy(thisNode->left, otherNode->left);
        treeCopy(thisNode->right, otherNode->right);
    } else if (thisNode) {
        delete thisNode;
        thisNode = nullptr;
    }
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTreeExpandLongest<N, ElemType, DT>::~KDTreeExpandLongest<N, ElemType, DT>() {
    if (pool)
        pool->destroy_and_free_all<TreeNode>();
}

// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTreeExpandLongest<N, ElemType, DT>::dimension() const {
    return N;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
typename Point<N>::DistType KDTreeExpandLongest<N, ElemType, DT>::distType() const {
    return DT;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTreeExpandLongest<N, ElemType, DT>::size() const {
    return treeSize;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTreeExpandLongest<N, ElemType, DT>::height() const {
    return heightHelper(root);
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTreeExpandLongest<N, ElemType, DT>::heightHelper(TreeNode *n) const {
    return n ? 1 + std::max(heightHelper(n->left),
                            heightHelper(n->right)) : -1;
}


template <size_t N, typename ElemType, typename Point<N>::DistType DT>
bool KDTreeExpandLongest<N, ElemType, DT>::empty() const {
    return treeSize == 0;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTreeExpandLongest<N, ElemType, DT>::printTreeInfo() const {
    std::cout << "Tree height is " << height()
    << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTreeExpandLongest<N, ElemType, DT>::clear() {
    pool->destroy_and_free_all<TreeNode>();
    root = nullptr;
    treeSize = 0;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTreeExpandLongest<N, ElemType, DT>::
insert(const Point<N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        (*ndPtr)->object = value;
    } else {
        *ndPtr = new TreeNode(0, pt, value);
        treeSize++;
    }
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
bool KDTreeExpandLongest<N, ElemType, DT>::contains(const Point<N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType& KDTreeExpandLongest<N, ElemType, DT>::operator[] (const Point<N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        treeSize++;
    }
    return (*ndPtr)->object;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType& KDTreeExpandLongest<N, ElemType, DT>::at(const Point<N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeExpandLongest&>(*this).at(pt));
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
const ElemType& KDTreeExpandLongest<N, ElemType, DT>::at(const Point<N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
typename KDTreeExpandLongest<N, ElemType, DT>::TreeNode**
KDTreeExpandLongest<N, ElemType, DT>::findNodePtr(const Point<N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTreeExpandLongest*>(this)->findNodePtr(pt));
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
typename KDTreeExpandLongest<N, ElemType, DT>::TreeNode*const*
KDTreeExpandLongest<N, ElemType, DT>::findNodePtr(const Point<N>& pt) const {
    TreeNode *const*n = &root;
    for (size_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1)
        n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType KDTreeExpandLongest<N, ElemType, DT>::
kNNValue(const Point<N>& pt, size_t k) const {
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType> bpq(k);
    kNNValueHelper(root, 0, pt, bpq);
    
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

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTreeExpandLongest<N, ElemType, DT>::kNNValueHelper(TreeNode *cur, size_t dim,
                                                   const Point<N>& pt, BoundedPQueue<ElemType> &bpq) const {
    bpq.enqueue(cur->object, Point<N>::template dist<DT>(cur->key, pt));
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

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <class Iter>
Iter KDTreeExpandLongest<N, ElemType, DT>::
rangeDiffKNNPairs(const Point<N>& pt, double fence, Iter returnIt) const {
    /*
     std::vector<std::pair<double, std::pair<Point<N>, ElemType>>> distKVPairs;
     distKVPairs.reserve(sqrt(treeSize));
     double bestSq = std::numeric_limits<double>::max(),
     bestDiffSq = std::numeric_limits<double>::max();
     rangeDiffKNNPairsHelper(root, 0, pt, fence, distKVPairs, bestSq, bestDiffSq);
     for (const auto &p : distKVPairs) {
     if (p.first < bestDiffSq)
     *returnIt++ = std::move(p.second);
     }
     return returnIt;
     */
    
    std::vector<std::tuple<double, const Point<N>&, const ElemType&>> distPtElemTuple;
    distPtElemTuple.reserve(sqrt(treeSize));
    double bestDistSq = std::numeric_limits<double>::max(),
    bestDistDiffSq = std::numeric_limits<double>::max(), fenceSq = fence*fence;
    const TreeNode *cur = root;
    bool hasNext = true;
    std::pair<double, const TreeNode*> st[static_cast<size_t>(log2(treeSize+1))],
                                       *it = st;
    while (it != st || hasNext) {
        if (!hasNext) {
            if ((--it)->first >= bestDistDiffSq || !(cur = it->second))
                continue;
            hasNext = true;
        }
        double curDistSq = Point<N>::template
                           dist<Point<N>::DistType::EUCSQ>(cur->key, pt);
        if (curDistSq < bestDistDiffSq) {
            if (curDistSq < bestDistSq) {
                bestDistSq = curDistSq;
                bestDistDiffSq = bestDistSq + fenceSq + 2*fence*sqrt(bestDistSq);
            }
            distPtElemTuple.emplace_back(curDistSq, cur->key, cur->object);
        }
        size_t dim = cur->dimToExpand;
        double diff = pt[dim] - cur->key[dim];
        const TreeNode* next = diff < 0.0 ? cur->left : cur->right;
        diff *= diff;
        if (next) {
            //*it++ = std::tie(diff, next == cur->left ? cur->right : cur->left);
            //auto p = std::pair{diff, next == cur->left ? cur->right : cur->left};
            *it++ = {diff, next == cur->left ? cur->right : cur->left};
           // st[it++-st] = {diff, next == cur->left ? cur->right : cur->left};
            cur = next;
        } else {
            if (diff < bestDistDiffSq) {
                next = next == cur->left ? cur->right : cur->left;
                if (next) {
                    cur = next;
                    continue;
                }
            }
            hasNext = false;
        }
    }
    
    for (const auto &[distSq, pt, elem] : distPtElemTuple) {
        if (distSq < bestDistDiffSq)
            *returnIt++ = {pt, elem};
    }
    return returnIt;
    
}


 template <size_t N, typename ElemType, typename Point<N>::DistType DT>
 void KDTreeExpandLongest<N, ElemType, DT>::
 rangeDiffKNNPairsHelper(TreeNode *cur, size_t dim, const Point<N>& pt,
 double diff, std::vector<std::pair<double,
 std::pair<Point<N>, ElemType>>> &distKVPairs,
 double& bestDistSq, double&bestDistDiffSq) const {
 auto distSq = Point<N>::template dist<Point<N>::DistType::EUCSQ>(cur->key, pt);
 if (distSq < bestDistDiffSq) {
 if (distSq < bestDistSq) {
 bestDistSq = distSq;
 bestDistDiffSq = bestDistSq + diff*diff + 2*diff*sqrt(bestDistSq);
 }
 distKVPairs.emplace_back(distSq, std::make_pair(cur->key, cur->object));
 }
 size_t nextDim = dim + 1 < N ? dim + 1 : 0;
 double thisDiff = pt[dim] - cur->key[dim];
 TreeNode *next = thisDiff < 0 ? cur->left : cur->right;
 if (next)
 rangeDiffKNNPairsHelper(next, nextDim, pt, diff, distKVPairs, bestDistSq, bestDistDiffSq);
 if (thisDiff*thisDiff < bestDistDiffSq) {
 next = next == cur->left ? cur->right : cur->left;
 if (next)
 rangeDiffKNNPairsHelper(next, nextDim, pt, diff, distKVPairs, bestDistSq, bestDistDiffSq);
 }
 }

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType KDTreeExpandLongest<N, ElemType, DT>::NNValue(const Point<N> &pt) const {
    
    
    double bestDist = std::numeric_limits<double>::max(), curDist, diff;
    size_t dim;
    const ElemType *bestValue = nullptr;
    TreeNode *cur = root, *next;
    bool hasNext = true;
    std::pair<double, TreeNode*> st[static_cast<size_t>(log2(treeSize+1))], *it = st;
    
    
    // LOGGGGGGGGGGGGGGG
    
    /*
     static size_t totalNumNodesSearches = 0, numNNSearches = 0, totalTreeSize = 0;
     static size_t numOfFullSearch = 0;
     size_t thisNumNodesSearches = 0;
     static bool logCondition;
     static constexpr size_t TREE_SIZE_LOWER_BOUND = 1200, TREE_SIZE_UPPER_BOUND = 60000;
     logCondition = treeSize <= TREE_SIZE_UPPER_BOUND && treeSize >= TREE_SIZE_LOWER_BOUND;
     */
    
    
    while (it != st || hasNext) {
        if (!hasNext) {
            if ((--it)->first >= bestDist || !(cur = it->second))
                continue;
            hasNext = true;
        }
        curDist = Point<N>::template dist<Point<N>::DistType::EUCSQ>(cur->key, pt);
        if (curDist < bestDist) {
            bestDist = curDist;
            bestValue = &cur->object;
        }
        
        // LOGGGGGGGGGGGGGGG
        // if (logCondition)
       //   thisNumNodesSearches++;
        
        dim = cur->dimToExpand;
        diff = pt[dim] - cur->key[dim];
        next = diff < 0.0 ? cur->left : cur->right;
        curDist = diff*diff;
        if (next) {
            //*it++ = std::tie(curDist, next == cur->left ? cur->right : cur->left);
            *it++ = {curDist, next == cur->left ? cur->right : cur->left};
            cur = next;
        } else {
            if (curDist < bestDist) {
                next = next == cur->left ? cur->right : cur->left;
                if (next) {
                    cur = next;
                    continue;
                }
            }
            hasNext = false;
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
     double bestDist = std::numeric_limits<double>::max();
     const ElemType *bestValue = nullptr;
     NNValueHelper(root, 0, pt, bestValue, bestDist);
     */
    
    return bestValue ? *bestValue : ElemType();
    
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <typename Point<N>::DistType thisDt,
typename std::enable_if<thisDt == Point<N>::DistType::EUC, int>::type>
void KDTreeExpandLongest<N, ElemType, DT>::
NNValueHelper(TreeNode *cur, size_t dim, const Point<N> &pt,
              const ElemType *&bestValue, double &bestDist) const {
    double curDist = Point<N>::template
    dist<Point<N>::DistType::EUCSQ>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    size_t nextDim = dim == N - 1 ? 0 : dim + 1;
    double diff = pt[dim] - cur->key[dim];
    TreeNode *next = diff < 0.0 ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    if (diff*diff < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, nextDim, pt, bestValue, bestDist);
    }
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <typename Point<N>::DistType thisDt,
typename std::enable_if<thisDt != Point<N>::DistType::EUC, int>::type>
void KDTreeExpandLongest<N, ElemType, DT>::
NNValueHelper(TreeNode *cur, size_t dim, const Point<N> &pt,
              const ElemType *&bestValue, double &bestDist) const {
    double curDist = Point<N>::template dist<DT>(cur->key, pt);
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

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
double KDTreeExpandLongest<N, ElemType, DT>::branchMin(const Point<N> &trPt,
                                                const Point<N> &searchPt, size_t idx) const {
    switch (DT) {
        case Point<N>::DistType::EUC:
        case Point<N>::DistType::MAN:
            return std::fabs(trPt[idx] - searchPt[idx]);
            /*
             case DistType::HAV:
             Point<N> pt;
             pt[idx] = searchPt[idx];
             pt[1-idx] = trPt[1-idx];
             return Point<N>::havDist(trPt, pt);
             */
    }
}





#endif /* KDTreeExpandLongest_hpp */
