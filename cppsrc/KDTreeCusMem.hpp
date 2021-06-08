//
//  KDTreeCusMemCusMem.hpp
//  274F16NearestSB
//
//  Created by nick on 1/12/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef KDTreeCusMem_hpp
#define KDTreeCusMem_hpp


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


template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT
= Point<_Tp, N>::DistType::EUC>
class KDTreeCusMem {
public:
    
    typedef _Tp                                   value_type;
    
    struct node_type {
        Point<value_type, N> key;
        ElemType value;
    };
    
    // Constructor: KDTreeCusMem();
    // Usage: KDTreeCusMem<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTreeCusMem.
    KDTreeCusMem() = default;
    
    // Constructor: KDTreeCusMem(FwdItType begin, FwdItType end);
    // Usage: KDTreeCusMem<3, int> myTree(vec.begin(), bec.end());
    // ----------------------------------------------------
    // Constructs a KDTreeCusMem from a collection. The tree will
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
    KDTreeCusMem(RAI, RAI);
    
    template <typename Const_RAI,
        typename std::enable_if<std::is_same<
        typename std::iterator_traits<typename
        std::remove_const_t<Const_RAI>>::iterator_category,
        std::random_access_iterator_tag>::value && std::is_const<typename
        std::remove_pointer<typename std::iterator_traits<Const_RAI>::pointer>
           ::type>::value, int>::type = 0>
    KDTreeCusMem(Const_RAI, Const_RAI);

    // Destructor: ~KDTreeCusMem()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTreeCusMem.
    ~KDTreeCusMem();
    
    // KDTreeCusMem(const KDTreeCusMem& rhs);
    // KDTreeCusMem& operator=(const KDTreeCusMem& rhs);
    // Usage: KDTreeCusMem<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Copy constructor and copy assignment operator.
    KDTreeCusMem(const KDTreeCusMem&);
    KDTreeCusMem& operator=(const KDTreeCusMem&)&;
    
    // KDTreeCusMem(const KDTreeCusMem& rhs);
    // KDTreeCusMem& operator=(const KDTreeCusMem& rhs);
    // Usage: KDTreeCusMem<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Move constructor and move assignment operator.
    KDTreeCusMem(KDTreeCusMem&&) noexcept;
    KDTreeCusMem& operator=(KDTreeCusMem&&)& noexcept;
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTreeCusMem.
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
    int height() const;
    bool empty() const;
    
    void clear();
    
    void printTreeInfo() const;
    
    // bool contains(const Point<_Tp, N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTreeCusMem.
    bool contains(const Point<value_type, N>&) const;
    
    // void insert(const Point<_Tp, N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTreeCusMem, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<value_type, N>&, const ElemType&);
    
    // ElemType& operator[](const Point<_Tp, N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTreeCusMem.
    // If the point does not exist, then it is added to the KDTreeCusMem using the
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
    // Given a point v and an integer k, finds the k points in the KDTreeCusMem
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<value_type, N>& key, size_t k) const;
    
    // Iter rangeDiffKNNPairs(const Point<_Tp, N>&, _Tp, Iter) const
    // Usage: Iter end = kd.rangeDiffKNNPairs(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a _Tp offset, return a set of points in the KDTreeCusMem
    // nearest to p such that the farthest one in the set is at least offset
    // distance close to p than the rest of the points in the tree.
    // The forward iterator is passed in and filled and the end will be returned.
    template <class OutputIter>
    OutputIter rangeDiffKNNPairs(const Point<value_type, N>&, value_type, OutputIter) const;
    
private:
    
    struct TreeNode {
        Point<value_type, N> key;
        TreeNode *left;
        TreeNode *right;
        ElemType object;
        
        TreeNode() = default;
        // TreeNode(const TreeNode&) = default;
        // TreeNode& operator = (const TreeNode&) = default;
        // TreeNode(TreeNode&&) = default;
        // TreeNode& operator = (TreeNode&&) = default;
        
        TreeNode(const Point<value_type, N>& k, const ElemType& obj)
        : key(k), left(nullptr), right(nullptr), object(obj) {}
        
        TreeNode(Point<value_type, N>&& k, ElemType&& obj)
        : key(std::move(k)), left(nullptr), right(nullptr),
        object(std::move(obj)) {}
        
        ~TreeNode() = default;
    };
    
    
    TreeNode *root;
    size_t treeSize;
    std::unique_ptr<PooledAllocator> pool;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    int heightHelper(TreeNode *n) const;
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void rangeCtorHelper(TreeNode*&, size_t, RAI, RAI, RAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(TreeNode *cur, size_t dim, const Point<value_type, N> &pt,
                        BoundedPQueue<ElemType, value_type> &bpq) const;
    
    //void rangeDiffKNNPairsHelper(TreeNode*, size_t, const Point<_Tp, N>&, _Tp,
    //                            std::vector<std::pair<_Tp,
    //                         std::pair<Point<_Tp, N>, ElemType>>>&, _Tp&, _Tp&) const;
    
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
    
    // ----------------------------------------------------
    // Helper meothod for deep copy
    void treeCopy(TreeNode*& thisNd, TreeNode *otherNd //TreeNode* ndPoolPtr
    );
    
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

/** KDTreeCusMem class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------


template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <typename Const_RAI,
    typename std::enable_if<std::is_same<typename std::iterator_traits<typename
    std::remove_const_t<Const_RAI>>::iterator_category,
    std::random_access_iterator_tag>::value &&
    std::is_const<typename std::remove_pointer<typename
    std::iterator_traits<Const_RAI>::pointer>::type>::value, int>::type>
KDTreeCusMem<_Tp, N, ElemType, DT>::KDTreeCusMem(Const_RAI cbegin, Const_RAI cend)
: treeSize(cend-cbegin), pool(std::make_unique<PooledAllocator>()) {
    std::vector<node_type> constructData(cbegin, cend);
    TreeNode* ndPoolPtr = root = pool->allocateExact<TreeNode>(treeSize);
    rangeCtorHelper(ndPoolPtr, 0, constructData.begin(), constructData.end(),
                    constructData.begin() +
                    (constructData.end() - constructData.begin())/2);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <typename RAI, typename std::enable_if<std::is_same<typename
    std::iterator_traits<RAI>::iterator_category,
    std::random_access_iterator_tag>::value && !std::is_const<typename
    std::remove_pointer< typename std::iterator_traits<RAI>::pointer>::type>::value, int>::type>
KDTreeCusMem<_Tp, N, ElemType, DT>::KDTreeCusMem(RAI begin, RAI end)
                                            //, std::array<_Tp, N> bboxHint)
: treeSize(end-begin), pool(std::make_unique<PooledAllocator>()) {
    TreeNode* ndPoolPtr = root = pool->allocateExact<TreeNode>(treeSize);
    rangeCtorHelper(ndPoolPtr, 0, begin, begin + (end - begin)/2, end);
    
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

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeCusMem<_Tp, N, ElemType, DT>::KDTreeCusMem(const KDTreeCusMem& rhs)
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

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeCusMem<_Tp, N, ElemType, DT>&
KDTreeCusMem<_Tp, N, ElemType, DT>::operator=(const KDTreeCusMem& rhs) & {
    if (this != &rhs) {
        delete root;
        treeSize = rhs.treeSize;
        root = new TreeNode();
        treeCopy(root, rhs.root);
    }
    return *this;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeCusMem<_Tp, N, ElemType, DT>::KDTreeCusMem(KDTreeCusMem&& rhs) noexcept
: root(rhs.root), treeSize(rhs.treeSize), pool(std::move(rhs.pool)) {
    rhs.root = nullptr;
    rhs.treeSize = 0;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeCusMem<_Tp, N, ElemType, DT>& KDTreeCusMem<_Tp, N, ElemType, DT>::
operator=(KDTreeCusMem&& rhs) & noexcept {
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
template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <class RAI>
void KDTreeCusMem<_Tp, N, ElemType, DT>::
rangeCtorHelper(TreeNode*& ndPoolPtr, RAI begin,
                RAI median, RAI end, std::array<_Tp, N>& bbox) {
    size_t dim = std::max_element(bbox.begin(), bbox.end()) - bbox.begin();
    std::nth_element(begin, median, end, [=](const auto& p1, const auto& p2) {
        return p1.first[dim] < p2.first[dim];});
    pool->construct(ndPoolPtr, std::move(median->first),
                    std::move(median->second));
    auto curNdPtr = ndPoolPtr;
    //size_t nextDim = dim == N - 1 ? 0 : dim + 1;
    _Tp dec;
    
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

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <class RAI>
void KDTreeCusMem<_Tp, N, ElemType, DT>::
rangeCtorHelper(TreeNode*& ndPoolPtr, size_t dim, RAI begin,
                RAI median, RAI end) {
    std::nth_element(begin, median, end, [=](const auto& nh1, const auto& nh2) {
        return nh1.key[dim] < nh2.key[dim];});
    pool->construct(ndPoolPtr, std::move(median->key),
                   std::move(median->value));
    auto curNdPtr = ndPoolPtr;
    size_t nextDim = dim == N - 1 ? 0 : dim + 1;
    
    if (begin == median - 1) {
        curNdPtr->left = ++ndPoolPtr;
        pool->construct(ndPoolPtr, std::move(begin->key),
                        std::move(begin->value));
    } else if (begin != median) {
        curNdPtr->left = ++ndPoolPtr;
        rangeCtorHelper(ndPoolPtr, nextDim, begin,
                        begin + (median-begin)/2, median);
    }
    
    if (median + 2 == end) {
        curNdPtr->right = ++ndPoolPtr;
        pool->construct(ndPoolPtr, std::move((median + 1)->key),
                        std::move((median+1)->value));
    } else if (median + 1 != end) {
        curNdPtr->right = ++ndPoolPtr;
        rangeCtorHelper(ndPoolPtr, nextDim, median+1,
                        median + (end-median+1)/2, end);
    }
}

 
template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeCusMem<_Tp, N, ElemType, DT>::treeCopy(TreeNode*& thisNode,
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

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
KDTreeCusMem<_Tp, N, ElemType, DT>::~KDTreeCusMem() {
    if (pool)
        pool->destroy_and_free_all<TreeNode>();
}

// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
constexpr size_t KDTreeCusMem<_Tp, N, ElemType, DT>::dimension() const {
    return N;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
typename Point<_Tp, N>::DistType KDTreeCusMem<_Tp, N, ElemType, DT>::distType() const {
    return DT;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
size_t KDTreeCusMem<_Tp, N, ElemType, DT>::size() const {
    return treeSize;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
int KDTreeCusMem<_Tp, N, ElemType, DT>::height() const {
    return heightHelper(root);
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
int KDTreeCusMem<_Tp, N, ElemType, DT>::heightHelper(TreeNode *n) const {
    return n ? 1 + std::max(heightHelper(n->left),
                            heightHelper(n->right)) : -1;
}


template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
bool KDTreeCusMem<_Tp, N, ElemType, DT>::empty() const {
    return treeSize == 0;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeCusMem<_Tp, N, ElemType, DT>::printTreeInfo() const {
    std::cout << "Tree height is " << height()
              << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeCusMem<_Tp, N, ElemType, DT>::clear() {
    pool->destroy_and_free_all<TreeNode>();
    root = nullptr;
    treeSize = 0;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeCusMem<_Tp, N, ElemType, DT>::
insert(const Point<_Tp, N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        (*ndPtr)->object = value;
    } else {
        *ndPtr = new TreeNode(pt, value);
        treeSize++;
    }
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
bool KDTreeCusMem<_Tp, N, ElemType, DT>::contains(const Point<_Tp, N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
ElemType& KDTreeCusMem<_Tp, N, ElemType, DT>::operator[] (const Point<_Tp, N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        treeSize++;
    }
    return (*ndPtr)->object;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
ElemType& KDTreeCusMem<_Tp, N, ElemType, DT>::at(const Point<_Tp, N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeCusMem&>(*this).at(pt));
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
const ElemType& KDTreeCusMem<_Tp, N, ElemType, DT>::at(const Point<_Tp, N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
typename KDTreeCusMem<_Tp, N, ElemType, DT>::TreeNode**
KDTreeCusMem<_Tp, N, ElemType, DT>::findNodePtr(const Point<_Tp, N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTreeCusMem*>(this)->findNodePtr(pt));
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
typename KDTreeCusMem<_Tp, N, ElemType, DT>::TreeNode*const*
KDTreeCusMem<_Tp, N, ElemType, DT>::findNodePtr(const Point<_Tp, N>& pt) const {
    TreeNode *const*n = &root;
    for (size_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1)
        n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
ElemType KDTreeCusMem<_Tp, N, ElemType, DT>::
kNNValue(const Point<_Tp, N>& pt, size_t k) const {
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType, _Tp> bpq(k);
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

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
void KDTreeCusMem<_Tp, N, ElemType, DT>::kNNValueHelper(TreeNode *cur, size_t dim,
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
template <class OutputIter>
OutputIter KDTreeCusMem<_Tp, N, ElemType, DT>::rangeDiffKNNPairs(const Point<_Tp, N>& pt,
                                                _Tp fence, OutputIter returnIt) const {
    /*
     std::vector<std::pair<_Tp, std::pair<Point<_Tp, N>, ElemType>>> distKVPairs;
     distKVPairs.reserve(sqrt(treeSize));
     _Tp bestSq = std::numeric_limits<_Tp>::max(),
     bestDiffSq = std::numeric_limits<_Tp>::max();
     rangeDiffKNNPairsHelper(root, 0, key, diff, distKVPairs, bestSq, bestDiffSq);
     for (const auto &p : distKVPairs) {
     if (p.first < bestDiffSq)
     *it++ = std::move(p.second);
     }
     return it;
     */
    
    std::vector<std::pair<_Tp, std::pair<const Point<_Tp, N>*, const ElemType*>>> distKVPairs;
    distKVPairs.reserve(sqrt(treeSize));
    _Tp bestDistSq = std::numeric_limits<_Tp>::max(),
    bestDistDiffSq = std::numeric_limits<_Tp>::max(),
    curDistSq, fenceSq = fence*fence;
    size_t dim = 0;
    TreeNode *cur = root, *next;
    bool hasNext = true;
    
    
    
    struct ActRecord {
        _Tp curDist;
        TreeNode *cur;
        size_t dim;
    };
    ActRecord st[static_cast<size_t>(log2(treeSize+1))],
    *it = st;
    while (it != st || hasNext) {
        if (!hasNext) {
            const auto &ar = *--it;
            if (!(ar.curDist < bestDistDiffSq && (cur = ar.cur)))
                continue;
            dim = ar.dim;
            hasNext = true;
        }
        curDistSq = Point<_Tp, N>::template
        dist<Point<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDistSq < bestDistDiffSq) {
            if (curDistSq < bestDistSq) {
                bestDistSq = curDistSq;
                bestDistDiffSq = bestDistSq + fenceSq + 2*fence*sqrt(bestDistSq);
            }
            distKVPairs.emplace_back(curDistSq,
                                     std::forward<std::pair<const Point<_Tp, N>*,
                                     const ElemType*>>({&cur->key, &cur->object}));
        }
        curDistSq = pt[dim] - cur->key[dim];
        next = curDistSq < 0.0 ? cur->left : cur->right;
        curDistSq *= curDistSq;
        dim = dim == N - 1 ? 0 : dim + 1;
        if (next) {
            *it++ = {curDistSq, next == cur->left ? cur->right : cur->left, dim};
            cur = next;
        } else {
            if (curDistSq < bestDistDiffSq) {
                next = next == cur->left ? cur->right : cur->left;
                if (next) {
                    cur = next;
                    continue;
                }
            }
            hasNext = false;
        }
    }
    
    for (const auto& [dist, kvPairs] : distKVPairs) {
        if (dist < bestDistDiffSq)
            *returnIt++ = {*kvPairs.first, *kvPairs.second};
    }
    return returnIt;
    
}

/*
 template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
 void KDTreeCusMem<_Tp, N, ElemType, DT>::
 rangeDiffKNNPairsHelper(TreeNode *cur, size_t dim, const Point<_Tp, N>& pt,
 _Tp diff, std::vector<std::pair<_Tp,
 std::pair<Point<_Tp, N>, ElemType>>> &distKVPairs,
 _Tp& bestDistSq, _Tp&bestDistDiffSq) const {
 auto distSq = Point<_Tp, N>::template dist<Point<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
 if (distSq < bestDistDiffSq) {
 if (distSq < bestDistSq) {
 bestDistSq = distSq;
 bestDistDiffSq = bestDistSq + diff*diff + 2*diff*sqrt(bestDistSq);
 }
 distKVPairs.emplace_back(distSq, std::make_pair(cur->key, cur->object));
 }
 size_t nextDim = dim + 1 < N ? dim + 1 : 0;
 _Tp thisDiff = pt[dim] - cur->key[dim];
 TreeNode *next = thisDiff < 0 ? cur->left : cur->right;
 if (next)
 rangeDiffKNNPairsHelper(next, nextDim, pt, diff, distKVPairs, bestDistSq, bestDistDiffSq);
 if (thisDiff*thisDiff < bestDistDiffSq) {
 next = next == cur->left ? cur->right : cur->left;
 if (next)
 rangeDiffKNNPairsHelper(next, nextDim, pt, diff, distKVPairs, bestDistSq, bestDistDiffSq);
 }
 }*/

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
ElemType KDTreeCusMem<_Tp, N, ElemType, DT>::NNValue(const Point<_Tp, N> &pt) const {
    
    
    _Tp bestDist = std::numeric_limits<_Tp>::max(), curDist, diff;
    size_t dim = 0, nextDim;
    const ElemType *bestValue = nullptr;
    TreeNode *cur = root, *next;
    bool hasNext = true;
    
    struct actRecord {
        _Tp curDist;
        TreeNode *cur;
        size_t dim;
    };
    actRecord st[static_cast<size_t>(log2(treeSize+1))], *it = st;
    
    
    // LOGGGGGGGGGGGGGGG
    
    
    static size_t totalNumNodesSearches = 0, numNNSearches = 0, totalTreeSize = 0;
    static size_t numOfFullSearch = 0;
    size_t thisNumNodesSearches = 0;
    static bool logCondition;
    static constexpr size_t TREE_SIZE_LOWER_BOUND = 1200, TREE_SIZE_UPPER_BOUND = 600000000;
    logCondition = treeSize <= TREE_SIZE_UPPER_BOUND && treeSize >= TREE_SIZE_LOWER_BOUND;
    


    while (it != st || hasNext) {
        if (!hasNext) {
            const auto &ar = *--it;
            if (!(ar.curDist < bestDist && (cur = ar.cur)))
                continue;
            dim = ar.dim;
            hasNext = true;
        }
        nextDim = dim == N - 1 ? 0 : dim + 1;
        curDist = Point<_Tp, N>::template dist<Point<_Tp, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDist < bestDist) {
            bestDist = curDist;
            bestValue = &cur->object;
        }
        
        // LOGGGGGGGGGGGGGGG
        if (logCondition)
            thisNumNodesSearches++;
        
        
        diff = pt[dim] - cur->key[dim];
        next = diff < 0.0 ? cur->left : cur->right;
        curDist = diff*diff;
        if (next) {
            *it++ = {curDist, next == cur->left ? cur->right : cur->left, nextDim};
            cur = next;
            dim = nextDim;
        } else {
            if (curDist < bestDist) {
                next = next == cur->left ? cur->right : cur->left;
                if (next) {
                    cur = next;
                    dim = nextDim;
                    continue;
                }
            }
            hasNext = false;
        }
    }
    
    // LOGGGGGGGGGGGGGGG
    
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
    
    
    /*
     
     _Tp bestDist = std::numeric_limits<_Tp>::max();
     const ElemType *bestValue = nullptr;
     NNValueHelper(root, 0, pt, bestValue, bestDist);
     */
    
    return bestValue ? *bestValue : ElemType();
    
}

template <typename _Tp, size_t N, typename ElemType, typename Point<_Tp, N>::DistType DT>
template <typename Point<_Tp, N>::DistType thisDt,
typename std::enable_if<thisDt == Point<_Tp, N>::DistType::EUC, int>::type>
void KDTreeCusMem<_Tp, N, ElemType, DT>::
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
void KDTreeCusMem<_Tp, N, ElemType, DT>::
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
_Tp KDTreeCusMem<_Tp, N, ElemType, DT>::branchMin(const Point<_Tp, N> &trPt,
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

#endif /* KDTreeCusMem_hpp */
