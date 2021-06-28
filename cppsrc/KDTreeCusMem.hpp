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


template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT
= PointND<FPType, N>::DistType::EUC>
class KDTreeCusMem {
public:
    
    typedef FPType                                   value_type;
    
    struct node_type {
        PointND<value_type, N> key;
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
    
    template <typename ConstRAI,
        typename std::enable_if<std::is_same<
        typename std::iterator_traits<typename
        std::remove_const_t<ConstRAI>>::iterator_category,
        std::random_access_iterator_tag>::value && std::is_const<typename
        std::remove_pointer<typename std::iterator_traits<ConstRAI>::pointer>
           ::type>::value, int>::type = 0>
    KDTreeCusMem(ConstRAI, ConstRAI);

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
    
    // std::uint8_t dimension() const;
    // Usage: std::uint8_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTreeCusMem.
    constexpr std::uint8_t dimension() const;
    typename PointND<value_type, N>::DistType distType() const;
    
    // std::size_t size() const;
    // std::size_t height() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree, the max height
    // and whether the tree is empty.
    std::size_t size() const;
    int height() const;
    bool empty() const;
    
    void Clear();
    
    void PrintTreeInfo() const;
    
    // bool contains(const PointND<FPType, N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTreeCusMem.
    bool contains(const PointND<value_type, N>&) const;
    
    // void insert(const PointND<FPType, N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTreeCusMem, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const PointND<value_type, N>&, const ElemType&);
    
    // ElemType& operator[](const PointND<FPType, N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTreeCusMem.
    // If the point does not exist, then it is added to the KDTreeCusMem using the
    // default value of ElemType as its key.
    ElemType& operator[](const PointND<value_type, N>& pt);
    
    // ElemType& at(const PointND<FPType, N>& pt);
    // const ElemType& at(const PointND<FPType, N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function //throws an out_of_range exception.
    ElemType& at(const PointND<value_type, N>& pt);
    const ElemType& at(const PointND<value_type, N>& pt) const;
    
    // ElemType kNNValue(const PointND<FPType, N>& key, std::size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTreeCusMem
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const PointND<value_type, N>& key, std::size_t k) const;
    
    // Iter NNsWithFence(const PointND<FPType, N>&, FPType, Iter) const
    // Usage: Iter end = kd.NNsWithFence(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a FPType offset, return a set of points in the KDTreeCusMem
    // nearest to p such that the farthest one in the set is at least offset
    // distance close to p than the rest of the points in the tree.
    // The forward iterator is passed in and filled and the end will be returned.
    template <class OutputIter>
    OutputIter NNsWithFence(const PointND<value_type, N>&, value_type, OutputIter) const;
    
private:
    
    struct TreeNode {
        PointND<value_type, N> key;
        TreeNode *left;
        TreeNode *right;
        ElemType object;
        
        TreeNode() = default;
        // TreeNode(const TreeNode&) = default;
        // TreeNode& operator = (const TreeNode&) = default;
        // TreeNode(TreeNode&&) = default;
        // TreeNode& operator = (TreeNode&&) = default;
        
        TreeNode(const PointND<value_type, N>& k, const ElemType& obj)
        : key(k), left(nullptr), right(nullptr), object(obj) {}
        
        TreeNode(PointND<value_type, N>&& k, ElemType&& obj)
        : key(std::move(k)), left(nullptr), right(nullptr),
        object(std::move(obj)) {}
        
        ~TreeNode() = default;
    };
    
    
    TreeNode *root;
    std::size_t treeSize;
    std::unique_ptr<PooledAllocator> pool;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    int heightHelper(TreeNode *n) const;
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void RangeCtorHelper(TreeNode*&, std::uint8_t dim, RAI, RAI, RAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(TreeNode *cur, std::uint8_t dim, const PointND<value_type, N> &pt,
                        BoundedPQueue<ElemType, value_type> &bpq) const;
    
    //void NNsWithFenceHelper(TreeNode*, std::size_t, const PointND<FPType, N>&, FPType,
    //                            std::vector<std::pair<FPType,
    //                         std::pair<PointND<FPType, N>, ElemType>>>&, FPType&, FPType&) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const PointND<value_type, N>& key) const;
    
    template <typename PointND<value_type, N>::DistType thisDt = DT,
              typename std::enable_if<thisDt == PointND<value_type, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, std::uint8_t dim, const PointND<value_type, N>&,
                       const ElemType *&, value_type&) const;
    
    template <typename PointND<value_type, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt != PointND<value_type, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, std::uint8_t dim, const PointND<value_type, N>&,
                       const ElemType*&, value_type&) const;
    
    // ----------------------------------------------------
    // Helper meothod for deep copy
    void treeCopy(TreeNode*& thisNd, TreeNode *otherNd //TreeNode* ndPoolPtr
    );
    
    // TreeNode** findNodePtr(const PointND<FPType, N>& pt);
    // TreeNode*const* findNodePtr(const PointND<FPType, N>& pt) const;
    // Usage: TreeNode **nodePtr = findNodePtr(pt);
    // ----------------------------------------------------
    // Returns the pointer pointing to the node address
    // corresponding to the given PointND. In this FPType pointing
    // fashion, we can construct a node at that location.
    TreeNode** findNodePtr(const PointND<value_type, N>& pt);
    TreeNode*const* findNodePtr(const PointND<value_type, N>& pt) const;
    
    value_type branchMin(const PointND<value_type, N>&, const PointND<value_type, N>&, std::size_t) const;
    
};

/** KDTreeCusMem class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------


template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename ConstRAI,
    typename std::enable_if<std::is_same<typename std::iterator_traits<typename
    std::remove_const_t<ConstRAI>>::iterator_category,
    std::random_access_iterator_tag>::value &&
    std::is_const<typename std::remove_pointer<typename
    std::iterator_traits<ConstRAI>::pointer>::type>::value, int>::type>
KDTreeCusMem<FPType, N, ElemType, DT>::KDTreeCusMem(ConstRAI cbegin, ConstRAI cend)
: treeSize(cend-cbegin), pool(std::make_unique<PooledAllocator>()) {
    std::vector<node_type> constructData(cbegin, cend);
    TreeNode* ndPoolPtr = root = pool->allocateExact<TreeNode>(treeSize);
    RangeCtorHelper(ndPoolPtr, 0, constructData.begin(), constructData.end(),
                    constructData.begin() +
                    (constructData.end() - constructData.begin())/2);
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename RAI, typename std::enable_if<std::is_same<typename
    std::iterator_traits<RAI>::iterator_category,
    std::random_access_iterator_tag>::value && !std::is_const<typename
    std::remove_pointer< typename std::iterator_traits<RAI>::pointer>::type>::value, int>::type>
KDTreeCusMem<FPType, N, ElemType, DT>::KDTreeCusMem(RAI begin, RAI end)
                                            //, std::array<FPType, N> bboxHint)
: treeSize(end-begin), pool(std::make_unique<PooledAllocator>()) {
    TreeNode* ndPoolPtr = root = pool->allocateExact<TreeNode>(treeSize);
    RangeCtorHelper(ndPoolPtr, 0, begin, begin + (end - begin)/2, end);
    
    /*
    struct actRecord {
        TreeNode** cur_ndPtr;
        std::uint8_t dim;
        RAI thisBeginIt, median, thisEndIt;
    };
    
    actRecord st[static_cast<std::size_t>(log2(treeSize))], *it = st;
    bool hasChild = true;
    TreeNode* cur_nd;
    std::uint8_t dim = 0;
    RAI thisBeginIt = begin, thisEndIt = end, median = thisBeginIt + (thisEndIt-thisBeginIt)/2;
    while (it != st || hasChild) {
        if (!hasChild) {
            hasChild = true;
            *(--it)->cur_ndPtr = ndPoolPtr;
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
        cur_nd = ndPoolPtr++;
        dim = dim == N - 1 ? 0 : dim + 1;
        
        if (thisBeginIt != median) {
            if (median + 1 != thisEndIt) {
                *it++ = {&cur_nd->right, dim, median + 1,
                         median + (thisEndIt-median+1)/2, thisEndIt};
            }
            thisEndIt = median;
            median = thisBeginIt + (median-thisBeginIt)/2;
            cur_nd->left = ndPoolPtr;
        } else if (median + 1 != thisEndIt) {
            thisBeginIt = median+1;
            median = median + (thisEndIt-median+1)/2;
            cur_nd->right = ndPoolPtr;
        } else {
            hasChild = false;
        }
    }
   */
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeCusMem<FPType, N, ElemType, DT>::KDTreeCusMem(const KDTreeCusMem& rhs)
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

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeCusMem<FPType, N, ElemType, DT>&
KDTreeCusMem<FPType, N, ElemType, DT>::operator=(const KDTreeCusMem& rhs) & {
    if (this != &rhs) {
        delete root;
        treeSize = rhs.treeSize;
        root = new TreeNode();
        treeCopy(root, rhs.root);
    }
    return *this;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeCusMem<FPType, N, ElemType, DT>::KDTreeCusMem(KDTreeCusMem&& rhs) noexcept
: root(rhs.root), treeSize(rhs.treeSize), pool(std::move(rhs.pool)) {
    rhs.root = nullptr;
    rhs.treeSize = 0;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeCusMem<FPType, N, ElemType, DT>& KDTreeCusMem<FPType, N, ElemType, DT>::
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
template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <class RAI>
void KDTreeCusMem<FPType, N, ElemType, DT>::
RangeCtorHelper(TreeNode*& ndPoolPtr, RAI begin,
                RAI median, RAI end, std::array<FPType, N>& bbox) {
    std::uint8_t dim = std::max_element(bbox.begin(), bbox.end()) - bbox.begin();
    std::nth_element(begin, median, end, [=](const auto& p1, const auto& p2) {
        return p1.first[dim] < p2.first[dim];});
    pool->construct(ndPoolPtr, std::move(median->first),
                    std::move(median->second));
    auto cur_ndPtr = ndPoolPtr;
    //std::uint8_t next_dim = dim == N - 1 ? 0 : dim + 1;
    FPType dec;
    
    if (begin == median - 1) {
        cur_ndPtr->left = ++ndPoolPtr;
        pool->construct(ndPoolPtr, std::move(begin->first),
                        std::move(begin->second));
    } else if (begin != median) {
        dec = (end-1)->first[dim] - (median-1)->first[dim];
        bbox[dim] -= dec;
        //next_dim = std::max_element(bbox, bbox+N)-bbox;
        cur_ndPtr->left = ++ndPoolPtr;
        RangeCtorHelper(ndPoolPtr, begin,
                        begin + (median-begin)/2, median, bbox);
        bbox[dim] += dec;
    }
    
    if (median + 2 == end) {
        cur_ndPtr->right = ++ndPoolPtr;
        pool->construct(ndPoolPtr, std::move((median + 1)->first),
                        std::move((median+1)->second));
    } else if (median + 1 != end) {
        dec = (median+1)->first[dim] - begin->first[dim];
        bbox[dim] -= dec;
        cur_ndPtr->right = ++ndPoolPtr;
        RangeCtorHelper(ndPoolPtr, median+1,
                        median + (end-median+1)/2, end, bbox);
        bbox[dim] += dec;
    }
}
*/

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <class RAI>
void KDTreeCusMem<FPType, N, ElemType, DT>::
RangeCtorHelper(TreeNode*& ndPoolPtr, std::uint8_t dim, RAI begin,
                RAI median, RAI end) {
    std::nth_element(begin, median, end, [=](const auto& nh1, const auto& nh2) {
        return nh1.key[dim] < nh2.key[dim];});
    pool->construct(ndPoolPtr, std::move(median->key),
                   std::move(median->value));
    auto cur_ndPtr = ndPoolPtr;
    std::uint8_t next_dim = dim == N - 1 ? 0 : dim + 1;
    
    if (begin == median - 1) {
        cur_ndPtr->left = ++ndPoolPtr;
        pool->construct(ndPoolPtr, std::move(begin->key),
                        std::move(begin->value));
    } else if (begin != median) {
        cur_ndPtr->left = ++ndPoolPtr;
        RangeCtorHelper(ndPoolPtr, next_dim, begin,
                        begin + (median-begin)/2, median);
    }
    
    if (median + 2 == end) {
        cur_ndPtr->right = ++ndPoolPtr;
        pool->construct(ndPoolPtr, std::move((median + 1)->key),
                        std::move((median+1)->value));
    } else if (median + 1 != end) {
        cur_ndPtr->right = ++ndPoolPtr;
        RangeCtorHelper(ndPoolPtr, next_dim, median+1,
                        median + (end-median+1)/2, end);
    }
}

 
template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeCusMem<FPType, N, ElemType, DT>::treeCopy(TreeNode*& thisNode,
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

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTreeCusMem<FPType, N, ElemType, DT>::~KDTreeCusMem() {
    if (pool)
        pool->destroy_and_free_all<TreeNode>();
}

// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
constexpr std::uint8_t KDTreeCusMem<FPType, N, ElemType, DT>::dimension() const {
    return N;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename PointND<FPType, N>::DistType KDTreeCusMem<FPType, N, ElemType, DT>::distType() const {
    return DT;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
std::size_t KDTreeCusMem<FPType, N, ElemType, DT>::size() const {
    return treeSize;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
int KDTreeCusMem<FPType, N, ElemType, DT>::height() const {
    return heightHelper(root);
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
int KDTreeCusMem<FPType, N, ElemType, DT>::heightHelper(TreeNode *n) const {
    return n ? 1 + std::max(heightHelper(n->left),
                            heightHelper(n->right)) : -1;
}


template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeCusMem<FPType, N, ElemType, DT>::empty() const {
    return treeSize == 0;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeCusMem<FPType, N, ElemType, DT>::PrintTreeInfo() const {
    std::cout << "Tree height is " << height()
              << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeCusMem<FPType, N, ElemType, DT>::Clear() {
    pool->destroy_and_free_all<TreeNode>();
    root = nullptr;
    treeSize = 0;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeCusMem<FPType, N, ElemType, DT>::
insert(const PointND<FPType, N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        (*ndPtr)->object = value;
    } else {
        *ndPtr = new TreeNode(pt, value);
        treeSize++;
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTreeCusMem<FPType, N, ElemType, DT>::contains(const PointND<FPType, N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType& KDTreeCusMem<FPType, N, ElemType, DT>::operator[] (const PointND<FPType, N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        treeSize++;
    }
    return (*ndPtr)->object;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType& KDTreeCusMem<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeCusMem&>(*this).at(pt));
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
const ElemType& KDTreeCusMem<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename KDTreeCusMem<FPType, N, ElemType, DT>::TreeNode**
KDTreeCusMem<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTreeCusMem*>(this)->findNodePtr(pt));
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename KDTreeCusMem<FPType, N, ElemType, DT>::TreeNode*const*
KDTreeCusMem<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) const {
    TreeNode *const*n = &root;
    for (std::uint8_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1)
        n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType KDTreeCusMem<FPType, N, ElemType, DT>::
kNNValue(const PointND<FPType, N>& pt, std::size_t k) const {
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType, FPType> bpq(k);
    kNNValueHelper(root, 0, pt, bpq);
    
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

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTreeCusMem<FPType, N, ElemType, DT>::kNNValueHelper(TreeNode *cur, std::uint8_t dim,
                                             const PointND<FPType, N>& pt, BoundedPQueue<ElemType, FPType> &bpq) const {
    bpq.enqueue(cur->object, PointND<FPType, N>::template dist<DT>(cur->key, pt));
    std::uint8_t next_dim = dim + 1 < N ? dim + 1 : 0;
    TreeNode *next = pt[dim] < cur->key[dim] ? cur->left : cur->right;
    if (next)
        kNNValueHelper(next, next_dim, pt, bpq);
    if (bpq.size() < bpq.maxSize()
        || branchMin(cur->key, pt, dim) < bpq.worst()) {
        TreeNode *other = next == cur->left ? cur->right : cur->left;
        if (other)
            kNNValueHelper(other, next_dim, pt, bpq);
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <class OutputIter>
OutputIter KDTreeCusMem<FPType, N, ElemType, DT>::NNsWithFence(const PointND<FPType, N>& pt,
                                                FPType fence, OutputIter returnIt) const {
    /*
     std::vector<std::pair<FPType, std::pair<PointND<FPType, N>, ElemType>>> distKVPairs;
     distKVPairs.reserve(sqrt(treeSize));
     FPType bestSq = std::numeric_limits<FPType>::max(),
     bestDiffSq = std::numeric_limits<FPType>::max();
     NNsWithFenceHelper(root, 0, key, diff, distKVPairs, bestSq, bestDiffSq);
     for (const auto &p : distKVPairs) {
     if (p.first < bestDiffSq)
     *it++ = std::move(p.second);
     }
     return it;
     */
    
    std::vector<std::pair<FPType, std::pair<const PointND<FPType, N>*, const ElemType*>>> distKVPairs;
    distKVPairs.reserve(sqrt(treeSize));
    FPType bestDistSq = std::numeric_limits<FPType>::max(),
    bestDistDiffSq = std::numeric_limits<FPType>::max(),
    curDistSq, fenceSq = fence*fence;
    std::uint8_t dim = 0;
    TreeNode *cur = root, *next;
    bool hasNext = true;
    
    
    
    struct ActRecord {
        FPType curDist;
        TreeNode *cur;
        std::uint8_t dim;
    };
    ActRecord st[32],
    *it = st;
    while (it != st || hasNext) {
        if (!hasNext) {
            const auto &ar = *--it;
            if (!(ar.curDist < bestDistDiffSq && (cur = ar.cur)))
                continue;
            dim = ar.dim;
            hasNext = true;
        }
        curDistSq = PointND<FPType, N>::template
        dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDistSq < bestDistDiffSq) {
            if (curDistSq < bestDistSq) {
                bestDistSq = curDistSq;
                bestDistDiffSq = bestDistSq + fenceSq + 2*fence*sqrt(bestDistSq);
            }
            distKVPairs.emplace_back(curDistSq,
                                     std::forward<std::pair<const PointND<FPType, N>*,
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
 template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
 void KDTreeCusMem<FPType, N, ElemType, DT>::
 NNsWithFenceHelper(TreeNode *cur, std::uint8_t dim, const PointND<FPType, N>& pt,
 FPType diff, std::vector<std::pair<FPType,
 std::pair<PointND<FPType, N>, ElemType>>> &distKVPairs,
 FPType& bestDistSq, FPType&bestDistDiffSq) const {
 auto distSq = PointND<FPType, N>::template dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
 if (distSq < bestDistDiffSq) {
 if (distSq < bestDistSq) {
 bestDistSq = distSq;
 bestDistDiffSq = bestDistSq + diff*diff + 2*diff*sqrt(bestDistSq);
 }
 distKVPairs.emplace_back(distSq, std::make_pair(cur->key, cur->object));
 }
 std::uint8_t next_dim = dim + 1 < N ? dim + 1 : 0;
 FPType thisDiff = pt[dim] - cur->key[dim];
 TreeNode *next = thisDiff < 0 ? cur->left : cur->right;
 if (next)
 NNsWithFenceHelper(next, next_dim, pt, diff, distKVPairs, bestDistSq, bestDistDiffSq);
 if (thisDiff*thisDiff < bestDistDiffSq) {
 next = next == cur->left ? cur->right : cur->left;
 if (next)
 NNsWithFenceHelper(next, next_dim, pt, diff, distKVPairs, bestDistSq, bestDistDiffSq);
 }
 }*/

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType KDTreeCusMem<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N> &pt) const {
    
    
    FPType bestDist = std::numeric_limits<FPType>::max(), curDist, diff;
    std::uint8_t dim = 0, next_dim;
    const ElemType *bestValue = nullptr;
    TreeNode *cur = root, *next;
    bool hasNext = true;
    
    struct actRecord {
        FPType curDist;
        TreeNode *cur;
        std::uint8_t dim;
    };
    actRecord st[32], *it = st;
    
    
    // LOGGGGGGGGGGGGGGG
    
    
    static std::size_t totalNumNodesSearches = 0, numNNSearches = 0, totalTreeSize = 0;
    static std::size_t numOfFullSearch = 0;
    std::size_t thisNumNodesSearches = 0;
    static bool logCondition;
    static constexpr std::size_t TREE_SIZE_LOWER_BOUND = 1200, TREE_SIZE_UPPER_BOUND = 600000000;
    logCondition = treeSize <= TREE_SIZE_UPPER_BOUND && treeSize >= TREE_SIZE_LOWER_BOUND;
    


    while (it != st || hasNext) {
        if (!hasNext) {
            const auto &ar = *--it;
            if (!(ar.curDist < bestDist && (cur = ar.cur)))
                continue;
            dim = ar.dim;
            hasNext = true;
        }
        next_dim = dim == N - 1 ? 0 : dim + 1;
        curDist = PointND<FPType, N>::template dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
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
            *it++ = {curDist, next == cur->left ? cur->right : cur->left, next_dim};
            cur = next;
            dim = next_dim;
        } else {
            if (curDist < bestDist) {
                next = next == cur->left ? cur->right : cur->left;
                if (next) {
                    cur = next;
                    dim = next_dim;
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
     
     FPType bestDist = std::numeric_limits<FPType>::max();
     const ElemType *bestValue = nullptr;
     NNValueHelper(root, 0, pt, bestValue, bestDist);
     */
    
    return bestValue ? *bestValue : ElemType();
    
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename PointND<FPType, N>::DistType thisDt,
typename std::enable_if<thisDt == PointND<FPType, N>::DistType::EUC, int>::type>
void KDTreeCusMem<FPType, N, ElemType, DT>::
NNValueHelper(TreeNode *cur, std::uint8_t dim, const PointND<FPType, N> &pt,
              const ElemType *&bestValue, FPType &bestDist) const {
    FPType curDist = PointND<FPType, N>::template
    dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    std::uint8_t next_dim = dim == N - 1 ? 0 : dim + 1;
    FPType diff = pt[dim] - cur->key[dim];
    TreeNode *next = diff < 0.0 ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    if (diff*diff < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename PointND<FPType, N>::DistType thisDt,
typename std::enable_if<thisDt != PointND<FPType, N>::DistType::EUC, int>::type>
void KDTreeCusMem<FPType, N, ElemType, DT>::
NNValueHelper(TreeNode *cur, std::uint8_t dim, const PointND<FPType, N> &pt,
              const ElemType *&bestValue, FPType &bestDist) const {
    FPType curDist = PointND<FPType, N>::template dist<DT>(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    std::uint8_t next_dim = dim == N - 1 ? 0 : dim + 1;
    TreeNode *next = pt[dim] < cur->key[dim] ? cur->left : cur->right;
    if (next)
        NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    if (branchMin(cur->key, pt, dim) < bestDist) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNValueHelper(next, next_dim, pt, bestValue, bestDist);
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
FPType KDTreeCusMem<FPType, N, ElemType, DT>::branchMin(const PointND<FPType, N> &trPt,
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

#endif /* KDTreeCusMem_hpp */
