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


template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT
= def::DistType::kEuc>
class KDTreeExpandLongest {
public:
    
    typedef FPType                                   value_type;
    
    struct node_type {
        PointND<value_type, N> key;
        ElemType value;
    };
    
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
    
    template <typename ConstRAI,
    typename std::enable_if<std::is_same<
    typename std::iterator_traits<typename
    std::remove_const_t<ConstRAI>>::iterator_category,
    std::random_access_iterator_tag>::value && std::is_const<typename
    std::remove_pointer<typename std::iterator_traits<ConstRAI>::pointer>
    ::type>::value, int>::type = 0>
    KDTreeExpandLongest(ConstRAI, ConstRAI);
    
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
    
    // std::uint8_t dimension() const;
    // Usage: std::uint8_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTreeExpandLongest.
    constexpr std::uint8_t dimension() const;
    typename def::DistType distType() const;
    
    // std::size_t size() const;
    // std::size_t height() const;
    // bool empty() const;
    // Usage: if (kd.empty())
    // ----------------------------------------------------
    // Returns the number of elements in the kd-tree, the max height
    // and whether the tree is empty.
    std::size_t size() const;
    std::uint32_t height() const;
    bool empty() const;
    
    void Clear();
    
    void PrintTreeInfo() const;
    
    // bool contains(const PointND<FPType, N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTreeExpandLongest.
    bool contains(const PointND<value_type, N>&) const;
    
    // void insert(const PointND<FPType, N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTreeExpandLongest, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const PointND<value_type, N>&, const ElemType&);
    
    // ElemType& operator[](const PointND<FPType, N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTreeExpandLongest.
    // If the point does not exist, then it is added to the KDTreeExpandLongest using the
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
    // Given a point v and an integer k, finds the k points in the KDTreeExpandLongest
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const PointND<value_type, N>& key, std::size_t k) const;
    
    // Iter NNsWithFence(const PointND<FPType, N>&, FPType, Iter) const
    // Usage: Iter end = kd.NNsWithFence(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a FPType offset, return a set of points in the KDTreeExpandLongest
    // nearest to p such that the farthest one in the set is at least offset
    // distance close to p than the rest of the points in the tree.
    // The forward iterator is passed in and filled and the end will be returned.
    template <class OutputIter>
    OutputIter NNsWithFence(const PointND<value_type, N>&, value_type, OutputIter) const;
    
private:
    
    
    struct TreeNode {
        TreeNode *right;
        std::uint8_t dim_to_expand;
        TreeNode *left;
        PointND<value_type, N> key;
        ElemType object;
        
        enum class dimEnum : std::size_t {
            LEFT_CHILD = N,
            NO_CHILDREN
        };
        
        TreeNode() = default;
        // TreeNode(const TreeNode&) = default;
        // TreeNode& operator = (const TreeNode&) = default;
        // TreeNode(TreeNode&&) = default;
        // TreeNode& operator = (TreeNode&&) = default;
        
        TreeNode(std::uint8_t dim_to_expand, const PointND<value_type, N>& k, const ElemType& obj)
        : right(nullptr), dim_to_expand(dim_to_expand), left(nullptr), key(k),
          object(obj) {}
        
        TreeNode(std::uint8_t dim_to_expand, PointND<value_type, N>&& k, ElemType&& obj)
        : right(nullptr), dim_to_expand(dim_to_expand), left(nullptr),
        key(std::move(k)), object(std::move(obj)) {}
        
        ~TreeNode() = default;
    };
    
    
    TreeNode *root;
    std::size_t treeSize;
    std::unique_ptr<PooledAllocator> pool;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    std::uint32_t heightHelper(TreeNode *n) const;
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void RangeCtorHelper(TreeNode*&, RAI, RAI, std::array<value_type, N*2>&);
    
    
    template <class RAI>
    static std::array<value_type, N*2> computeInitBBox(RAI, RAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(TreeNode *cur, std::uint8_t dim, const PointND<value_type, N> &pt,
                        BoundedPQueue<ElemType, value_type> &bpq) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const PointND<value_type, N>& key) const;
    
    template <typename def::DistType thisDt = DT,
    typename std::enable_if<thisDt == def::DistType::kEuc, int>::type = 0>
    void NNValueHelper(TreeNode*, std::uint8_t dim, const PointND<value_type, N>&,
                       const ElemType *&, value_type&) const;
    
    template <typename def::DistType thisDt = DT,
    typename std::enable_if<thisDt != def::DistType::kEuc, int>::type = 0>
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

/** KDTreeExpandLongest class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------


template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <typename ConstRAI,
typename std::enable_if<std::is_same<typename std::iterator_traits<typename
std::remove_const_t<ConstRAI>>::iterator_category,
std::random_access_iterator_tag>::value &&
std::is_const<typename std::remove_pointer<typename
std::iterator_traits<ConstRAI>::pointer>::type>::value, int>::type>
KDTreeExpandLongest<FPType, N, ElemType, DT>::KDTreeExpandLongest(ConstRAI cbegin, ConstRAI cend)
: treeSize(cend-cbegin), pool(std::make_unique<PooledAllocator>()) {
    std::vector<node_type> constructData(cbegin, cend);
    TreeNode* ndPoolPtr = root = pool->allocateExact<TreeNode>(treeSize);
    auto bbox = computeInitBBox(cbegin, cend);
    RangeCtorHelper(ndPoolPtr, constructData.begin(), constructData.end(), bbox);
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <typename RAI, typename std::enable_if<std::is_same<typename
std::iterator_traits<RAI>::iterator_category,
std::random_access_iterator_tag>::value && !std::is_const<typename
std::remove_pointer< typename std::iterator_traits<RAI>::pointer>::type>::value, int>::type>
KDTreeExpandLongest<FPType, N, ElemType, DT>::KDTreeExpandLongest(RAI begin, RAI end)
//, std::array<FPType, N> bboxHint)
: treeSize(end-begin), pool(std::make_unique<PooledAllocator>()) {
    auto bbox = computeInitBBox(begin, end);
    TreeNode* ndPoolPtr = root = pool->allocateExact<TreeNode>(treeSize);
    RangeCtorHelper(ndPoolPtr, begin, end, bbox);
    
    /*
     struct actRecord {
         TreeNode** cur_ndPtr;
         RAI thisBeginIt, median, thisEndIt;
     };
     
     actRecord st[static_cast<std::size_t>(
     (treeSize))], *it = st;
     bool hasChild = true;
     RAI thisBeginIt = begin, thisEndIt = end, median = thisBeginIt + (thisEndIt-thisBeginIt)/2;
     while (it != st || hasChild) {
         if (!hasChild) {
             hasChild = true;
             *(--it)->cur_ndPtr = ndPoolPtr;
             thisBeginIt = it->thisBeginIt;
             thisEndIt = it->thisEndIt;
             median = it->median;
         }
     
         std::uint8_t dim = 0;
         FPType maxSpan = bbox[1] - bbox[0];
         for (std::size_t i = 1; i < N; ++i) {
             std::size_t bboxLowIdx = i * 2, bboxHighIdx = bboxLowIdx + 1;
             FPType span = bbox[bboxHighIdx] - bbox[bboxLowIdx];
             if (span > maxSpan) {
                 maxSpan = span;
                 dim = i;
             }
         }
         std::nth_element(thisBeginIt, median, thisEndIt,
                          [=](const auto& p1, const auto& p2) {
                              return p1.first[dim] < p2.first[dim];});
         pool->construct(ndPoolPtr, dim, std::move(median->first),
                         std::move(median->second));
         TreeNode* cur_nd = ndPoolPtr++;
     
         if (thisBeginIt != median) {
             if (median + 1 != thisEndIt) {
                 *it++ = {&cur_nd->right, median + 1,
                          median + (thisEndIt-median+1)/2, thisEndIt};
             }
             bbox[dim*2+1] = cur_nd->key[dim];
             thisEndIt = median;
             median = thisBeginIt + (median-thisBeginIt)/2;
             cur_nd->left = ndPoolPtr;
         } else if (median + 1 != thisEndIt) {
             bbox[dim*2] = cur_nd->key[dim];
             thisBeginIt = median+1;
             median += (thisEndIt-median+1)/2;
             cur_nd->right = ndPoolPtr;
         } else {
             hasChild = false;
         }
     } */
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongest<FPType, N, ElemType, DT>::KDTreeExpandLongest(const KDTreeExpandLongest& rhs)
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

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongest<FPType, N, ElemType, DT>&
KDTreeExpandLongest<FPType, N, ElemType, DT>::operator=(const KDTreeExpandLongest& rhs) & {
    if (this != &rhs) {
        delete root;
        treeSize = rhs.treeSize;
        root = new TreeNode();
        treeCopy(root, rhs.root);
    }
    return *this;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongest<FPType, N, ElemType, DT>::KDTreeExpandLongest(KDTreeExpandLongest&& rhs) noexcept
: root(rhs.root), treeSize(rhs.treeSize), pool(std::move(rhs.pool)) {
    rhs.root = nullptr;
    rhs.treeSize = 0;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongest<FPType, N, ElemType, DT>& KDTreeExpandLongest<FPType, N, ElemType, DT>::
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

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <class RAI>
std::array<FPType, N*2> KDTreeExpandLongest<FPType, N, ElemType, DT>::
computeInitBBox(RAI begin, RAI end) {
    std::array<FPType, N*2> bbox;
    std::generate(bbox.begin(), bbox.end(),
                  [lowHighToggle = 0lu, lowHigh = std::array<FPType, 2>{
                   std::numeric_limits<FPType>::max(),
                   std::numeric_limits<FPType>::min()}]() mutable {
                      return lowHigh[lowHighToggle++%2];});
    std::for_each(begin, end, [&](const auto &nh) mutable {
        for (std::size_t i = 0; i < N; ++i) {
            FPType ptValOnithDim = nh.key[i];
            auto &bboxLow = bbox[i*2], &bboxHigh = bbox[i*2+1];
            bboxLow = std::min(bboxLow, ptValOnithDim);
            bboxHigh = std::max(bboxHigh, ptValOnithDim);
        }
    });
    return bbox;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <class RAI>
void KDTreeExpandLongest<FPType, N, ElemType, DT>::
RangeCtorHelper(TreeNode*& ndPoolPtr, RAI begin, RAI end,
                std::array<FPType, N*2> &bbox) {
    std::uint8_t dim = 0;
    FPType maxSpan = bbox[1] - bbox[0];
    for (std::size_t i = 1; i != N; ++i) {
        std::size_t bboxLowIdx = i * 2, bboxHighIdx = bboxLowIdx + 1;
        FPType span = bbox[bboxHighIdx] - bbox[bboxLowIdx];
        if (span > maxSpan) {
            maxSpan = span;
            dim = i;
        }
    }
    
    RAI median = begin + (end - begin)/2;
    std::nth_element(begin, median, end, [=](const auto& nh1, const auto& nh2) {
        return nh1.key[dim] < nh2.key[dim];});
    pool->construct(ndPoolPtr, dim, std::move(median->key),
                    std::move(median->value));
    TreeNode* cur_ndPtr = ndPoolPtr;
    
    if (begin == median - 1) {
        //cur_ndPtr->dim_to_expand = static_cast<std::underlying_type_t<typename TreeNode::dimEnum>>(TreeNode::dimEnum::LEFT_CHILD);
        cur_ndPtr->left = ++ndPoolPtr;
        pool->construct(ndPoolPtr, static_cast<std::underlying_type_t<typename TreeNode::dimEnum>>(TreeNode::dimEnum::NO_CHILDREN), std::move(begin->key),
                        std::move(begin->value));
    } else if (begin != median) {
        cur_ndPtr->left = ++ndPoolPtr;
        auto prev_high_on_dim = bbox[dim*2+1];
        bbox[dim*2+1] = cur_ndPtr->key[dim];
        //bbox[dim*2+1] = std::max_element(begin, median, [dim](const auto &p1, const auto&p2){return p1.first[dim] < p2.first[dim];})->first[dim];
        RangeCtorHelper(ndPoolPtr, begin, median, bbox);
        bbox[dim*2+1] = prev_high_on_dim;
    }
    
    if (median + 2 == end) {
        //cur_ndPtr->dim_to_expand = dim;
        cur_ndPtr->right = ++ndPoolPtr;
        pool->construct(ndPoolPtr, static_cast<std::underlying_type_t<typename TreeNode::dimEnum>>(TreeNode::dimEnum::NO_CHILDREN), std::move((median+1)->key),
                        std::move((median+1)->value));
    } else if (median + 1 != end) {
        cur_ndPtr->right = ++ndPoolPtr;
        auto prev_low_on_dim = bbox[dim*2];
        bbox[dim*2] = cur_ndPtr->key[dim];
        //bbox[dim*2] = std::min_element(median + 1, end, [dim](const auto &p1, const auto &p2){return p1.first[dim] < p2.first[dim];})->first[dim];
        RangeCtorHelper(ndPoolPtr, median+1, end, bbox);
        bbox[dim*2] = prev_low_on_dim;
    }
}


template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongest<FPType, N, ElemType, DT>::treeCopy(TreeNode*& thisNode,
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
            thisNode = new TreeNode(otherNode->dim_to_expand, otherNode->key, otherNode->object);
        }
        treeCopy(thisNode->left, otherNode->left);
        treeCopy(thisNode->right, otherNode->right);
    } else if (thisNode) {
        delete thisNode;
        thisNode = nullptr;
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
KDTreeExpandLongest<FPType, N, ElemType, DT>::~KDTreeExpandLongest() {
    if (pool)
        pool->destroy_and_free_all<TreeNode>();
}

// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
constexpr std::uint8_t KDTreeExpandLongest<FPType, N, ElemType, DT>::dimension() const {
    return N;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
typename def::DistType KDTreeExpandLongest<FPType, N, ElemType, DT>::distType() const {
    return DT;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
std::size_t KDTreeExpandLongest<FPType, N, ElemType, DT>::size() const {
    return treeSize;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
std::uint32_t KDTreeExpandLongest<FPType, N, ElemType, DT>::height() const {
    return heightHelper(root);
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
std::uint32_t KDTreeExpandLongest<FPType, N, ElemType, DT>::heightHelper(TreeNode *n) const {
    return n ? 1 + std::max(heightHelper(n->left), heightHelper(n->right)) : -1;
}


template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
bool KDTreeExpandLongest<FPType, N, ElemType, DT>::empty() const {
    return treeSize == 0;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongest<FPType, N, ElemType, DT>::PrintTreeInfo() const {
    std::cout << "Tree height is " << height()
    << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongest<FPType, N, ElemType, DT>::Clear() {
    pool->destroy_and_free_all<TreeNode>();
    root = nullptr;
    treeSize = 0;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongest<FPType, N, ElemType, DT>::
insert(const PointND<FPType, N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        (*ndPtr)->object = value;
    } else {
        *ndPtr = new TreeNode(0, pt, value);
        treeSize++;
    }
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
bool KDTreeExpandLongest<FPType, N, ElemType, DT>::contains(const PointND<FPType, N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType& KDTreeExpandLongest<FPType, N, ElemType, DT>::operator[] (const PointND<FPType, N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        treeSize++;
    }
    return (*ndPtr)->object;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType& KDTreeExpandLongest<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeExpandLongest&>(*this).at(pt));
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
const ElemType& KDTreeExpandLongest<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
typename KDTreeExpandLongest<FPType, N, ElemType, DT>::TreeNode**
KDTreeExpandLongest<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTreeExpandLongest*>(this)->findNodePtr(pt));
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
typename KDTreeExpandLongest<FPType, N, ElemType, DT>::TreeNode*const*
KDTreeExpandLongest<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) const {
    TreeNode *const*n = &root;
    for (std::uint8_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1)
        n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongest<FPType, N, ElemType, DT>::
kNNValue(const PointND<FPType, N>& pt, std::size_t k) const {
    //if (empty())
    //    return ElemType();
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

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
void KDTreeExpandLongest<FPType, N, ElemType, DT>::kNNValueHelper(TreeNode *cur, std::uint8_t dim,
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

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <class OutputIter>
OutputIter KDTreeExpandLongest<FPType, N, ElemType, DT>::
NNsWithFence(const PointND<FPType, N>& pt, FPType fence, OutputIter returnIt) const {
    std::vector<std::tuple<FPType, const PointND<FPType, N>&, const ElemType&>> distPtElemTuple;
    distPtElemTuple.reserve(std::sqrt(treeSize));
    std::pair<FPType, const TreeNode*> st[32],
                                       *it = st;
    FPType bestDistSq = std::numeric_limits<FPType>::max(),
           bestDistDiffSq = std::numeric_limits<FPType>::max(),
           fenceSq = fence*fence;
    const TreeNode *cur = root;
   
    while (true) {
        FPType curDistSq = PointND<FPType, N>::template dist<def::DistType::kEucSq>(cur->key, pt);
        if (curDistSq < bestDistDiffSq) {
            if (curDistSq < bestDistSq) {
                bestDistSq = curDistSq;
                bestDistDiffSq = bestDistSq + fenceSq + 2*fence*sqrt(bestDistSq);
            }
            distPtElemTuple.emplace_back(curDistSq, cur->key, cur->object);
        }
        
        if (cur->right) {
            std::uint8_t dim = cur->dim_to_expand;
            FPType diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(it++) std::pair<FPType, const TreeNode*>(diff*diff, cur->right);
                cur = cur->left;
            } else {
                new(it++) std::pair<FPType, const TreeNode*>(diff*diff, cur->left);
                cur = cur->right;
            }
        } else if (!(cur = cur->left)) {
            do {
                if (it == st)
                    goto FINAL;
            } while ((--it)->first > bestDistDiffSq || !(cur = it->second));
        } 
    }
    
FINAL:
    for (const auto &[distSq, pt, elem] : distPtElemTuple) {
        if (distSq < bestDistDiffSq)
            *returnIt++ = {pt, elem};
    }
    return returnIt;
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
ElemType KDTreeExpandLongest<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N> &pt) const {
    
    std::pair<FPType, const TreeNode*> st[32],
                                       *it = st;
    FPType bestDist = std::numeric_limits<FPType>::max();
    const ElemType *bestValue = nullptr;
    const TreeNode *cur = root;
    
    // LOGGGGGGGGGGGGGGG
    
    /*
     static std::size_t totalNumNodesSearches = 0, numNNSearches = 0, totalTreeSize = 0;
     static std::uint8_t NumOfFullSearch = 0;
     std::size_t thisNumNodesSearches = 0;
     static bool logCondition;
     static constexpr std::size_t TREE_SIZE_LOWER_BOUND = 1200, TREE_SIZE_UPPER_BOUND = 60000000;
     logCondition = treeSize <= TREE_SIZE_UPPER_BOUND && treeSize >= TREE_SIZE_LOWER_BOUND;
    
    */
    
    while (true) {
        FPType curDist = PointND<FPType, N>::template dist<def::DistType::kEucSq>(cur->key, pt);
        if (curDist < bestDist) {
            bestDist = curDist;
            bestValue = &cur->object;
        }
        
        // LOGGGGGGGGGGGGGGG
         //if (logCondition)
        // thisNumNodesSearches++;
        
        if (cur->right) {
            std::uint8_t dim = cur->dim_to_expand;
            FPType diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(it++) std::pair<FPType, const TreeNode*>(diff*diff, cur->right);
                cur = cur->left;
            } else {
                new(it++) std::pair<FPType, const TreeNode*>(diff*diff, cur->left);
                cur = cur->right;
            }
        } else if (!(cur = cur->left)) {
            do {
                if (it == st)
                    return bestValue ? *bestValue : ElemType();
            } while ((--it)->first > bestDist || !(cur = it->second));
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
     FPType bestDist = std::numeric_limits<FPType>::max();
     const ElemType *bestValue = nullptr;
     NNValueHelper(root, 0, pt, bestValue, bestDist);
     */
    
    return bestValue ? *bestValue : ElemType();
    
}

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <typename def::DistType thisDt,
typename std::enable_if<thisDt == def::DistType::kEuc, int>::type>
void KDTreeExpandLongest<FPType, N, ElemType, DT>::
NNValueHelper(TreeNode *cur, std::uint8_t dim, const PointND<FPType, N> &pt,
              const ElemType *&bestValue, FPType &bestDist) const {
    FPType curDist = PointND<FPType, N>::template
    dist<def::DistType::kEucSq>(cur->key, pt);
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

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
template <typename def::DistType thisDt,
typename std::enable_if<thisDt != def::DistType::kEuc, int>::type>
void KDTreeExpandLongest<FPType, N, ElemType, DT>::
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

template <typename FPType, std::uint8_t N, typename ElemType, typename def::DistType DT>
FPType KDTreeExpandLongest<FPType, N, ElemType, DT>::branchMin(const PointND<FPType, N> &trPt,
                                                const PointND<FPType, N> &searchPt, std::size_t idx) const {
    switch (DT) {
        case def::DistType::kEuc:
        case def::DistType::kMan:
            return std::fabs(trPt[idx] - searchPt[idx]);
            /*
             case DistType::kHav:
             PointND<FPType, N> pt;
             pt[idx] = searchPt[idx];
             pt[1-idx] = trPt[1-idx];
             return PointND<FPType, N>::havDist(trPt, pt);
             */
    }
}





#endif /* KDTreeExpandLongest_hpp */
