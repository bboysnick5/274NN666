/**
 * File: KDTree.h
 * Author: Yunlong Nick Liu
 * ------------------------
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 */

#ifndef KDTREE_INCLUDED
#define KDTREE_INCLUDED

#include "Point.hpp"
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


template <typename FPType, std::size_t N, typename ElemType,
          typename PointND<FPType, N>::DistType DT = PointND<FPType, N>::DistType::EUC>
class KDTree {
public:
    
    typedef FPType                                   value_type;

    struct node_type {
        PointND<value_type, N> key;
        ElemType value;
    };
    
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    KDTree() = default;
    
    // Constructor: KDTree(FwdItType begin, FwdItType end);
    // Usage: KDTree<3, int> myTree(vec.begin(), bec.end());
    // ----------------------------------------------------
    // Constructs a KDTree from a collection. The tree will
    // be balanced using median constructing method
    // NOTE: The tree will not eliminate duplicates and the
    //       intended behavior will not be comprimised, tho
    //       less efficient with extra wasteful space.
    template <typename RAI, typename std::enable_if<std::is_same<typename
    std::iterator_traits<RAI>::iterator_category,
    std::random_access_iterator_tag>::value &&
    !std::is_const<typename std::remove_pointer<
    typename std::iterator_traits<RAI>::pointer>::type>::value, int>::type = 0>
    KDTree(RAI, RAI);
    
    template <typename ConstRAI,
    typename std::enable_if<std::is_same<typename std::iterator_traits<typename
    std::remove_const_t<ConstRAI>>::iterator_category,
    std::random_access_iterator_tag>::value &&
    std::is_const<typename std::remove_pointer<
    typename std::iterator_traits<ConstRAI>::pointer>::type>::value, int>::type = 0>
    KDTree(ConstRAI cbegin, ConstRAI cend);
    
    // Destructor: ~KDTree()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTree.
    ~KDTree();
    
    // KDTree(const KDTree& rhs);
    // KDTree& operator=(const KDTree& rhs);
    // Usage: KDTree<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Copy constructor and copy assignment operator.
    KDTree(const KDTree&);
    KDTree& operator=(const KDTree&)&;
    
    // KDTree(const KDTree& rhs);
    // KDTree& operator=(const KDTree& rhs);
    // Usage: KDTree<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Move constructor and move assignment operator.
    KDTree(KDTree&&) noexcept;
    KDTree& operator=(KDTree&&)& noexcept;
    
    // std::size_t dimension() const;
    // Usage: std::size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
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
    int height() const;
    bool empty() const;
    
    void Clear();
    
    void PrintTreeInfo() const;
    
    // bool contains(const PointND<FPType, N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTree.
    bool contains(const PointND<value_type, N>&) const;
    
    // void insert(const PointND<FPType, N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTree, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const PointND<value_type, N>&, const ElemType&);
    
    // ElemType& operator[](const PointND<FPType, N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTree.
    // If the point does not exist, then it is added to the KDTree using the
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
    // Given a point v and an integer k, finds the k points in the KDTree
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const PointND<value_type, N>& key, std::size_t k) const;
    
    // Iter NNsWithFence(const PointND<FPType, N>&, FPType, Iter) const
    // Usage: Iter end = kd.NNsWithFence(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a FPType offset, return a set of points in the KDTree
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
        :  key(std::move(k)), left(nullptr), right(nullptr),
          object(std::move(obj)) {}
        
        ~TreeNode() {
            delete left;
            delete right;
        }
    };
    
    
    TreeNode *root;
    std::size_t treeSize;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    int heightHelper(TreeNode *n) const;
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void RangeCtorHelper(TreeNode*&, std::size_t, RAI, RAI, RAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(TreeNode *cur, std::size_t dim, const PointND<value_type, N> &pt,
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
    void NNValueHelper(TreeNode*, std::size_t, const PointND<value_type, N>&,
                       const ElemType *&, value_type&) const;
    
    template <typename PointND<value_type, N>::DistType thisDt = DT,
    typename std::enable_if<thisDt != PointND<value_type, N>::DistType::EUC, int>::type = 0>
    void NNValueHelper(TreeNode*, std::size_t, const PointND<value_type, N>&,
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

/** KDTree class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------



template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename ConstRAI,
typename std::enable_if<std::is_same<typename std::iterator_traits<typename
std::remove_const_t<ConstRAI>>::iterator_category,
std::random_access_iterator_tag>::value && std::is_const<typename
std::remove_pointer< typename std::iterator_traits<ConstRAI>::pointer>::type>::value, int>::type>
KDTree<FPType, N, ElemType, DT>::KDTree(ConstRAI cbegin, ConstRAI cend)
: treeSize(cend-cbegin) {
    std::vector<node_type> constructData(cbegin, cend);
//std::copy(cbegin, cend, std::back_inserter(constructData));
    RangeCtorHelper(root, 0, constructData.begin(), constructData.end(),
                    constructData.begin() +
                    (constructData.end() - constructData.begin())/2);
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename RAI, typename std::enable_if<std::is_same<typename
std::iterator_traits<RAI>::iterator_category,
std::random_access_iterator_tag>::value &&
!std::is_const<typename std::remove_pointer<
typename std::iterator_traits<RAI>::pointer>::type>::value, int>::type>
KDTree<FPType, N, ElemType, DT>::KDTree(RAI begin, RAI end) : treeSize(end-begin) {
    
    if (treeSize == 1) {
        root = new TreeNode(std::move(begin->key),
                            std::move(begin->value));
    } else if (treeSize > 1) {
        RangeCtorHelper(root, 0, begin, begin + (end - begin)/2, end);
    } else {
        this->root = nullptr;
        treeSize = 0;
        std::cerr << "invalid iterators to construct KDTree.\n";
    }
     

    
    /*
    struct actRecord {
        TreeNode** cur_ndPtr;
        std::size_t dim;
        RAI thisBeginIt, median, thisEndIt;
    };
    
    actRecord st[static_cast<std::size_t>(log2(treeSize+1))];
    actRecord *it = st;
    bool hasChild = true;
    TreeNode** cur_ndPtr = &root;
    std::size_t dim = 0, nextDim = dim == N - 1 ? 0 : dim + 1;
    RAI thisBeginIt = begin, thisEndIt = end, median = thisBeginIt + (thisEndIt-thisBeginIt)/2;
    while (it != st || hasChild) {
        if (!hasChild) {
            hasChild = true;
            --it;
            cur_ndPtr = it->cur_ndPtr;
            dim = it->dim;
            thisBeginIt = it->thisBeginIt;
            thisEndIt = it->thisEndIt;
            median = it->median;
        }
        
        std::nth_element(thisBeginIt, median, thisEndIt,
                         [=](const std::pair<PointND<FPType, N>, ElemType>& p1,
                                   const std::pair<PointND<FPType, N>, ElemType>& p2) {
                                        return p1.first[dim] < p2.first[dim];});
        *cur_ndPtr = new TreeNode(std::move(median->first),
                                 std::move(median->second));
        nextDim = dim == N - 1 ? 0 : dim + 1;
        
       if (median != thisBeginIt) {
            if (median + 1 != thisEndIt) {
                *it++ = {&(*cur_ndPtr)->left, nextDim, thisBeginIt, thisBeginIt + (median-thisBeginIt)/2, median};
                cur_ndPtr = &(*cur_ndPtr)->right;
                thisBeginIt = median+1;
                median = median + (thisEndIt-median+1)/2;
                dim = nextDim;
            } else {
                cur_ndPtr = &(*cur_ndPtr)->left;
                thisEndIt = median;
                median = thisBeginIt + (median-thisBeginIt)/2;
                dim = nextDim;
            }
       } else if (median + 1 != thisEndIt) {
           cur_ndPtr = &(*cur_ndPtr)->right;
           thisBeginIt = median+1;
           median = median + (thisEndIt-median+1)/2;
           dim = nextDim;
       } else {
            hasChild = false;
       }
    } */
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTree<FPType, N, ElemType, DT>::KDTree(const KDTree& rhs)
: root(new TreeNode()), treeSize(rhs.treeSize) {
    // wrong logic.
    // should be check whether this size is greater than other.
    // if not allocate additional space. if yes
    //std::destroy_n(root, treeSize);
    //pool.free_all();
    //root = pool.allocate<TreeNode>(rhs.treeSize);
    //auto *ndPoolIt = addressof(*root);
    //pool.construct(root, {rhs.root->key, rhs.root->obj});
    treeCopy(root, rhs.root);
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTree<FPType, N, ElemType, DT>&
KDTree<FPType, N, ElemType, DT>::operator=(const KDTree& rhs) & {
    if (this != &rhs) {
        delete root;
        treeSize = rhs.treeSize;
        root = new TreeNode();
        treeCopy(root, rhs.root);
    }
    return *this;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTree<FPType, N, ElemType, DT>::KDTree(KDTree&& rhs) noexcept
: root(rhs.root), treeSize(rhs.treeSize) {
    rhs.root = nullptr;
    rhs.treeSize = 0;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTree<FPType, N, ElemType, DT>& KDTree<FPType, N, ElemType, DT>::
operator = (KDTree&& rhs) & noexcept {
    if (this != &rhs) {
        delete root;
        root = rhs.root;
        treeSize = rhs.treeSize;
        rhs.root = nullptr;
        rhs.treeSize = 0;
    }
    return *this;
}


template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <class RAI>
void KDTree<FPType, N, ElemType, DT>::
RangeCtorHelper(TreeNode*& cur_ndPtr, std::size_t dim, RAI begin,
                RAI median, RAI end) {
    std::nth_element(begin, median, end, [=](const auto& nh1, const auto& nh2) {
        return nh1.key[dim] < nh2.key[dim];});
    cur_ndPtr = new TreeNode(std::move(median->key),
                            std::move(median->value));
    std::size_t nextDim = dim == N - 1 ? 0 : dim + 1;
    if (begin == median - 1) {
        cur_ndPtr->left = new TreeNode(std::move(begin->key),
                                      std::move(begin->value));
    } else if (begin != median) {
        RangeCtorHelper(cur_ndPtr->left, nextDim, begin,
                        begin + (median-begin)/2, median);
    }
    if (median + 1 == end - 1) {
        cur_ndPtr->right = new TreeNode(std::move((median + 1)->key),
                                       std::move((median+1)->value));
    } else if (median + 1 != end) {
        RangeCtorHelper(cur_ndPtr->right, nextDim, median+1,
                        median + (end-median+1)/2, end);
    }
} 

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTree<FPType, N, ElemType, DT>::treeCopy(TreeNode*& thisNode,
                                       TreeNode* otherNode) {
    if (otherNode) {
        if (thisNode) {
            thisNode->object = otherNode->object;
            thisNode->key = otherNode->key;
        } else {
            thisNode = new TreeNode(otherNode->key, otherNode->object);
        }
        treeCopy(thisNode->left, otherNode->left);
        treeCopy(thisNode->right, otherNode->right);
    } else if (thisNode) {
        delete thisNode;
        thisNode = nullptr;
    }
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
KDTree<FPType, N, ElemType, DT>::~KDTree() {
    delete root;
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
constexpr std::size_t KDTree<FPType, N, ElemType, DT>::dimension() const {
    return N;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename PointND<FPType, N>::DistType KDTree<FPType, N, ElemType, DT>::distType() const {
    return DT;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
std::size_t KDTree<FPType, N, ElemType, DT>::size() const {
    return treeSize;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
int KDTree<FPType, N, ElemType, DT>::height() const {
    return heightHelper(root);
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
int KDTree<FPType, N, ElemType, DT>::heightHelper(TreeNode *n) const {
    return n ? 1 + std::max(heightHelper(n->left),
                            heightHelper(n->right)) : -1;
}


template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTree<FPType, N, ElemType, DT>::empty() const {
    return treeSize == 0;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTree<FPType, N, ElemType, DT>::PrintTreeInfo() const {
    std::cout << "Tree height is " << height()
              << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTree<FPType, N, ElemType, DT>::Clear() {
    delete root;
    treeSize = 0;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTree<FPType, N, ElemType, DT>::
insert(const PointND<FPType, N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        (*ndPtr)->object = value;
    } else {
        *ndPtr = new TreeNode(pt, value);
        treeSize++;
    }
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
bool KDTree<FPType, N, ElemType, DT>::contains(const PointND<FPType, N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType& KDTree<FPType, N, ElemType, DT>::operator[] (const PointND<FPType, N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        treeSize++;
    }
    return (*ndPtr)->object;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType& KDTree<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTree&>(*this).at(pt));
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
const ElemType& KDTree<FPType, N, ElemType, DT>::at(const PointND<FPType, N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename KDTree<FPType, N, ElemType, DT>::TreeNode**
KDTree<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTree*>(this)->findNodePtr(pt));
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
typename KDTree<FPType, N, ElemType, DT>::TreeNode*const*
KDTree<FPType, N, ElemType, DT>::findNodePtr(const PointND<FPType, N>& pt) const {
    TreeNode *const*n = &root;
    for (std::size_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1)
        n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType KDTree<FPType, N, ElemType, DT>::
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

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTree<FPType, N, ElemType, DT>::kNNValueHelper(TreeNode *cur, std::size_t dim,
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
template <class OutputIter>
OutputIter KDTree<FPType, N, ElemType, DT>::NNsWithFence(const PointND<FPType, N>& pt,
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
    
    std::vector<std::pair<FPType, std::pair<const PointND<FPType, N>&, const ElemType&>>> distKVPairs;
    distKVPairs.reserve(sqrt(treeSize));
    FPType bestDistSq = std::numeric_limits<FPType>::max(),
           bestDistDiffSq = std::numeric_limits<FPType>::max(),
           curDistSq, fenceSq = fence*fence;
    std::size_t dim = 0;
    TreeNode *cur = root, *next;
    bool hasNext = true;
    
    
    
    struct ActRecord {
        FPType curDist;
        TreeNode *cur;
        std::size_t dim;
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
                                     std::forward<std::pair<const PointND<FPType, N>&,
                                     const ElemType&>>({cur->key, cur->object}));
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
            *returnIt++ = {kvPairs.first, kvPairs.second};
    }
    return returnIt;
    
}

/*
template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
void KDTree<FPType, N, ElemType, DT>::
NNsWithFenceHelper(TreeNode *cur, std::size_t dim, const PointND<FPType, N>& pt,
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
    std::size_t nextDim = dim + 1 < N ? dim + 1 : 0;
    FPType thisDiff = pt[dim] - cur->key[dim];
    TreeNode *next = thisDiff < 0 ? cur->left : cur->right;
    if (next)
        NNsWithFenceHelper(next, nextDim, pt, diff, distKVPairs, bestDistSq, bestDistDiffSq);
    if (thisDiff*thisDiff < bestDistDiffSq) {
        next = next == cur->left ? cur->right : cur->left;
        if (next)
            NNsWithFenceHelper(next, nextDim, pt, diff, distKVPairs, bestDistSq, bestDistDiffSq);
    }
}*/

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
ElemType KDTree<FPType, N, ElemType, DT>::NNValue(const PointND<FPType, N> &pt) const {
    
    
    FPType bestDist = std::numeric_limits<FPType>::max(), curDist, diff;
    std::size_t dim = 0, nextDim;
    const ElemType *bestValue = nullptr;
    TreeNode *cur = root, *next;
    bool hasNext = true;
    
    struct actRecord {
        FPType curDist;
        TreeNode *cur;
        std::size_t dim;
    };
    actRecord st[32], *it = st;
    while (it != st || hasNext) {
        if (!hasNext) {
            const auto &ar = *--it;
            if (!(ar.curDist < bestDist && (cur = ar.cur)))
                continue;
            dim = ar.dim;
            hasNext = true;
        }
        nextDim = dim == N - 1 ? 0 : dim + 1;
        curDist = PointND<FPType, N>::template
                  dist<PointND<FPType, N>::DistType::EUCSQ>(cur->key, pt);
        if (curDist < bestDist) {
            bestDist = curDist;
            bestValue = &cur->object;
        }
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
    
    /*
    
    FPType bestDist = std::numeric_limits<FPType>::max();
    const ElemType *bestValue = nullptr;
    NNValueHelper(root, 0, pt, bestValue, bestDist);
    */
    
    return bestValue ? *bestValue : ElemType();
    
}

template <typename FPType, std::size_t N, typename ElemType, typename PointND<FPType, N>::DistType DT>
template <typename PointND<FPType, N>::DistType thisDt,
typename std::enable_if<thisDt == PointND<FPType, N>::DistType::EUC, int>::type>
void KDTree<FPType, N, ElemType, DT>::
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
void KDTree<FPType, N, ElemType, DT>::
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
FPType KDTree<FPType, N, ElemType, DT>::branchMin(const PointND<FPType, N> &trPt,
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

#endif // KDTREE_INCLUDED
