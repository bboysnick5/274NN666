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
#include <cmath>
#include <iostream>
#include <unordered_map>
#include <map>
#include <utility>
#include <set>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>


template <size_t N, typename ElemType>
class KDTree {
public:
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    KDTree();
    
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    
    template <class FwdItType>
    KDTree(FwdItType begin, FwdItType end);
    
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
    // Deep-copies the contents of another KDTree into this one.
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
    size_t dimension() const;
    
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
    
    // bool contains(const Point<N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTree.
    bool contains(const Point<N>& pt) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTree, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>& pt, const ElemType& value);
    
    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTree.
    // If the point does not exist, then it is added to the KDTree using the
    // default value of ElemType as its key.
    ElemType& operator[](const Point<N>& pt);
    
    // ElemType& at(const Point<N>& pt);
    // const ElemType& at(const Point<N>& pt) const;
    // Usage: cout << kd.at(v) << endl;
    // ----------------------------------------------------
    // Returns a reference to the key associated with the point pt. If the point
    // is not in the tree, this function throws an out_of_range exception.
    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;
    
    // ElemType kNNValue(const Point<N>& key, size_t k) const
    // Usage: cout << kd.kNNValue(v, 3) << endl;
    // ----------------------------------------------------
    // Given a point v and an integer k, finds the k points in the KDTree
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, size_t k) const;
    
    
private:
    
    struct TreeNode {
        TreeNode *left;
        TreeNode *right;
        Point<N> key;
        ElemType object;
        
        TreeNode(const Point<N>& k = Point<N>(), const ElemType& obj = ElemType())
        : left(nullptr), right(nullptr), key(k), object(obj) {}
        
        ~TreeNode() {
            delete left;
            delete right;
        }
    };
    
    TreeNode *root;
    size_t treeSize;
    
    size_t heightHelper(TreeNode *n) const;
    
    template <class FwdItType>
    void rangeCtorHelper(TreeNode*& cur, int level, FwdItType begin, FwdItType end);
    
    void kNNValueHelper(TreeNode *cur, int level, const Point<N> &pt,
                        BoundedPQueue<ElemType> &bpq) const;
    
    // Identical to kNNValue method with k equals 1. This is used to speed up
    // the search when only need to find the nearest neighbor.
    ElemType NNValue(const Point<N>& key) const;
    void NNValueHelper(TreeNode *cur, int level, const Point<N> &pt,
                       double &bestDist, ElemType *&bestValue) const;
    void treeCopy(TreeNode *thisNode, TreeNode *otherNode);
    
    TreeNode** findNodePtr(const Point<N>& pt);
    TreeNode*const* findNodePtr(const Point<N>& pt) const;


};

/** KDTree class implementation details */

// ------------ CONSTRUCTOR AND DESTRUCTOR ------------------

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree() : root(nullptr), treeSize(0) {}

template <size_t N, typename ElemType>
template <class FwdItType>
KDTree<N, ElemType>::KDTree(FwdItType begin, FwdItType end) {
    rangeCtorHelper(root, 0, begin, end);
}

template <size_t N, typename ElemType>
template <class FwdItType>
void KDTree<N, ElemType>::rangeCtorHelper(TreeNode*& curNdPtr, int level,
                                          FwdItType begin, FwdItType end) {
    if (begin != end) {
        FwdItType median = begin + (end-begin)/2;
        std::nth_element(begin, median, end, [=](const std::pair<Point<N>,
            ElemType>& p1, const std::pair<Point<N>, ElemType>& p2) {
            return p1.first[level%N] < p2.first[level%N];});
        curNdPtr = new TreeNode(median->first, median->second);
        rangeCtorHelper(curNdPtr->left, level+1, begin, median);
        rangeCtorHelper(curNdPtr->right, level+1, median+1, end);
    }
}


template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree() {
    delete root;
    root = nullptr;
}


template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const {
    return N;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs)
    : root(new TreeNode()), treeSize(rhs.treeSize) {
    treeCopy(root, rhs.root);
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs) {
    if (this != &rhs) {
        treeSize = rhs.treeSize;
        root = new TreeNode();
        treeCopy(root, rhs.root);
    }
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const {
    return treeSize;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::height() const {
    return heightHelper(root);
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::heightHelper(TreeNode *n) const {
    return n? 1 + std::max(heightHelper(n->left), heightHelper(n->right)) : -1;
}


template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const {
    return treeSize == 0;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        (*ndPtr)->object = value;
    } else {
        *ndPtr = new TreeNode(pt, value);
        treeSize++;
    }
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[] (const Point<N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        treeSize++;
    }
    return (*ndPtr)->object;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTree&>(*this).at(pt));
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::TreeNode**
KDTree<N, ElemType>::findNodePtr(const Point<N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTree*>(this)->findNodePtr(pt));
}

template <size_t N, typename ElemType>
typename KDTree<N, ElemType>::TreeNode*const*
KDTree<N, ElemType>::findNodePtr(const Point<N>& pt) const {
    int level = N-1;
    TreeNode *const*n = &root;
    while (*n && (*n)->key != pt && ++level) {
        n = pt[level%N] < (*n)->key[level%N] ? &(*n)->left : &(*n)->right;
    }
    return n;
}


template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::kNNValue(const Point<N>& pt, size_t k) const {
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType> bpq(k);
    kNNValueHelper(root, 0, pt, bpq);
    
    std::multimap<int, ElemType, std::greater<int>> freqMap;
    while (!bpq.empty()) {
        ElemType elem = bpq.dequeueMin();
        for (typename std::multimap<int, ElemType>::iterator it = freqMap.begin();
             it != freqMap.end(); ++it) {
            if (it->second == elem) {
                freqMap.insert(std::make_pair(it->first+1, it->second));
                freqMap.erase(it);
                goto outer;
            }
        }
        freqMap.insert(std::make_pair(1, elem));
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

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::kNNValueHelper(TreeNode *cur, int level, const Point<N>& pt,
                                         BoundedPQueue<ElemType> &bpq) const {
    if (!cur)
        return;
    bpq.enqueue(cur->object, Distance(cur->key, pt));
    TreeNode *next = pt[level%N] < cur->key[level%N] ? cur->left : cur->right;
    kNNValueHelper(next, level + 1, pt, bpq);
    if (bpq.size() < bpq.maxSize()
        || std::fabs(pt[level%N] - cur->key[level%N]) < bpq.worst()) {
        TreeNode *other = next == cur->left ? cur->right : cur->left;
        kNNValueHelper(other, level + 1, pt, bpq);
    }
}

template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::NNValue(const Point<N> &pt) const {
    double bestDist = std::numeric_limits<double>::max();
    ElemType *valuePtr = nullptr;
    NNValueHelper(root, 0, pt, bestDist, valuePtr);
    return valuePtr == nullptr ? ElemType() : *valuePtr;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::NNValueHelper(TreeNode *cur, int level,
            const Point<N> &pt, double &bestDist, ElemType *&bestValue) const {
    if (!cur)
        return;
    double curDist = Distance(cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    TreeNode *next = pt[level%N] < cur->key[level%N] ? cur->left : cur->right;
    NNValueHelper(next, level+1, pt, bestDist, bestValue);
    if (std::fabs(pt[level%N] - cur->key[level%N]) < bestDist) {
        TreeNode *other = next == cur->left ? cur->right : cur->left;
        NNValueHelper(other, level+1, pt, bestDist, bestValue);
    }
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::treeCopy(TreeNode *thisNode, TreeNode *otherNode) {
    thisNode->object = otherNode->object;
    thisNode->key = otherNode->key;
    if (otherNode->left) {
        if (!thisNode->left)
            thisNode->left = new TreeNode();
        treeCopy(thisNode->left, otherNode->left);
    } else if (thisNode->left) {
        delete thisNode->left;
        thisNode->left = nullptr;
    }
    if (otherNode->right) {
        if (!thisNode->right)
            thisNode->right = new TreeNode();
        treeCopy(thisNode->right, otherNode->right);
    } else if (thisNode->right) {
        delete thisNode->right;
        thisNode->left = nullptr;
    }
}




#endif // KDTREE_INCLUDED
