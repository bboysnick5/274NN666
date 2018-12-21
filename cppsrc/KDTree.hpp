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
#include <functional>


enum class DistType {
    EUC = 0,
    MAN,
    HAV
};

template <size_t N, typename ElemType, DistType dType>
class KDTree {
public:
    // Constructor: KDTree();
    // Usage: KDTree<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTree.
    KDTree();
    
    // Constructor: KDTree(FwdItType begin, FwdItType end);
    // Usage: KDTree<3, int> myTree(vec.begin(), bec.end());
    // ----------------------------------------------------
    // Constructs a KDTree from a collection. The tree will
    // be balanced using median constructing method
    // NOTE: The tree will not eliminate duplicates and the
    //       intended behavior will not be comprimised, tho
    //       less efficient with extra wasteful space.
    template <class RAI>
    KDTree(RAI begin, RAI end);
    
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
    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);
    
    // KDTree(const KDTree& rhs);
    // KDTree& operator=(const KDTree& rhs);
    // Usage: KDTree<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Move constructor and move assignment operator.
    KDTree(KDTree&& rhs);
    KDTree& operator=(KDTree&& rhs);
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
    size_t dimension() const;
    DistType distType() const;
    
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
    
    // Iter rangeDiffKNNPairs(const Point<N>&, double, Iter) const
    // Usage: Iter end = kd.rangeDiffKNNPairs(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a double offset, return a set of points in the KDTree
    // nearest to p such that the farthest one in the set is at least offset
    // distance close to p than the rest of the points in the tree.
    // The forward iterator is passed in and filled and the end will be returned.
    template <class Iter>
    Iter rangeDiffKNNPairs(const Point<N>&, double, Iter) const;
    
private:
    
    struct TreeNode {
        TreeNode *left;
        TreeNode *right;
        Point<N> key;
        ElemType object;
        
        TreeNode(const Point<N>& k, const ElemType& obj)
        : left(nullptr), right(nullptr), key(k), object(obj) {}
        
        TreeNode(const Point<N>&& k = Point<N>(), const ElemType&& obj = ElemType())
        : left(nullptr), right(nullptr), key(std::move(k)), object(std::move(obj)) {}
        
        ~TreeNode() {
            delete left;
            delete right;
        }
    };
    
    TreeNode *root;
    size_t treeSize;
    
    
    const static std::function<double (const Point<N>& p1, const Point<N> &p2)>
        distFuncs[];
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    size_t heightHelper(TreeNode *n) const;
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void rangeCtorHelper(TreeNode*&, int, RAI, RAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(TreeNode *cur, int level, const Point<N> &pt,
                        BoundedPQueue<ElemType> &bpq) const;
    
    void rangeDiffKNNPairsHelper(TreeNode*, int, const Point<N> &, double,
                                 std::vector<std::pair<double,
                                 std::pair<Point<N>, ElemType>>>&, double&) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const Point<N>& key) const;
    void NNValueHelper(TreeNode *cur, int level, const Point<N> &pt,
                       double &bestDist, ElemType *&bestValue) const;
    
    // ----------------------------------------------------
    // Helper meothod for deep copy
    void treeCopy(TreeNode*& thisNode, TreeNode *otherNode);
    
    // TreeNode** findNodePtr(const Point<N>& pt);
    // TreeNode*const* findNodePtr(const Point<N>& pt) const;
    // Usage: TreeNode **nodePtr = findNodePtr(pt);
    // ----------------------------------------------------
    // Returns the pointer pointing to the node address
    // corresponding to the given Point. In this double pointing
    // fashion, we can construct a node at that location.
    TreeNode** findNodePtr(const Point<N>& pt);
    TreeNode*const* findNodePtr(const Point<N>& pt) const;
    
    inline double branchMin(const Point<N>&, const Point<N>&, int) const;

};

/** KDTree class implementation details */


template <size_t N, typename ElemType, DistType dType>
const std::function<double (const Point<N> &p1, const Point<N> &p2)>
KDTree<N, ElemType, dType>::distFuncs[] = {Point<N>::euclDist,
    Point<N>::manhDist, Point<N>::havDist};


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, DistType dType>
KDTree<N, ElemType, dType>::KDTree() : root(nullptr), treeSize(0) {}

template <size_t N, typename ElemType, DistType dType>
template <class RAI>
KDTree<N, ElemType, dType>::KDTree(RAI begin, RAI end) : treeSize(0) {
    rangeCtorHelper(root, 0, begin, end);
}

template <size_t N, typename ElemType, DistType dType>
KDTree<N, ElemType, dType>::KDTree(const KDTree& rhs)
: root(new TreeNode()), treeSize(rhs.treeSize) {
    treeCopy(root, rhs.root);
}

template <size_t N, typename ElemType, DistType dType>
KDTree<N, ElemType, dType>&
KDTree<N, ElemType, dType>::operator=(const KDTree& rhs) {
    if (this != &rhs) {
        treeSize = rhs.treeSize;
        root = new TreeNode();
        treeCopy(root, rhs.root);
    }
    return *this;
}

template <size_t N, typename ElemType, DistType dType>
KDTree<N, ElemType, dType>::KDTree(KDTree&& rhs)
: root(rhs.root), treeSize(rhs.treeSize) {
    rhs.root = nullptr;
    rhs.treeSize = 0;
}

template <size_t N, typename ElemType, DistType dType>
KDTree<N, ElemType, dType>& KDTree<N, ElemType, dType>::operator=(KDTree&& rhs){
    if (this != &rhs) {
        delete root;
        root = rhs.root;
        treeSize = rhs.treeSize;
        rhs.root = nullptr;
        rhs.treeSize = 0;
    }
    return *this;
}

template <size_t N, typename ElemType, DistType dType>
template <class RAI>
void KDTree<N, ElemType, dType>::rangeCtorHelper(TreeNode*& curNdPtr, int level,
                                                 RAI begin, RAI end) {
    RAI median;
    if (begin != end) {
        median = begin + (end-begin)/2;
        std::nth_element(begin, median, end, [=](const std::pair<Point<N>,
            ElemType>& p1, const std::pair<Point<N>, ElemType>& p2) {
            return p1.first[level%N] < p2.first[level%N];});
        curNdPtr = new TreeNode(std::move(median->first),
                                std::move(median->second));
        treeSize++;
        rangeCtorHelper(curNdPtr->left, level+1, begin, median);
        rangeCtorHelper(curNdPtr->right, level+1, median+1, end);
    }
}

template <size_t N, typename ElemType, DistType dType>
void KDTree<N, ElemType, dType>::treeCopy(TreeNode*& thisNode,
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

template <size_t N, typename ElemType, DistType dType>
KDTree<N, ElemType, dType>::~KDTree() {
    delete root;
    root = nullptr;
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, DistType dType>
size_t KDTree<N, ElemType, dType>::dimension() const {
    return N;
}

template <size_t N, typename ElemType, DistType dType>
DistType KDTree<N, ElemType, dType>::distType() const {
    return dType;
}

template <size_t N, typename ElemType, DistType dType>
size_t KDTree<N, ElemType, dType>::size() const {
    return treeSize;
}

template <size_t N, typename ElemType, DistType dType>
size_t KDTree<N, ElemType, dType>::height() const {
    return heightHelper(root);
}

template <size_t N, typename ElemType, DistType dType>
size_t KDTree<N, ElemType, dType>::heightHelper(TreeNode *n) const {
    return n? 1 + std::max(heightHelper(n->left), heightHelper(n->right)) : -1;
}


template <size_t N, typename ElemType, DistType dType>
bool KDTree<N, ElemType, dType>::empty() const {
    return treeSize == 0;
}


// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, DistType dType>
void KDTree<N, ElemType, dType>::insert(const Point<N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        (*ndPtr)->object = value;
    } else {
        *ndPtr = new TreeNode(pt, value);
        treeSize++;
    }
}

template <size_t N, typename ElemType, DistType dType>
bool KDTree<N, ElemType, dType>::contains(const Point<N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <size_t N, typename ElemType, DistType dType>
ElemType& KDTree<N, ElemType, dType>::operator[] (const Point<N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        treeSize++;
    }
    return (*ndPtr)->object;
}

template <size_t N, typename ElemType, DistType dType>
ElemType& KDTree<N, ElemType, dType>::at(const Point<N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTree&>(*this).at(pt));
}

template <size_t N, typename ElemType, DistType dType>
const ElemType& KDTree<N, ElemType, dType>::at(const Point<N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <size_t N, typename ElemType, DistType dType>
typename KDTree<N, ElemType, dType>::TreeNode**
KDTree<N, ElemType, dType>::findNodePtr(const Point<N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTree*>(this)->findNodePtr(pt));
}

template <size_t N, typename ElemType, DistType dType>
typename KDTree<N, ElemType, dType>::TreeNode*const*
KDTree<N, ElemType, dType>::findNodePtr(const Point<N>& pt) const {
    int dim = N-1;
    TreeNode *const*n = &root;
    while (*n && (*n)->key != pt && ++dim) {
        n = pt[dim%N] < (*n)->key[dim%N] ? &(*n)->left : &(*n)->right;
    }
    return n;
}


template <size_t N, typename ElemType, DistType dType>
ElemType KDTree<N, ElemType, dType>::kNNValue(const Point<N>& pt, size_t k) const {
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType> bpq(k);
    kNNValueHelper(root, 0, pt, bpq);
    
    std::multimap<int, ElemType, std::greater<int>> freqMap;
    while (!bpq.empty()) {
        ElemType elem = bpq.dequeueMin();
        for (typename std::multimap<int, ElemType>::iterator
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

template <size_t N, typename ElemType, DistType dType>
void KDTree<N, ElemType, dType>::kNNValueHelper(TreeNode *cur, int level,
const Point<N>& pt, BoundedPQueue<ElemType> &bpq) const {
    if (!cur)
        return;
    bpq.enqueue(cur->object, distFuncs[(int)dType](cur->key, pt));
    TreeNode *next = pt[level%N] < cur->key[level%N] ? cur->left : cur->right;
    kNNValueHelper(next, level + 1, pt, bpq);
    if (bpq.size() < bpq.maxSize()
        || branchMin(cur->key, pt, level%N) < bpq.worst()) {
        TreeNode *other = next == cur->left ? cur->right : cur->left;
        kNNValueHelper(other, level + 1, pt, bpq);
    }
}

template <size_t N, typename ElemType, DistType dType>
template <class Iter>
Iter KDTree<N, ElemType, dType>::rangeDiffKNNPairs(const Point<N>& key,
                                            double diff, Iter it) const {
    std::vector<std::pair<double, std::pair<Point<N>, ElemType>>> distKVPairs;
    distKVPairs.reserve(sqrt(treeSize));
    double best = std::numeric_limits<double>::max();
    rangeDiffKNNPairsHelper(root, 0, key, diff, distKVPairs, best);
    for (const auto& p : distKVPairs) {
        if (p.first < best + diff)
            *it++ = std::move(p.second);
    }
    return it;
}

template <size_t N, typename ElemType, DistType dType>
void KDTree<N, ElemType, dType>::rangeDiffKNNPairsHelper(TreeNode *cur, int level,
const Point<N>& pt, double diff, std::vector<std::pair<double,
std::pair<Point<N>, ElemType>>> &distKVPairs, double& bestDist) const {
    if (!cur)
        return;
    auto dist = distFuncs[static_cast<int>(dType)](cur->key, pt);
    if (dist < bestDist) {
        bestDist = dist;
        distKVPairs.push_back(std::move(distKVPairs[0]));
        distKVPairs[0] = std::make_pair(dist,
                                        std::make_pair(cur->key, cur->object));
    } else if (dist < bestDist + diff) {
        distKVPairs.emplace_back(dist, std::make_pair(cur->key, cur->object));
    }
    TreeNode *next = pt[level%N] < cur->key[level%N] ? cur->left : cur->right;
    rangeDiffKNNPairsHelper(next, level+1, pt, diff, distKVPairs, bestDist);
    if (branchMin(cur->key, pt, level%N) < bestDist+diff) {
        TreeNode *other = next == cur->left ? cur->right : cur->left;
        rangeDiffKNNPairsHelper(other, level+1, pt, diff, distKVPairs, bestDist);
    }
}

template <size_t N, typename ElemType, DistType dType>
ElemType KDTree<N, ElemType, dType>::NNValue(const Point<N> &pt) const {
    double bestDist = std::numeric_limits<double>::max();
    ElemType *valuePtr = nullptr;
    NNValueHelper(root, 0, pt, bestDist, valuePtr);
    return valuePtr == nullptr ? ElemType() : *valuePtr;
}

template <size_t N, typename ElemType, DistType dType>
void KDTree<N, ElemType, dType>::NNValueHelper(TreeNode *cur, int level,
            const Point<N> &pt, double &bestDist, ElemType *&bestValue) const {
    if (!cur)
        return;
    double curDist = distFuncs[(int)dType](cur->key, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->object;
    }
    TreeNode *next = pt[level%N] < cur->key[level%N] ? cur->left : cur->right;
    NNValueHelper(next, level+1, pt, bestDist, bestValue);
    if (branchMin(cur->key, pt, level%N) < bestDist) {
        TreeNode *other = next == cur->left ? cur->right : cur->left;
        NNValueHelper(other, level+1, pt, bestDist, bestValue);
    }
}

template <size_t N, typename ElemType, DistType dType>
double KDTree<N, ElemType, dType>::branchMin(const Point<N> &trPt,
const Point<N> &searchPt, int idx) const {
    switch (dType) {
        case DistType::EUC:
        case DistType::MAN:
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

#endif // KDTREE_INCLUDED
