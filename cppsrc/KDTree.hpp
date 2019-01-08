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


template <size_t N, typename ElemType, typename Point<N>::DistType DT
                                                = Point<N>::DistType::EUC>
class KDTree {
public:
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
    template <class RAI>
    KDTree(RAI, RAI);
    
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
    KDTree& operator=(KDTree&&)&;
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTree.
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
    // Returns whether the specified point is contained in the KDTree.
    bool contains(const Point<N>&) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTree, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>&, const ElemType&);
    
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
        
        TreeNode() = default;
       // TreeNode(const TreeNode&) = default;
       // TreeNode& operator = (const TreeNode&) = default;
       // TreeNode(TreeNode&&) = default;
       // TreeNode& operator = (TreeNode&&) = default;

        TreeNode(const Point<N>& k, const ElemType& obj)
        : left(nullptr), right(nullptr), key(k), object(obj) {}
        
        TreeNode(const Point<N>&& k, const ElemType&& obj)
        : left(nullptr), right(nullptr), key(std::move(k)),
          object(std::move(obj)) {}
        
        ~TreeNode() {
            delete left;
            delete right;
        }
    };
    
    TreeNode *root;
    size_t treeSize;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    size_t heightHelper(TreeNode *n) const;
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void rangeCtorHelper(TreeNode*&, size_t, RAI, RAI, RAI);
    
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
    
    double branchMin(const Point<N>&, const Point<N>&, size_t) const;

};

/** KDTree class implementation details */


// ----------------------------------------------------------
// ----------------------- BIG FIVE -------------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <class RAI>
KDTree<N, ElemType, DT>::KDTree(RAI begin, RAI end) : treeSize(end-begin) {
    
    if (treeSize == 1) {
        root = new TreeNode(std::move(begin->first),
                            std::move(begin->second));
    } else if (treeSize > 1) {
        rangeCtorHelper(root, 0, begin, end, begin + (end - begin)/2);
    }
     

    
    /*
    struct actRecord {
        TreeNode** curNdPtr;
        size_t dim;
        RAI thisBeginIt, thisEndIt, median;
    };
    
    actRecord st[static_cast<size_t>(log2(treeSize+1))];
    actRecord *it = st;
    bool hasChild = true;
    TreeNode** curNdPtr = &root;
    size_t dim = 0, nextDim = dim == N - 1 ? 0 : dim + 1;
    RAI thisBeginIt = begin, thisEndIt = end, median = thisBeginIt + (thisEndIt-thisBeginIt)/2;
    while (it != st || hasChild) {
        if (!hasChild) {
            hasChild = true;
            --it;
            curNdPtr = it->curNdPtr;
            dim = it->dim;
            thisBeginIt = it->thisBeginIt;
            thisEndIt = it->thisEndIt;
            median = it->median;
            //std::tie(curNdPtr, dim, thisBeginIt, thisEndIt, median) = *--it;
        }
        
        std::nth_element(thisBeginIt, median, thisEndIt,
                         [=](const std::pair<Point<N>, ElemType>& p1,
                                   const std::pair<Point<N>, ElemType>& p2) {
                                        return p1.first[dim] < p2.first[dim];});
        *curNdPtr = new TreeNode(std::move(median->first),
                                 std::move(median->second));
        nextDim = dim == N - 1 ? 0 : dim + 1;
        
       if (median != thisBeginIt) {
            if (median + 1 != thisEndIt) {
                *it++ = {&(*curNdPtr)->left, nextDim, thisBeginIt, median, thisBeginIt + (median-thisBeginIt)/2};
                curNdPtr = &(*curNdPtr)->right;
                thisBeginIt = median+1;
                median = median + (thisEndIt-median+1)/2;
                dim = nextDim;
            } else {
                curNdPtr = &(*curNdPtr)->left;
                thisEndIt = median;
                median = thisBeginIt + (median-thisBeginIt)/2;
                dim = nextDim;
            }
       } else if (median + 1 != thisEndIt) {
           curNdPtr = &(*curNdPtr)->right;
           thisBeginIt = median+1;
           median = median + (thisEndIt-median+1)/2;
           dim = nextDim;
       } else {
            hasChild = false;
       }
    } */
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTree<N, ElemType, DT>::KDTree(const KDTree& rhs)
: root(new TreeNode()), treeSize(rhs.treeSize) {
    treeCopy(root, rhs.root);
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTree<N, ElemType, DT>&
KDTree<N, ElemType, DT>::operator=(const KDTree& rhs) & {
    if (this != &rhs) {
        delete root;
        treeSize = rhs.treeSize;
        root = new TreeNode();
        treeCopy(root, rhs.root);
    }
    return *this;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTree<N, ElemType, DT>::KDTree(KDTree&& rhs) noexcept
: root(rhs.root), treeSize(rhs.treeSize) {
    rhs.root = nullptr;
    rhs.treeSize = 0;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTree<N, ElemType, DT>& KDTree<N, ElemType, DT>::operator=(KDTree&& rhs) & {
    if (this != &rhs) {
        delete root;
        root = rhs.root;
        treeSize = rhs.treeSize;
        rhs.root = nullptr;
        rhs.treeSize = 0;
    }
    return *this;
}


template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <class RAI>
void KDTree<N, ElemType, DT>::
rangeCtorHelper(TreeNode*& curNdPtr, size_t dim, RAI begin,
                RAI end, RAI median) {
    std::nth_element(begin, median, end, [=](const auto& p1, const auto& p2) {
        return p1.first[dim] < p2.first[dim];});
    curNdPtr = new TreeNode(std::move(median->first),
                            std::move(median->second));
    size_t nextDim = dim == N - 1 ? 0 : dim + 1;
    if (begin == median - 1) {
        curNdPtr->left = new TreeNode(std::move(begin->first),
                                      std::move(begin->second));
    } else if (begin != median) {
        rangeCtorHelper(curNdPtr->left, nextDim, begin,
                        median, begin + (median-begin)/2);
    }
    if (median + 1 == end - 1) {
        curNdPtr->right = new TreeNode(std::move((median + 1)->first),
                                       std::move((median+1)->second));
    } else if (median + 1 != end) {
        rangeCtorHelper(curNdPtr->right, nextDim, median+1,
                        end, median + (end-median+1)/2);
    }
} 

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTree<N, ElemType, DT>::treeCopy(TreeNode*& thisNode,
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

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
KDTree<N, ElemType, DT>::~KDTree() {
    delete root;
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTree<N, ElemType, DT>::dimension() const {
    return N;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
typename Point<N>::DistType KDTree<N, ElemType, DT>::distType() const {
    return DT;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTree<N, ElemType, DT>::size() const {
    return treeSize;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTree<N, ElemType, DT>::height() const {
    return heightHelper(root);
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTree<N, ElemType, DT>::heightHelper(TreeNode *n) const {
    return n ? 1 + std::max(heightHelper(n->left),
                            heightHelper(n->right)) : -1;
}


template <size_t N, typename ElemType, typename Point<N>::DistType DT>
bool KDTree<N, ElemType, DT>::empty() const {
    return treeSize == 0;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTree<N, ElemType, DT>::printTreeInfo() const {
    std::cout << "Tree height is " << height() << std::endl
              << "Tree size is " << size() << std::endl;
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTree<N, ElemType, DT>::clear() {
    delete root;
    treeSize = 0;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTree<N, ElemType, DT>::
insert(const Point<N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        (*ndPtr)->object = value;
    } else {
        *ndPtr = new TreeNode(pt, value);
        treeSize++;
    }
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
bool KDTree<N, ElemType, DT>::contains(const Point<N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType& KDTree<N, ElemType, DT>::operator[] (const Point<N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        treeSize++;
    }
    return (*ndPtr)->object;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType& KDTree<N, ElemType, DT>::at(const Point<N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTree&>(*this).at(pt));
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
const ElemType& KDTree<N, ElemType, DT>::at(const Point<N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
typename KDTree<N, ElemType, DT>::TreeNode**
KDTree<N, ElemType, DT>::findNodePtr(const Point<N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTree*>(this)->findNodePtr(pt));
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
typename KDTree<N, ElemType, DT>::TreeNode*const*
KDTree<N, ElemType, DT>::findNodePtr(const Point<N>& pt) const {
    TreeNode *const*n = &root;
    for (size_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1)
        n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType KDTree<N, ElemType, DT>::
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
void KDTree<N, ElemType, DT>::kNNValueHelper(TreeNode *cur, size_t dim,
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
Iter KDTree<N, ElemType, DT>::rangeDiffKNNPairs(const Point<N>& pt,
                                                double fence, Iter returnIt) const {
    /*
    std::vector<std::pair<double, std::pair<Point<N>, ElemType>>> distKVPairs;
    distKVPairs.reserve(sqrt(treeSize));
    double bestSq = std::numeric_limits<double>::max(),
           bestDiffSq = std::numeric_limits<double>::max();
    rangeDiffKNNPairsHelper(root, 0, key, diff, distKVPairs, bestSq, bestDiffSq);
    for (const auto &p : distKVPairs) {
        if (p.first < bestDiffSq)
            *it++ = std::move(p.second);
    }
    return it;
    */
    
    std::vector<std::pair<double, std::pair<const Point<N>*, const ElemType*>>> distKVPairs;
    distKVPairs.reserve(sqrt(treeSize));
    double bestDistSq = std::numeric_limits<double>::max(),
           bestDistDiffSq = std::numeric_limits<double>::max(),
           curDistSq, fenceSq = fence*fence;
    size_t dim = 0;
    TreeNode *cur = root, *next;
    bool hasNext = true;
    
    
    struct actRecord {
        TreeNode *cur;
        size_t dim;
        double curDist;
    };
    actRecord st[static_cast<size_t>(log2(treeSize+1))], *it = st;
    while (it != st || hasNext) {
        if (!hasNext) {
            const auto &ar = *--it;
            if (!(ar.curDist < bestDistDiffSq && (cur = ar.cur)))
                continue;
            dim = ar.dim;
            hasNext = true;
        }
        curDistSq = Point<N>::template
        dist<Point<N>::DistType::EUCSQ>(cur->key, pt);
        if (curDistSq < bestDistDiffSq) {
            if (curDistSq < bestDistSq) {
                bestDistSq = curDistSq;
                bestDistDiffSq = bestDistSq + fenceSq + 2*fence*sqrt(bestDistSq);
            }
            distKVPairs.emplace_back(curDistSq, std::make_pair(&cur->key, &cur->object));
        }
        curDistSq = pt[dim] - cur->key[dim];
        next = curDistSq < 0.0 ? cur->left : cur->right;
        curDistSq *= curDistSq;
        dim = dim == N - 1 ? 0 : dim + 1;
        if (next) {
            *it++ = {next == cur->left ? cur->right : cur->left, dim, curDistSq};
            cur = next;
            //dim = nextDim;
        } else {
            if (curDistSq < bestDistDiffSq) {
                next = next == cur->left ? cur->right : cur->left;
                if (next) {
                    cur = next;
                    //dim = nextDim;
                    continue;
                }
            }
            hasNext = false;
        }
    }
    
    
    
    for (const auto &p : distKVPairs) {
        if (p.first < bestDistDiffSq) {
            const auto &pSec = p.second;
            *returnIt++ = {*pSec.first, *pSec.second};
        }
    }
    return returnIt;
    
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTree<N, ElemType, DT>::
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
ElemType KDTree<N, ElemType, DT>::NNValue(const Point<N> &pt) const {
    
    
    double bestDist = std::numeric_limits<double>::max(), curDist, diff;
    size_t dim = 0, nextDim;
    const ElemType *bestValue = nullptr;
    TreeNode *cur = root, *next;
    bool hasNext = true;
    
    struct actRecord {
        TreeNode *cur;
        size_t dim;
        double curDist;
    };
    actRecord st[static_cast<size_t>(log2(treeSize+1))], *it = st;
    while (it != st || hasNext) {
        if (!hasNext) {
            const auto &ar = *--it;
            if (!(ar.curDist < bestDist && (cur = ar.cur)))
                continue;
            dim = ar.dim;
            hasNext = true;
        }
        nextDim = dim == N - 1 ? 0 : dim + 1;
        curDist = Point<N>::template
                  dist<Point<N>::DistType::EUCSQ>(cur->key, pt);
        if (curDist < bestDist) {
            bestDist = curDist;
            bestValue = &cur->object;
        }
        diff = pt[dim] - cur->key[dim];
        next = diff < 0.0 ? cur->left : cur->right;
        curDist = diff*diff;
        if (next) {
            *it++ = {next == cur->left ? cur->right : cur->left, nextDim, curDist};
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
    
    double bestDist = std::numeric_limits<double>::max();
    const ElemType *bestValue = nullptr;
    NNValueHelper(root, 0, pt, bestValue, bestDist);
    */
    
    return bestValue ? *bestValue : ElemType();
    
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <typename Point<N>::DistType thisDt,
typename std::enable_if<thisDt == Point<N>::DistType::EUC, int>::type>
void KDTree<N, ElemType, DT>::
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
void KDTree<N, ElemType, DT>::
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
double KDTree<N, ElemType, DT>::branchMin(const Point<N> &trPt,
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

#endif // KDTREE_INCLUDED
