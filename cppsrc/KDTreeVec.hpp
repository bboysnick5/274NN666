/**
 * File: KDTree.h
 * Author: Yunlong Nick Liu
 * ------------------------
 * An interface representing a kd-tree in some number of dimensions. The tree
 * can be constructed from a set of data and then queried for membership and
 * nearest neighbors.
 */

#ifndef KDTreeVec_INCLUDED
#define KDTreeVec_INCLUDED

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
#include <memory>
#include <vector>


template <size_t N, typename ElemType, DistType dType>
class KDTreeVec {
public:
    // Constructor: KDTreeVec();
    // Usage: KDTreeVec<3, int> myTree;
    // ----------------------------------------------------
    // Constructs an empty KDTreeVec.
    KDTreeVec() = default;
    
    // Constructor: KDTreeVec(FwdItType begin, FwdItType end);
    // Usage: KDTreeVec<3, int> myTree(vec.begin(), bec.end());
    // ----------------------------------------------------
    // Constructs a KDTreeVec from a collection. The tree will
    // be balanced using median constructing method
    // NOTE: The tree will no eliminate duplicates and the
    //       intended behavior will not be comprimised, tho
    //       less efficient with extra wasteful space.
    template <class RAI>
    KDTreeVec(RAI begin, RAI end);
    
    // Destructor: ~KDTreeVec()
    // Usage: (implicit)
    // ----------------------------------------------------
    // Cleans up all resources used by the KDTreeVec.
    virtual ~KDTreeVec() = default;
    
    // KDTreeVec(const KDTreeVec& rhs);
    // KDTreeVec& operator=(const KDTreeVec& rhs);
    // Usage: KDTreeVec<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Copy constructor and copy assignment operator.
    KDTreeVec(const KDTreeVec& rhs) = default;
    KDTreeVec& operator=(const KDTreeVec& rhs) = default;
    
    // KDTreeVec(const KDTreeVec& rhs);
    // KDTreeVec& operator=(const KDTreeVec& rhs);
    // Usage: KDTreeVec<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Move constructor and move assignment operator.
    KDTreeVec(KDTreeVec&& rhs) = default;
    KDTreeVec& operator=(KDTreeVec&& rhs) = default;
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTreeVec.
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
    // Returns whether the specified point is contained in the KDTreeVec.
    bool contains(const Point<N>& pt) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTreeVec, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>& pt, const ElemType& value);
    
    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTreeVec.
    // If the point does not exist, then it is added to the KDTreeVec using the
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
    // Given a point v and an integer k, finds the k points in the KDTreeVec
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, size_t k) const;
    
    template <class Iter>
    Iter rangeDiffKNNPairs(const Point<N>&, double, Iter) const;
    
private:
    
    using TrIt = typename std::vector<std::shared_ptr
                 <std::pair<Point<N>, ElemType>>>::iterator;
    
    size_t treeSize;
    std::vector<std::shared_ptr<std::pair<Point<N>, ElemType>>> _tv;
    
    const static std::function<double (const Point<N>& p1, const Point<N> &p2)>
    distFuncs[];
    
    template <class RAI>
    void rangeCtorHelper(int idx, int level, RAI begin, RAI end);
    
    size_t heightHelper(int idx) const;

    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(int, int, const Point<N>&,
                        BoundedPQueue<ElemType>&) const;
    
    void rangeDiffKNNPairsHelper(int, int, const Point<N> &, double,
                                 std::vector<std::pair<double,
                                 std::pair<Point<N>, ElemType>>>&, double&) const;
    
    // ----------------------------------------------------
    // Identical to kNNValue method with k equals 1. NNValue
    // and its helper method are used to speed up the search
    // when finding the nearest neighbor only.
    ElemType NNValue(const Point<N>& key) const;
    void NNValueHelper(int, int, const Point<N>&, double&, ElemType*&) const;
    
    // TreeNode** findNodePtr(const Point<N>& pt);
    // TreeNode*const* findNodePtr(const Point<N>& pt) const;
    // Usage: TreeNode **nodePtr = findNodePtr(pt);
    // ----------------------------------------------------
    // Returns the pointer pointing to the node address
    // corresponding to the given Point. In this double pointing
    // fashion, we can construct a node at that location.
    inline const TrIt findNdItr(const Point<N>& pt) const;
    inline TrIt findNdItr(const Point<N>& pt);

    
    inline double branchMin(const Point<N>&, const Point<N>&, int) const;
};

/** KDTreeVec class implementation details */


template <size_t N, typename ElemType, DistType dType>
const std::function<double (const Point<N> &p1, const Point<N> &p2)>
KDTreeVec<N, ElemType, dType>::distFuncs[] = {Point<N>::euclDist,
    Point<N>::manhDist, Point<N>::havDist};


// ----------------------------------------------------------
// --------------------- CONSTRUCTORS -----------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, DistType dType>
template <class RAI>
KDTreeVec<N, ElemType, dType>::KDTreeVec(RAI begin, RAI end) : treeSize(0) {
    _tv.resize((end - begin)*2);
    rangeCtorHelper(0, 0, begin, end);
    size_t i = _tv.size() - 1;
    for (; !_tv[i]; --i) ;
    _tv.resize(i+1);
}

template <size_t N, typename ElemType, DistType dType>
template <class RAI>
void KDTreeVec<N, ElemType, dType>::
rangeCtorHelper(int idx, int level, RAI begin, RAI end) {
    RAI median;
    if (begin != end) {
        median = begin + (end - begin)/2;
        std::nth_element(begin, median, end, [=](const std::pair<Point<N>,
            ElemType>& p1, const std::pair<Point<N>, ElemType>& p2) {
            return p1.first[level%N] < p2.first[level%N];});
        _tv[idx] = std::make_shared<std::pair<Point<N>, ElemType>>(std::move(*median));
        treeSize++;
        rangeCtorHelper(2*idx+1, level+1, begin, median);
        rangeCtorHelper(2*idx+2, level+1, median+1, end);
    }
}


// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, DistType dType>
size_t KDTreeVec<N, ElemType, dType>::dimension() const {
    return N;
}

template <size_t N, typename ElemType, DistType dType>
DistType KDTreeVec<N, ElemType, dType>::distType() const {
    return dType;
}

template <size_t N, typename ElemType, DistType dType>
size_t KDTreeVec<N, ElemType, dType>::size() const {
    return treeSize;
}

template <size_t N, typename ElemType, DistType dType>
size_t KDTreeVec<N, ElemType, dType>::height() const {
    return heightHelper(0);
}

template <size_t N, typename ElemType, DistType dType>
size_t KDTreeVec<N, ElemType, dType>::heightHelper(int idx) const {
    return _tv[idx] ? 1 + std::max(heightHelper(2*idx+1), heightHelper(2*idx+2))
                    : -1;
}


template <size_t N, typename ElemType, DistType dType>
bool KDTreeVec<N, ElemType, dType>::empty() const {
    return treeSize == 0;
}


// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, DistType dType>
void KDTreeVec<N, ElemType, dType>::insert(const Point<N>& pt, const ElemType& value) {
    auto ndPtr = findNdItr(pt);
    if (*ndPtr) {
        (*ndPtr)->second = value;
    } else {
        *ndPtr = std::make_shared<std::pair<Point<N>, ElemType>>
                 (std::forward_as_tuple(pt, value));
        treeSize++;
    }
}

template <size_t N, typename ElemType, DistType dType>
bool KDTreeVec<N, ElemType, dType>::contains(const Point<N>& pt) const {
    return *findNdItr(pt) != nullptr;
}

template <size_t N, typename ElemType, DistType dType>
ElemType& KDTreeVec<N, ElemType, dType>::operator[] (const Point<N>& pt) {
    auto ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = std::make_shared<std::pair<Point<N>, ElemType>>
                 (std::forward_as_tuple(pt, ElemType()));;
        treeSize++;
    }
    return (*ndPtr)->object;
}

template <size_t N, typename ElemType, DistType dType>
ElemType& KDTreeVec<N, ElemType, dType>::at(const Point<N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeVec&>(*this).at(pt));
}

template <size_t N, typename ElemType, DistType dType>
const ElemType& KDTreeVec<N, ElemType, dType>::at(const Point<N>& pt) const {
    auto *const*ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        throw out_of_range("The point is out of range");
    }
    return (*ndPtr)->object;
}

template <size_t N, typename ElemType, DistType dType>
typename KDTreeVec<N, ElemType, dType>::TrIt
KDTreeVec<N, ElemType, dType>::findNdItr(const Point<N>& pt) {
    return const_cast<TrIt>(static_cast<const KDTreeVec*>(this)->findNdIt(pt));
}

template <size_t N, typename ElemType, DistType dType>
const typename KDTreeVec<N, ElemType, dType>::TrIt
KDTreeVec<N, ElemType, dType>::findNdItr(const Point<N>& pt) const {
    int level = N-1, idx = 0;
    std::shared_ptr<std::pair<Point<N>, ElemType>> cur;
    while ((cur = _tv[idx]) && cur && cur->first != pt && ++level) {
        idx = pt[level%N] < cur->first[level%N] ? 2*idx+1 : 2*idx+2;
    }
    return _tv.begin() + idx;
}


template <size_t N, typename ElemType, DistType dType>
ElemType KDTreeVec<N, ElemType, dType>::kNNValue(const Point<N>& pt, size_t k) const {
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType> bpq(k);
    kNNValueHelper(0, 0, pt, bpq);
    
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
void KDTreeVec<N, ElemType, dType>::kNNValueHelper(int curIdx, int level,
const Point<N>& pt, BoundedPQueue<ElemType> &bpq) const {
    if (curIdx >= _tv.size() || !_tv[curIdx])
        return;
    std::shared_ptr<std::pair<Point<N>, ElemType>> cur = _tv[curIdx];
    bpq.enqueue(cur->second, distFuncs[(int)dType](cur->first, pt));
    int next = pt[level%N] < cur->first[level%N] ? curIdx*2+1 : curIdx*2+2;
    kNNValueHelper(next, level + 1, pt, bpq);
    if (bpq.size() < bpq.maxSize()
        || branchMin(cur->first, pt, level%N) < bpq.worst()) {
        int other = next == curIdx*2+1 ? next+1 : next-1;
        kNNValueHelper(other, level + 1, pt, bpq);
    }
}

template <size_t N, typename ElemType, DistType dType>
template <class Iter>
Iter KDTreeVec<N, ElemType, dType>::rangeDiffKNNPairs(const Point<N>& key,
                                                   double diff, Iter it) const {
    std::vector<std::pair<double, std::pair<Point<N>, ElemType>>> distKVPairs;
    distKVPairs.reserve(sqrt(treeSize));
    double best = std::numeric_limits<double>::max();
    rangeDiffKNNPairsHelper(0, 0, key, diff, distKVPairs, best);
    for (const auto& p : distKVPairs) {
        if (p.first < best + diff)
            *it++ = std::move(p.second);
    }
    return it;
}

template <size_t N, typename ElemType, DistType dType>
void KDTreeVec<N, ElemType, dType>::rangeDiffKNNPairsHelper(int curIdx,
int level, const Point<N>& pt, double diff, std::vector<std::pair<double,
std::pair<Point<N>, ElemType>>> &distKVPairs, double& bestDist) const {
    if (curIdx >= _tv.size() || !_tv[curIdx])
        return;
    std::shared_ptr<std::pair<Point<N>, ElemType>> cur = _tv[curIdx];
    auto dist = distFuncs[static_cast<int>(dType)](cur->first, pt);
    if (dist < bestDist + diff) {
        if (dist < bestDist) {
            bestDist = dist;
            distKVPairs.push_back(std::move(distKVPairs[0]));
            distKVPairs[0] = std::make_pair(dist, std::make_pair(cur->first,
                                                                 cur->second));
        } else {
            distKVPairs.emplace_back(dist, std::make_pair(cur->first,
                                                          cur->second));
        }
    }
    int next = pt[level%N] < cur->first[level%N] ? curIdx*2+1 : curIdx*2+2;
    rangeDiffKNNPairsHelper(next, level+1, pt, diff, distKVPairs, bestDist);
    if (branchMin(cur->first, pt, level%N) < bestDist+diff) {
        int other = next == curIdx*2+1 ? next+1 : next-1;
        rangeDiffKNNPairsHelper(other, level+1, pt, diff, distKVPairs, bestDist);
    }
}

template <size_t N, typename ElemType, DistType dType>
ElemType KDTreeVec<N, ElemType, dType>::NNValue(const Point<N> &pt) const {
    double bestDist = std::numeric_limits<double>::max();
    ElemType *valuePtr = nullptr;
    NNValueHelper(0, 0, pt, bestDist, valuePtr);
    return valuePtr == nullptr ? ElemType() : *valuePtr;
}

template <size_t N, typename ElemType, DistType dType>
void KDTreeVec<N, ElemType, dType>::NNValueHelper(int curIdx, int level,
const Point<N> &pt, double &bestDist, ElemType *&bestValue) const {
    if (curIdx >= _tv.size() || !_tv[curIdx])
        return;
    std::shared_ptr<std::pair<Point<N>, ElemType>> cur = _tv[curIdx];
    double curDist = distFuncs[(int)dType](cur->first, pt);
    if (curDist < bestDist) {
        bestDist = curDist;
        bestValue = &cur->second;
    }
    int next = pt[level%N] < cur->first[level%N] ? curIdx*2+1 : curIdx*2+2;
    NNValueHelper(next, level+1, pt, bestDist, bestValue);
    if (branchMin(cur->first, pt, level%N) < bestDist) {
        int other = next == curIdx*2+1 ? next+1 : next-1;
        NNValueHelper(other, level+1, pt, bestDist, bestValue);
    }
}

template <size_t N, typename ElemType, DistType dType>
double KDTreeVec<N, ElemType, dType>::branchMin(const Point<N> &trPt,
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

#endif // KDTreeVec_INCLUDED
