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


template <size_t N, typename ElemType, typename Point<N>::DistType DT
= Point<N>::DistType::EUC>
class KDTreeExpandLongestVec {
public:
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
    ~KDTreeExpandLongestVec() = default;
    
    // KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs);
    // KDTreeExpandLongestVec& operator=(const KDTreeExpandLongestVec& rhs);
    // Usage: KDTreeExpandLongestVec<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Copy constructor and copy assignment operator.
    KDTreeExpandLongestVec(const KDTreeExpandLongestVec&) = default;
    KDTreeExpandLongestVec& operator=(const KDTreeExpandLongestVec&)& = default;
    
    // KDTreeExpandLongestVec(const KDTreeExpandLongestVec& rhs);
    // KDTreeExpandLongestVec& operator=(const KDTreeExpandLongestVec& rhs);
    // Usage: KDTreeExpandLongestVec<3, int> one = two;
    // Usage: one = two;
    // -----------------------------------------------------
    // Move constructor and move assignment operator.
    KDTreeExpandLongestVec(KDTreeExpandLongestVec&&) noexcept = default;
    KDTreeExpandLongestVec& operator=(KDTreeExpandLongestVec&&)& noexcept = default;
    
    // size_t dimension() const;
    // Usage: size_t dim = kd.dimension();
    // ----------------------------------------------------
    // Returns the dimension of the points stored in this KDTreeExpandLongestVec.
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
    int height() const;
    bool empty() const;
    
    void clear();
    
    void printTreeInfo() const;
    
    // bool contains(const Point<N>& pt) const;
    // Usage: if (kd.contains(pt))
    // ----------------------------------------------------
    // Returns whether the specified point is contained in the KDTreeExpandLongestVec.
    bool contains(const Point<N>&) const;
    
    // void insert(const Point<N>& pt, const ElemType& value);
    // Usage: kd.insert(v, "This value is associated with v.");
    // ----------------------------------------------------
    // Inserts the point pt into the KDTreeExpandLongestVec, associating it with the specified
    // value. If the element already existed in the tree, the new value will
    // overwrite the existing one.
    void insert(const Point<N>&, const ElemType&);
    
    // ElemType& operator[](const Point<N>& pt);
    // Usage: kd[v] = "Some Value";
    // ----------------------------------------------------
    // Returns a reference to the value associated with point pt in the KDTreeExpandLongestVec.
    // If the point does not exist, then it is added to the KDTreeExpandLongestVec using the
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
    // Given a point v and an integer k, finds the k points in the KDTreeExpandLongestVec
    // nearest to v and returns the most common value associated with those
    // points. In the event of a tie, one of the most frequent value will be
    // chosen.
    ElemType kNNValue(const Point<N>& key, size_t k) const;
    
    // Iter rangeDiffKNNPairs(const Point<N>&, double, Iter) const
    // Usage: Iter end = kd.rangeDiffKNNPairs(pt, 0.33, begin);
    // ----------------------------------------------------
    // Given a point p and a double offset, return a set of points in the KDTreeExpandLongestVec
    // nearest to p such that the farthest one in the set is at least offset
    // distance close to p than the rest of the points in the tree.
    // The forward iterator is passed in and filled and the end will be returned.
    template <class Iter>
    Iter rangeDiffKNNPairs(const Point<N>&, double, Iter) const;
    
private:
    
    
    struct TreeNode {
        int dimToExpand;
        Point<N> key;
        unsigned int rightIdx;
        //ElemType object;
        
        TreeNode() = default;
        // TreeNode(const TreeNode&) = default;
        // TreeNode& operator = (const TreeNode&) = default;
        // TreeNode(TreeNode&&) = default;
        // TreeNode& operator = (TreeNode&&) = default;
        
        TreeNode(int dimToExpand, const Point<N>& k)
        : dimToExpand(dimToExpand), key(k), rightIdx(0) {}
        
        TreeNode(int dimToExpand, Point<N>&& k)
        : dimToExpand(dimToExpand), key(std::move(k)), rightIdx(0) {}
        
        ~TreeNode() = default;
    };
    
    std::vector<TreeNode> ndVec;
    std::vector<ElemType> objectVec;
    
    // ----------------------------------------------------
    // Helper method for finding the height of a tree
    int heightHelper(const TreeNode *n) const;
    
    // ----------------------------------------------------
    // Helper method for range constructor
    template <class RAI>
    void rangeCtorHelper(RAI, RAI, RAI, std::array<double, N*2>&);
    
    
    template <class RAI>
    static std::array<double, N*2> computeInitBBox(RAI, RAI);
    
    // ----------------------------------------------------
    // Helper method for kNNValue search
    void kNNValueHelper(TreeNode *cur, size_t dim, const Point<N> &pt,
                        BoundedPQueue<ElemType> &bpq) const;
    
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

/** KDTreeExpandLongestVec class implementation details */


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
KDTreeExpandLongestVec<N, ElemType, DT>::KDTreeExpandLongestVec(Const_RAI cbegin, Const_RAI cend) {
    ndVec.reserve(cend - cbegin);
    objectVec.reserve(cend - cbegin);
    std::vector<std::pair<Point<N>, ElemType>> constructData(cbegin, cend);
    auto bbox = computeInitBBox(cbegin, cend);
    rangeCtorHelper(constructData.begin(), constructData.begin() +
                    (constructData.end() - constructData.begin())/2,
                    constructData.end(), bbox);
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <typename RAI, typename std::enable_if<std::is_same<typename
std::iterator_traits<RAI>::iterator_category,
std::random_access_iterator_tag>::value && !std::is_const<typename
std::remove_pointer< typename std::iterator_traits<RAI>::pointer>::type>::value, int>::type>
KDTreeExpandLongestVec<N, ElemType, DT>::KDTreeExpandLongestVec(RAI begin, RAI end) {
    ndVec.reserve(end - begin);
    objectVec.reserve(end - begin);
    auto bbox = computeInitBBox(begin, end);
    rangeCtorHelper(begin, begin + (end - begin)/2, end, bbox);
    
    /*
     struct actRecord {
     TreeNode** curNdPtr;
     RAI thisBeginIt, median, thisEndIt;
     };
     
     actRecord st[static_cast<size_t>(log2(treeSize+1))], *it = st;
     bool hasChild = true;
     RAI thisBeginIt = begin, thisEndIt = end, median = thisBeginIt + (thisEndIt-thisBeginIt)/2;
     while (it != st || hasChild) {
     if (!hasChild) {
     hasChild = true;
     *(--it)->curNdPtr = ndPoolPtr;
     thisBeginIt = it->thisBeginIt;
     thisEndIt = it->thisEndIt;
     median = it->median;
     }
     
     size_t dim = 0;
     double maxSpan = bbox[1] - bbox[0];
     for (size_t i = 1; i < N; ++i) {
     size_t bboxLowIdx = i * 2, bboxHighIdx = bboxLowIdx + 1;
     double span = bbox[bboxHighIdx] - bbox[bboxLowIdx];
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
     TreeNode* curNd = ndPoolPtr++;
     
     if (thisBeginIt != median) {
     if (median + 1 != thisEndIt) {
     *it++ = {&curNd->right, median + 1,
     median + (thisEndIt-median+1)/2, thisEndIt};
     }
     bbox[dim*2+1] = curNd->key[dim];
     thisEndIt = median;
     median = thisBeginIt + (median-thisBeginIt)/2;
     curNd->left = ndPoolPtr;
     } else if (median + 1 != thisEndIt) {
     bbox[dim*2] = curNd->key[dim];
     thisBeginIt = median+1;
     median += (thisEndIt-median+1)/2;
     curNd->right = ndPoolPtr;
     } else {
     hasChild = false;
     }
     } */
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <class RAI>
std::array<double, N*2> KDTreeExpandLongestVec<N, ElemType, DT>::
computeInitBBox(RAI begin, RAI end) {
    std::array<double, N*2> bbox;
    std::generate(bbox.begin(), bbox.end(),
                  [lowHighToggle = 0lu, lowHigh = std::array<double, 2>{
                                                  std::numeric_limits<double>::max(),
                                                  std::numeric_limits<double>::min()}]() mutable {
                      return lowHigh[lowHighToggle++%2];});
    std::for_each(begin, end, [&](const auto &p) mutable {
        for (size_t i = 0; i < N; ++i) {
            double ptValOnithDim = p.first[i];
            auto &bboxLow = bbox[i*2], &bboxHigh = bbox[i*2+1];
            bboxLow = std::min(bboxLow, ptValOnithDim);
            bboxHigh = std::max(bboxHigh, ptValOnithDim);
        }
    });
    return bbox;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
template <class RAI>
void KDTreeExpandLongestVec<N, ElemType, DT>::
rangeCtorHelper(RAI begin, RAI median, RAI end,
                std::array<double, N*2> &bbox) {
    size_t dim = 0;
    double maxSpan = bbox[1] - bbox[0];
    for (size_t i = 1; i != N; ++i) {
        size_t bboxLowIdx = i * 2, bboxHighIdx = bboxLowIdx + 1;
        double span = bbox[bboxHighIdx] - bbox[bboxLowIdx];
        if (span > maxSpan) {
            maxSpan = span;
            dim = i;
        }
    }
    
    std::nth_element(begin, median, end, [=](const auto& p1, const auto& p2) {
        return p1.first[dim] < p2.first[dim];});
    ndVec.emplace_back(dim, std::move(median->first));
    objectVec.emplace_back(std::move(median->second));
    TreeNode* curNdPtr = &ndVec[ndVec.size()-1];
    
    if (begin == median - 1) {
        //curNdPtr->left = &*ndVec.end();
        ndVec.emplace_back(-1, std::move(begin->first));
        objectVec.emplace_back(std::move(begin->second));
    } else if (begin != median) {
        auto prevDimHigh = bbox[dim*2+1];
        bbox[dim*2+1] = curNdPtr->key[dim];
        //bbox[dim*2+1] = std::max_element(begin, median, [dim](const auto &p1, const auto&p2){return p1.first[dim] < p2.first[dim];})->first[dim];
        rangeCtorHelper(begin, begin + (median-begin)/2, median, bbox);
        bbox[dim*2+1] = prevDimHigh;
    }
    
    if (median + 2 == end) {
        curNdPtr->rightIdx = static_cast<unsigned int>(ndVec.size());
        ndVec.emplace_back(-1, std::move((median+1)->first));
        objectVec.emplace_back(std::move((median+1)->second));
    } else if (median + 1 != end) {
        curNdPtr->rightIdx = static_cast<unsigned int>(ndVec.size());
        auto prevDimLow = bbox[dim*2];
        bbox[dim*2] = curNdPtr->key[dim];
        //bbox[dim*2] = std::min_element(median + 1, end, [dim](const auto &p1, const auto &p2){return p1.first[dim] < p2.first[dim];})->first[dim];
        rangeCtorHelper(median+1, median + (end-median+1)/2, end, bbox);
        bbox[dim*2] = prevDimLow;
    }
}

// ----------------------------------------------------------
// ----------------- TREE INFORMATION  ----------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTreeExpandLongestVec<N, ElemType, DT>::dimension() const {
    return N;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
typename Point<N>::DistType KDTreeExpandLongestVec<N, ElemType, DT>::distType() const {
    return DT;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
size_t KDTreeExpandLongestVec<N, ElemType, DT>::size() const {
    return ndVec.size();
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
int KDTreeExpandLongestVec<N, ElemType, DT>::height() const {
    const TreeNode *root = &ndVec[0];
    return heightHelper(root);
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
int KDTreeExpandLongestVec<N, ElemType, DT>::heightHelper(const TreeNode *n) const {
    //return n ? 1 + std::max(heightHelper(n->left), heightHelper(n->right)) : -1;
    return 1;
}


template <size_t N, typename ElemType, typename Point<N>::DistType DT>
bool KDTreeExpandLongestVec<N, ElemType, DT>::empty() const {
    return ndVec.empty();
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTreeExpandLongestVec<N, ElemType, DT>::printTreeInfo() const {
    std::cout << "Tree height is " << height()
    << "\nTree size is " << size() << "\n";
}

// ----------------------------------------------------------
// ----------------- MODIFIERS AND ACCESS -------------------
// ----------------------------------------------------------

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTreeExpandLongestVec<N, ElemType, DT>::clear() {
    ndVec.clear();
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTreeExpandLongestVec<N, ElemType, DT>::
insert(const Point<N>& pt, const ElemType& value) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (*ndPtr) {
        //(*ndPtr)->object = value;
    } else {
        //*ndPtr = new TreeNode(0, pt, value);
        //treeSize++;
    }
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
bool KDTreeExpandLongestVec<N, ElemType, DT>::contains(const Point<N>& pt) const {
    return *findNodePtr(pt) != nullptr;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType& KDTreeExpandLongestVec<N, ElemType, DT>::operator[] (const Point<N>& pt) {
    TreeNode **ndPtr = findNodePtr(pt);
    if (!*ndPtr) {
        *ndPtr = new TreeNode(pt, ElemType());
        //treeSize++;
    }
    return (*ndPtr)->object;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType& KDTreeExpandLongestVec<N, ElemType, DT>::at(const Point<N>& pt) {
    return const_cast<ElemType&>(static_cast<const KDTreeExpandLongestVec&>(*this).at(pt));
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
const ElemType& KDTreeExpandLongestVec<N, ElemType, DT>::at(const Point<N>& pt) const {
    TreeNode *const*n = findNodePtr(pt);
    if (!*n) {
        //throw out_of_range("The point is out of range");
    }
    return (*n)->object;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
typename KDTreeExpandLongestVec<N, ElemType, DT>::TreeNode**
KDTreeExpandLongestVec<N, ElemType, DT>::findNodePtr(const Point<N>& pt) {
    return const_cast<TreeNode**>(static_cast<const KDTreeExpandLongestVec*>(this)->findNodePtr(pt));
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
typename KDTreeExpandLongestVec<N, ElemType, DT>::TreeNode*const*
KDTreeExpandLongestVec<N, ElemType, DT>::findNodePtr(const Point<N>& pt) const {
    //TreeNode *const*n = &root;
    TreeNode *const*n;
    for (size_t dim = 0; *n && (*n)->key != pt; dim = dim == N - 1 ? 0 : dim+1){}
        //n = pt[dim] < (*n)->key[dim] ? &(*n)->left : &(*n)->right;
    return n;
}


template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType KDTreeExpandLongestVec<N, ElemType, DT>::
kNNValue(const Point<N>& pt, size_t k) const {
    //if (empty())
    //    return ElemType();
    if (k == 1)
        return NNValue(pt);
    
    BoundedPQueue<ElemType> bpq(k);
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

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
void KDTreeExpandLongestVec<N, ElemType, DT>::kNNValueHelper(TreeNode *cur, size_t dim,
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
Iter KDTreeExpandLongestVec<N, ElemType, DT>::
rangeDiffKNNPairs(const Point<N>& pt, double fence, Iter returnIt) const {
    std::vector<std::tuple<double, const Point<N>&, const ElemType&>> distPtElemTuple;
    distPtElemTuple.reserve(std::sqrt(ndVec.size()));
    std::pair<double, const TreeNode*> st[static_cast<size_t>(log2(ndVec.size()+1))],
                                       *it = st;
    double bestDistSq = std::numeric_limits<double>::max(),
           bestDistDiffSq = bestDistSq, fenceSq = fence*fence;
    const TreeNode *cur = ndVec.data();
    
    while (true) {
        double curDistSq = Point<N>::template
        dist<Point<N>::DistType::EUCSQ>(cur->key, pt);
        if (curDistSq < bestDistDiffSq) {
            if (curDistSq < bestDistSq) {
                bestDistSq = curDistSq;
                bestDistDiffSq = bestDistSq + fenceSq + 2*fence*sqrt(bestDistSq);
            }
            distPtElemTuple.emplace_back(curDistSq, cur->key, objectVec[cur-ndVec.data()]);
        }
        
        if (unsigned int rightIdx = cur->rightIdx) {
            int dim = cur->dimToExpand;
            double diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(it++) std::pair<double, const TreeNode*>(diff*diff, ndVec.data() + rightIdx);
                ++cur;
            } else {
                new(it++) std::pair<double, const TreeNode*>(diff*diff, cur+1);
                cur = ndVec.data() + rightIdx;
            }
        } else if (cur++->dimToExpand == -1) {
            do {
                if (it == st)
                    goto FINAL;
            } while ((--it)->first >= bestDistDiffSq || !(cur = it->second));
        }
    }
    
FINAL:
    for (const auto &[distSq, pt, elem] : distPtElemTuple) {
        if (distSq < bestDistDiffSq)
            *returnIt++ = {pt, elem};
    }
    return returnIt;
}

template <size_t N, typename ElemType, typename Point<N>::DistType DT>
ElemType KDTreeExpandLongestVec<N, ElemType, DT>::NNValue(const Point<N> &pt) const {
    
    std::pair<double, const TreeNode*> st[static_cast<size_t>(log2(ndVec.size()+1))],
                                       *it = st;
    double bestDist = std::numeric_limits<double>::max();
    const ElemType *bestValue = nullptr;
    const TreeNode *cur = ndVec.data();
    
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
        double curDist = Point<N>::template dist<Point<N>::DistType::EUCSQ>(cur->key, pt);
        if (curDist < bestDist) {
            bestDist = curDist;
            bestValue = &objectVec[cur-ndVec.data()];
        }
        
        // LOGGGGGGGGGGGGGGG
        //if (logCondition)
        // thisNumNodesSearches++;
        
        if (unsigned int rightIdx = cur->rightIdx) {
            int dim = cur->dimToExpand;
            double diff = pt[dim] - cur->key[dim];
            if (diff < 0.0) {
                new(it++) std::pair<double, const TreeNode*>(diff*diff, ndVec.data() + rightIdx);
                ++cur;
            } else {
                new(it++) std::pair<double, const TreeNode*>(diff*diff, cur+1);
                cur = ndVec.data() + rightIdx;
            }
        } else if (cur++->dimToExpand == -1) {
            do {
                if (it == st)
                    return bestValue ? *bestValue : ElemType();
            } while ((--it)->first >= bestDist || !(cur = it->second));
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
void KDTreeExpandLongestVec<N, ElemType, DT>::
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
void KDTreeExpandLongestVec<N, ElemType, DT>::
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
double KDTreeExpandLongestVec<N, ElemType, DT>::branchMin(const Point<N> &trPt,
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




#endif /* KDTreeExpandLongestVec_hpp */
