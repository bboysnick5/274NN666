//
//  Utility.hpp
//  274F16NearestSB
//
//  Created by nick on 3/1/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef Utility_hpp
#define Utility_hpp

#include "Definition.hpp"
#include <stdio.h>
//#include <execution>
#include <omp.h>
#include <utility>


/*
 Allow pass-in value calculation formula and store the computed data
 to save up one value computation as it were two times
 as if provided via comparator.
 */

namespace utility {

// cycle swap
template<typename T>
constexpr void CycleSwap(T& t1, T& t2, T& t3) {
    T temp(std::move(t1));
    t1 = std::move(t2);
    t2 = std::move(t3);
    t3 = std::move(temp);
}

template <class Compare, class ForwardIterator>
void Sort3(ForwardIterator x, ForwardIterator y, ForwardIterator z, Compare c)
{
    if (!c(*y, *x))             // if x <= y
    {
        if (!c(*z, *y))         // if y <= z
            return;             // x <= y && y <= z
                                // x <= y && y > z
        std::swap(*y, *z);      // x <= z && y < z
        if (c(*y, *x))          // if x > y
            std::swap(*x, *y);  // x < y && y <= z
        return;
    }
    if (c(*z, *y)) {            // x > y, if y > z
        std::swap(*x, *z);      // x < y && y < z
        return;
    }
    std::swap(*x, *y);          // x > y && y <= z
                                // x < y && x <= z
    if (c(*z, *y))              // if y > z
        std::swap(*y, *z);      // x <= y && y < z
}


template <class ForwardIterator, class GetDist, class Compare>
static ForwardIterator
MinElementGivenDistFunc(ForwardIterator first, ForwardIterator last, GetDist distFunc, Compare comp,
                        def::PolicyTag<def::ThreadingPolicy::kSingle>) {
    if (first != last) [[likely]] {
        ForwardIterator i = first;
        auto bestDist = distFunc(*first);
        while (++i != last) {
            if (auto dist = distFunc(*i); comp(dist, bestDist)) {
                bestDist = dist;
                first = i;
            }
        }
    }
    return first;
}

template <class ForwardIterator, class GetDist, class Compare>
static ForwardIterator
MinElementGivenDistFunc(ForwardIterator first, ForwardIterator last, GetDist distFunc, Compare comp,
                        def::PolicyTag<def::ThreadingPolicy::kMultiOmp>) {
    if (first == last) [[unlikely]]
        return first;
    auto result = first;
    auto bestDist = distFunc(*first++);
    #pragma omp parallel
    {
        auto this_thread_best_dist = bestDist;
        ForwardIterator this_thread_best_it = first;
        #pragma omp for
        for (ForwardIterator i = first; i < last; ++i) {
            if (auto this_dist = distFunc(*i); comp(this_dist, this_thread_best_dist)) {
                this_thread_best_dist = this_dist;
                this_thread_best_it = i;
            }
        }
        #pragma omp critical
        {
            if (comp(this_thread_best_dist, bestDist)) {
                bestDist = this_thread_best_dist;
                result = this_thread_best_it;
            }
        }
    }
    return result;
    /*
    if (first == last) [[unlikely]]
        return first;
    std::ptrdiff_t result_idx = 0;
    auto bestDist = distFunc(*first++);
    auto last_idx = last - first;
    #pragma omp parallel
    {
        auto this_thread_best_dist = bestDist;
        std::ptrdiff_t this_thread_best_idx = 0;
        #pragma omp for
        for (std::ptrdiff_t i = 0; i < last_idx; ++i) {
            if (auto this_dist = distFunc(*(first+i)); comp(this_dist, this_thread_best_dist)) {
                this_thread_best_dist = this_dist;
                this_thread_best_idx = i;
            }
        }
        #pragma omp critical
        {
            if (comp(this_thread_best_dist, bestDist)) {
                bestDist = this_thread_best_dist;
                result_idx = this_thread_best_idx;
            }
        }
    }
    return first + result_idx; // this version is compatible with openmp < 4.0, but requires RAI.
    */
}


template <def::ThreadingPolicy policy = def::ThreadingPolicy::kSingle,
          class ForwardIterator, class GetDist, class Compare>
static ForwardIterator
MinElementGivenDistFunc(ForwardIterator first, ForwardIterator last, GetDist distFunc, Compare comp) {
    return MinElementGivenDistFunc(first, last, distFunc, comp, def::PolicyTag<policy>{});
}

};
#endif /* Utility_hpp */
