//
//  Utility.hpp
//  274F16NearestSB
//
//  Created by nick on 3/1/19.
//  Copyright © 2019 Yunlong Liu. All rights reserved.
//

#ifndef Utility_hpp
#define Utility_hpp

#include <stdio.h>
//#include <execution>
#include <omp.h>

/*
 Allow pass-in value calculation formula and store the computed data
 to save up one value computation as it were two times
 as if provided via comparator.
 */

namespace utility {

// cycle swap
template<typename T>
constexpr void swap(T& t1, T& t2, T& t3) {
    T temp(std::move(t1));
    t1 = std::move(t2);
    t2 = std::move(t3);
    t3 = std::move(temp);
}

template <class _ForwardIterator, class _GetDist, class _Compare>
static _ForwardIterator
MinElementGivenDistFunc(_ForwardIterator __first, _ForwardIterator __last, _GetDist __distFunc, _Compare __comp) {
    if (__first != __last) [[likely]] {
        _ForwardIterator __i = __first;
        auto __bestDist = __distFunc(*__first);
        while (++__i != __last) {
            if (auto __dist = __distFunc(*__i); __comp(__dist, __bestDist)) {
                __bestDist = __dist;
                __first = __i;
            }
        }
    }
    return __first;
}

template <class _ForwardIterator, class _GetDist, class _Compare>
static _ForwardIterator
MinElementGivenDistFunc_p(_ForwardIterator __first, _ForwardIterator __last, _GetDist __distFunc, _Compare __comp) {
    if (__first == __last) [[unlikely]]
        return __first;
    auto __result = __first;
    auto __bestDist = __distFunc(*__first++);
    #pragma omp parallel
    {
        auto this_thread_best_dist = __bestDist;
        _ForwardIterator this_thread_best_it = __first;
        #pragma omp for
        for (_ForwardIterator __i = __first; __i < __last; ++__i) {
            if (auto this_dist = __distFunc(*__i); __comp(this_dist, this_thread_best_dist)) {
                this_thread_best_dist = this_dist;
                this_thread_best_it = __i;
            }
        }
        #pragma omp critical
        {
            if (__comp(this_thread_best_dist, __bestDist)) {
                __bestDist = this_thread_best_dist;
                __result = this_thread_best_it;
            }
        }
    }
    return __result;
        /*
        if (__first == __last) [[unlikely]]
            return __first;
        std::ptrdiff_t __result_idx = 0;
        auto __bestDist = __distFunc(*__first++);
        auto __last_idx = __last - __first;
        #pragma omp parallel
        {
            auto this_thread_best_dist = __bestDist;
            std::ptrdiff_t this_thread_best_idx = 0;
            #pragma omp for
            for (std::ptrdiff_t __i = 0; __i < __last_idx; ++__i) {
                if (auto this_dist = __distFunc(*(__first+__i)); __comp(this_dist, this_thread_best_dist)) {
                    this_thread_best_dist = this_dist;
                    this_thread_best_idx = __i;
                }
            }
            #pragma omp critical
            {
                if (__comp(this_thread_best_dist, __bestDist)) {
                    __bestDist = this_thread_best_dist;
                    __result_idx = this_thread_best_idx;
                }
            }
        }
        return __first + __result_idx; // this version is compatible with openmp < 4.0, but requires RAI.
         */
    }
};
#endif /* Utility_hpp */
