//
//  Algorithm.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 7/23/21.
//  Copyright © 2021 Yunlong Liu. All rights reserved.
//

#ifndef Algorithm_hpp
#define Algorithm_hpp

#include "Definition.hpp"
#include "Point.hpp"

#include <Vc/Vc>

#include <stdio.h>
//#include <execution>
#include <omp.h>
#include <utility>
#include <type_traits>


template <def::ThreadingPolicy Policy = def::ThreadingPolicy::kSingle,
          def::DistType DT = def::DistType::kEucSq>
class Algo {
    
public:
    // concept refactor
    template <class ForwardIterator, class PTType, def::ThreadingPolicy P = Policy,
              std::enable_if_t<P == def::ThreadingPolicy::kSingle, bool> = true>
    inline static ForwardIterator
    LinearNNSearch(ForwardIterator begin, ForwardIterator end, const PTType& pt) {
        if (begin != end) [[likely]] {
            ForwardIterator i = begin;
            auto bestDist = pt.template dist<DT>(*begin);
            while (++i != end) {
                if (auto dist = pt.template dist<DT>(*i); dist < bestDist) {
                    bestDist = dist;
                    begin = i;
                }
            }
        }
        return begin; 
    }
    
    template <class ForwardIterator, class PTType, class DistFunc, def::ThreadingPolicy P = Policy,
              std::enable_if_t<P == def::ThreadingPolicy::kSingle, bool> = true>
    inline static ForwardIterator
    LinearNNSearch(ForwardIterator begin, ForwardIterator end, const PTType& pt, DistFunc df) {
        if (begin != end) [[likely]] {
            ForwardIterator i = begin;
            auto bestDist = df(pt, *begin);
            while (++i != end) {
                if (auto dist = df(pt, *i); dist < bestDist) {
                    bestDist = dist;
                    begin = i;
                }
            }
        }
        return begin;
    }
    
    /*
    // Used to indicate if a passed in iterator is simdized
    // With C++ 20 we can do requires {T.Size;}
    template <typename T, typename = int>
    struct HasSizeVar : std::false_type { };
    
    template <typename T>
    struct HasSizeVar <T, decltype((void) T::Size, 0)> : std::true_type { };
     */
    
    /*
    template <class SimdizedRAI, class PTType,
    std::enable_if_t<Policy == def::ThreadingPolicy::kSimdAoSoA && DT == def::DistType::kEucSq, bool> = true>
    inline static SimdizedRAI
    LinearNNSearch(SimdizedRAI begin, SimdizedRAI end, const PTType& pt) {
        if (begin != end) [[likely]] {
            SimdizedRAI i = begin;
            auto bestDist = pt.template dist<DT>(*begin);
            while (++i != end) {
                if (auto dist = pt.template dist<DT>(*i); dist < bestDist) {
                    bestDist = dist;
                    begin = i;
                }
            }
        }
        return begin;
    } */
    
    
    template <class RAI, class PTType, def::ThreadingPolicy P = Policy,
    std::enable_if_t<P == def::ThreadingPolicy::kSimdSoA &&
    DT == def::DistType::kEucSq &&
    std::is_same_v<typename std::iterator_traits<RAI>::iterator_category,
    std::random_access_iterator_tag>, bool> = true>
    inline static RAI
    LinearNNSearch(RAI begin, RAI end, const PTType& pt) {
        if (begin != end) [[likely]] {
            Vc::simdize<RAI> i = begin;
            auto bestDist = pt.template dist<DT>(*i);;
            while (++i != end) {
                if (auto dist = pt.template dist<DT>(*i); dist < bestDist) {
                    bestDist = dist;
                    begin = i;
                }
            }
        }
        return begin;
    }
    
    template <class RAI, class PTType, def::ThreadingPolicy P = Policy,
    std::enable_if_t<P == def::ThreadingPolicy::kSimdSoA &&
    DT == def::DistType::kHavComp &&
    std::is_same_v<typename std::iterator_traits<RAI>::iterator_category,
    std::random_access_iterator_tag>, bool> = true>
    inline static RAI
    LinearNNSearch(RAI begin, RAI end, const PTType& pt) {
        if (begin != end) [[likely]] {
            Vc::simdize<RAI> i = begin;
            auto bestDist = pt.template dist<DT>(*i);;
            while (++i != end) {
                if (auto dist = pt.template dist<DT>(*i); dist < bestDist) {
                    bestDist = dist;
                    begin = i;
                }
            }
        }
        return begin;
    }
    
    template <class ForwardIterator, class PTType, def::ThreadingPolicy P = Policy,
              std::enable_if_t<P == def::ThreadingPolicy::kMultiOmp, bool> = true>
    inline static ForwardIterator
    LinearNNSearch(ForwardIterator begin, ForwardIterator end, const PTType& pt) {
        if (begin == end) [[unlikely]]
            return begin;
        auto result = begin;
        auto bestDist = pt.template dist<DT>(*begin++);
        #pragma omp parallel
        {
            auto this_thread_best_dist = bestDist;
            ForwardIterator this_thread_best_it = begin;
            #pragma omp for
            for (ForwardIterator i = begin; i < end; ++i) {
                if (auto this_dist = pt.template dist<DT>(*i); this_dist < this_thread_best_dist) {
                    this_thread_best_dist = this_dist;
                    this_thread_best_it = i;
                }
            }
            #pragma omp critical
            {
                if (this_thread_best_dist < bestDist) {
                    bestDist = this_thread_best_dist;
                    result = this_thread_best_it;
                }
            }
        }
        return result;
        /*
         if (begin == end) [[unlikely]]
         return begin;
         std::ptrdiff_t result_idx = 0;
         auto bestDist = dist_func(*begin++);
         auto end_idx = end - begin;
         #pragma omp parallel
         {
         auto this_thread_best_dist = bestDist;
         std::ptrdiff_t this_thread_best_idx = 0;
         #pragma omp for
         for (std::ptrdiff_t i = 0; i < end_idx; ++i) {
         if (auto this_dist = dist_func(*(begin+i)); comp(this_dist, this_thread_best_dist)) {
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
         return begin + result_idx; // this version is compatible with openmp < 4.0, but requires RAI.
         */
    }
};


#endif /* Algorithm_hpp */
