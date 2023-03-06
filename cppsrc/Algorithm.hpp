/*
 * @Author: Nick Liu
 * @Date: 2021-07-23 16:02:39
 * @LastEditTime: 2022-08-08 20:46:51
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/Algorithm.hpp
 */

#ifndef Algorithm_hpp
#define Algorithm_hpp

#include "Definition.hpp"
#include "Point.hpp"
#include "SBLoc.hpp"
#include "Utility.hpp"

// #include <Vc/Vc>

#include <stdio.h>

#include <concepts>
#include <iterator>
// #include <execution>
#include <omp.h>

#include <type_traits>
#include <utility>

namespace algo {

/*
template <def::ParPolicy ParPolicy, std::floating_point FPType,
          def::DistType DT = def::DistType::kEucSq,
          std::forward_iterator SBlocIt>
requires(std::same_as<typename std::iterator_traits<SBlocIt>::value_type,
                      SBLoc<FPType>>) SBlocIt
    LinearNNSearch(SBlocIt begin, SBlocIt end, const auto &pt) {
    if (begin != end) [[likely]] {
        decltype(begin) i = begin;
        auto bestDist = pt.template dist<DT>(*begin);
        while (++i != end) {
            if (const auto dist = pt.template dist<DT>(*i); dist < bestDist) {
                bestDist = dist;
                begin = i;
            }
        }
    }
    return begin;
}

template <def::ParPolicy ParPolicy, def::DistType DT = def::DistType::kEucSq>
requires(Config.thread_policy == def::ThreadingPolicy::kSingle &&
         Config.mem_layout == def::MemoryLayout::kSoA) std::forward_iterator
    auto LinearNNSearch(std::forward_iterator auto begin,
                        std::forward_iterator auto end, const auto &pt) {
    if (begin != end) [[likely]] {
        decltype(begin) i = begin;
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
*/
template <std::forward_iterator PointIt>
struct Deref {
    typename std::iterator_traits<PointIt>::value_type operator()(PointIt it) {
        return *it;
    }
};

template <def::ParPolicy ParPolicy, def::DistType DT = def::DistType::kEucSq,
          std::forward_iterator ContainerIt, class ToPtFn = Deref<ContainerIt>,
          class Point>
    requires requires(ContainerIt it, ToPtFn to_pt_fn) {
                 { to_pt_fn(it) } -> std::convertible_to<Point>;
             }
inline ContainerIt LinearNNSearch(ContainerIt begin, ContainerIt end,
                                  const Point &pt,
                                  ToPtFn &&to_pt_fn = Deref<ContainerIt>{}) {
    if (begin != end) [[likely]] {
        decltype(begin) i = begin;
        auto bestDist = pt.template dist<DT>(to_pt_fn(begin));
        while (++i != end) {
            if (const auto dist = pt.template dist<DT>(to_pt_fn(i));
                dist < bestDist) {
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
std::enable_if_t<Policy == def::ThreadingPolicy::kSimdAosoa && DT ==
def::DistType::kEucSq, bool> = true> inline static SimdizedRAI
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
/*
template <class RAI, class PTType, def::ThreadingPolicy P = Policy,
          std::enable_if_t<P == def::ThreadingPolicy::kSimdSoA &&
                               DT == def::DistType::kEucSq &&
                               std::is_same_v<typename std::iterator_traits<
                                                  RAI>::iterator_category,
                                              std::random_access_iterator_tag>,
                           bool> = true>
inline static RAI LinearNNSearch(RAI begin, RAI end, const PTType &pt) {
    if (begin != end) [[likely]] {
        Vc::simdize<RAI> i = begin;
        auto bestDist = pt.template dist<DT>(*i);
        ;
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
                               std::is_same_v<typename std::iterator_traits<
                                                  RAI>::iterator_category,
                                              std::random_access_iterator_tag>,
                           bool> = true>
inline static RAI LinearNNSearch(RAI begin, RAI end, const PTType &pt) {
    if (begin != end) [[likely]] {
        Vc::simdize<RAI> i = begin;
        auto bestDist = pt.template dist<DT>(*i);
        ;
        while (++i != end) {
            if (auto dist = pt.template dist<DT>(*i); dist < bestDist) {
                bestDist = dist;
                begin = i;
            }
        }
    }
    return begin;
} */

template <def::ParPolicy ParPolicy, def::DistType DT = def::DistType::kEucSq>
    requires(ParPolicy.thread_policy == def::ThreadingPolicy::kMultiOmp &&
             ParPolicy.mem_layout == def::MemoryLayout::kSoA)
std::forward_iterator
    auto LinearNNSearch(std::forward_iterator auto begin,
                        std::forward_iterator auto end, const auto &pt) {
    if (begin == end) [[unlikely]]
        return begin;
    auto result = begin;
    auto bestDist = pt.template dist<DT>(*begin++);
#pragma omp parallel
    {
        auto this_thread_best_dist = bestDist;
        decltype(begin) this_thread_best_it = begin;
#pragma omp for
        for (decltype(begin) i = begin; i < end; ++i) {
            if (auto this_dist = pt.template dist<DT>(*i);
                this_dist < this_thread_best_dist) {
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
     auto bestDist = dist_fn(*begin++);
     auto end_idx = end - begin;
     #pragma omp parallel
     {
     auto this_thread_best_dist = bestDist;
     std::ptrdiff_t this_thread_best_idx = 0;
     #pragma omp for
     for (std::ptrdiff_t i = 0; i < end_idx; ++i) {
     if (auto this_dist = dist_fn(*(begin+i)); comp(this_dist,
     this_thread_best_dist)) { this_thread_best_dist = this_dist;
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
     return begin + result_idx; // this version is compatible with openmp
     < 4.0, but requires RAI.
     */
}

template <def::ParPolicy ParPolicy, def::DistType DT = def::DistType::kEucSq>
    requires(ParPolicy.thread_policy == def::ThreadingPolicy::kMultiHand &&
             ParPolicy.mem_layout == def::MemoryLayout::kSoA)
std::forward_iterator
    auto LinearNNSearch(std::forward_iterator auto begin,
                        std::forward_iterator auto end, const auto &pt) {
    if (begin == end) [[unlikely]]
        return begin;
    auto result = begin;
    auto bestDist = pt.template dist<DT>(*begin++);
#pragma omp parallel
    {
        auto this_thread_best_dist = bestDist;
        decltype(begin) this_thread_best_it = begin;
#pragma omp for
        for (decltype(begin) i = begin; i < end; ++i) {
            if (auto this_dist = pt.template dist<DT>(*i);
                this_dist < this_thread_best_dist) {
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
     auto bestDist = dist_fn(*begin++);
     auto end_idx = end - begin;
     #pragma omp parallel
     {
     auto this_thread_best_dist = bestDist;
     std::ptrdiff_t this_thread_best_idx = 0;
     #pragma omp for
     for (std::ptrdiff_t i = 0; i < end_idx; ++i) {
     if (auto this_dist = dist_fn(*(begin+i)); comp(this_dist,
     this_thread_best_dist)) { this_thread_best_dist = this_dist;
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
     return begin + result_idx; // this version is compatible with openmp
     < 4.0, but requires RAI.
     */
}

};  // namespace algo

#endif /* Algorithm_hpp */
