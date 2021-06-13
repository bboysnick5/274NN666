//
//  Definition.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 5/25/21.
//  Copyright © 2021 Yunlong Liu. All rights reserved.
//

#ifndef Definition_hpp
#define Definition_hpp

#include <stdio.h>
#include <numbers>
#include <iterator>

namespace def {

using kDefaultDistType = double;
inline constexpr std::size_t kDefaultMaxCacheCellVecSize = 896;
// to better fit into 32kb L1 with each vec element being 3*sizeof(double) + pointer

inline constexpr double kDefalutAveActualLocsPerCell = 0.4;

inline constexpr bool kDefaultToTestAccuracy = false;
inline constexpr std::size_t kDefaultAccuracyTestDurationInSecs = 10;
inline constexpr std::size_t kDefaultSearchBenchDurationInSecs = 10;
inline constexpr std::size_t kMaxTestLocs = 1 << 24;

template <class ConstIt>
concept const_iterator = std::is_const_v<typename std::remove_pointer_t<typename std::iterator_traits<ConstIt>::pointer>>;

template <class Non_Const_It>
concept non_const_iterator = !std::is_const_v<typename std::remove_pointer_t<typename std::iterator_traits<Non_Const_It>::pointer>>;



template <typename T>
constexpr auto &kMathPi = std::numbers::pi_v<T>;

enum class ThreadingPolicy {
    kSingle,
    kMultiOmp,
    kMultiHand
};

template <def::ThreadingPolicy = def::ThreadingPolicy::kSingle>
struct Policy_Tag {};
    
template <>
struct Policy_Tag<def::ThreadingPolicy::kMultiOmp> {};
    
template <>
struct Policy_Tag<def::ThreadingPolicy::kMultiHand> {};

}

namespace debug {
inline static bool mis_match = false;
}


#endif /* Definition_hpp */
