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

namespace def {

using kDefaultDistType = double;
inline constexpr std::size_t kDefaultMaxCacheCellVecSize = 896;
// to better fit into 32kb L1 with each vec element being 3*sizeof(double) + pointer

inline constexpr double kDefalutAveActualLocsPerCell = 0.4;

inline constexpr bool kDefaultToTestAccuracy = false;
inline constexpr std::size_t kDefaultAccuracyTestDurationInSecs = 10;
inline constexpr std::size_t kDefaultSearchBenchDurationInSecs = 10;
inline constexpr std::size_t kMaxTestLocs = 1 << 24;



template <typename T>
constexpr auto &kMathPi = std::numbers::pi_v<T>;

enum class ThreadingPolicy {
    kSingle,
    kMultiOmp,
    kMultiHand
};
}

namespace debug {
inline static bool mis_match = false;
}


#endif /* Definition_hpp */
