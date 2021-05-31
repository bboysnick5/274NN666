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

namespace Def {

using DEFAULT_DIST_TYPE = double;
inline constexpr std::size_t DEFAULT_MAX_CACHE_CELL_VEC_SIZE = 1200;

inline constexpr double DEFAULT_AVE_ACTUAL_LOCS_PER_CELL = 0.4;

inline constexpr bool DEFAULT_TEST_ACCURACY = false;
inline constexpr std::size_t DEFAULT_TEST_ACCURACY_DURATION_IN_SECONDS = 10;
inline constexpr std::size_t DEFAULT_TIME_SEARCH_DURATION_IN_SECONDS = 10;
inline constexpr std::size_t MAX_SEARCH_LOCS = 1 << 24;



template <typename T>
constexpr auto &PI = std::numbers::pi_v<T>;

enum class Threading_Policy {
    SINGLE,
    MULTI_OMP,
    MULTI_HAND
};
}

#endif /* Definition_hpp */
