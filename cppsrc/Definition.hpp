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
    
    template <typename T>
    constexpr auto &PI = std::numbers::pi_v<T>;

    enum class Threading_Policy {
        SINGLE,
        MULTI_OMP,
        MULTI_HAND
    };
}

#endif /* Definition_hpp */
