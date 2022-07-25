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

#include <Vc/Vc>

#include <stdio.h>
//#include <execution>
#include <omp.h>
#include <utility>
#include <type_traits>



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

};
#endif /* Utility_hpp */
