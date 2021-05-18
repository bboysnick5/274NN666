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

/*
 Allow pass-in value calculation formula and store the computed data
 to save up one value computation as it were two times
 as if provided via comparator. 
 */
template <class _ForwardIterator, class _GetDist, class _Compare>
_ForwardIterator
custom_min_element(_ForwardIterator __first, _ForwardIterator __last, _GetDist __gd, _Compare __comp) {
    if (__first != __last) {
        _ForwardIterator __i = __first;
        auto __bestDist = __gd(*__first);
        while (++__i != __last) {
            if (auto __dist = __gd(*__i); __comp(__dist, __bestDist)) {
                __bestDist = __dist;
                __first = __i;
            }
        }
    }
    return __first;
}

#endif /* Utility_hpp */
