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

template <class _ForwardIterator, class _GetDist, class _Compare>
_ForwardIterator
custom_min_element(_ForwardIterator __first, _ForwardIterator __last, _GetDist __gd, _Compare __comp) {
    //if (__first == __last)
      //  return __first;
    
    _ForwardIterator __result = __first;
    auto __bestDist = __gd(*++__first);
    for(; __first != __last; ++__first) {
        if (auto __dist = __gd(*__first); __comp(__dist, __bestDist)) {
            __bestDist = __dist;
            __result = __first;
        }
    }
    return __result;
}

#endif /* Utility_hpp */
