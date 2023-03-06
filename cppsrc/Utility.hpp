/*
 * @Author: Nick Liu
 * @Date: 2021-06-21 10:18:59
 * @LastEditTime: 2022-09-05 14:41:41
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/Utility.hpp
 */

#ifndef Utility_hpp
#define Utility_hpp

#include "Definition.hpp"

// #include <Vc/Vc>

#include <stdio.h>
// #include <execution>
#include <omp.h>

#include <boost/mp11/algorithm.hpp>
#include <type_traits>
#include <utility>

/*
 Allow pass-in value calculation formula and store the computed data
 to save up one value computation as it were two times
 as if provided via comparator.
 */

namespace utility {

using namespace boost::mp11;

template <typename Variant, typename T>
constexpr size_t IndexInVariant = mp_find<Variant, T>::value;

// cycle swap
template <typename T>
constexpr void CycleSwap(T &t1, T &t2, T &t3) {
    T temp(std::move(t1));
    t1 = std::move(t2);
    t2 = std::move(t3);
    t3 = std::move(temp);
}

struct NoOp {
    template <class... Args>
    constexpr void operator()(Args &&...) const {}
};

template <class Compare, class ForwardIterator>
void Sort3(ForwardIterator x, ForwardIterator y, ForwardIterator z, Compare c) {
    if (!c(*y, *x))  // if x <= y
    {
        if (!c(*z, *y))         // if y <= z
            return;             // x <= y && y <= z
                                // x <= y && y > z
        std::swap(*y, *z);      // x <= z && y < z
        if (c(*y, *x))          // if x > y
            std::swap(*x, *y);  // x < y && y <= z
        return;
    }
    if (c(*z, *y)) {        // x > y, if y > z
        std::swap(*x, *z);  // x < y && y < z
        return;
    }
    std::swap(*x, *y);      // x > y && y <= z
                            // x < y && x <= z
    if (c(*z, *y))          // if y > z
        std::swap(*y, *z);  // x <= y && y < z
}

template <class T>
struct movable_il {
    mutable T t;
    constexpr operator T() const && { return std::move(t); }
    constexpr movable_il(T &&in) : t(std::move(in)) {}
};

template <class VT>
struct fix_vt {
    using type = VT;
};
template <class VT>
using fix_vt_t = typename fix_vt<VT>::type;
template <class VT>
struct fix_vt<const VT> : fix_vt<VT> {};
template <class K, class V>
struct fix_vt<std::pair<K, V>> {
    using type = std::pair<std::remove_cv_t<K>, std::remove_cv_t<V>>;
};

template <class C, class T = fix_vt_t<typename C::value_type>>
constexpr C container_from_il(std::initializer_list<movable_il<T>> il) {
    C r(std::make_move_iterator(il.begin()), std::make_move_iterator(il.end()));
    return r;
}

};     // namespace utility
#endif /* Utility_hpp */
