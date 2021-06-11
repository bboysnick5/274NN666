/**
 * File: PointND.h
 * -------------
 * A class representing a point in N-dimensional space. Unlike the other class
 * templates you've seen before, PointND is parameterized over an integer rather
 * than a type. This allows the compiler to verify that the type is being used
 * correctly.
 */
#ifndef POINTND_INCLUDED
#define POINTND_INCLUDED

//#include <oneapi/dpl/execution>
#include <cmath>
#include <stdlib.h>
#include <array>
//#include <oneapi/dpl/numeric>
#include <numeric>
#include <cstddef>
#include <numbers>
#include <utility>
//#include <execution>



template <typename _Tp, std::size_t _N>
class PointND {
public:

    typedef PointND                                                        __self;
    typedef _Tp                                                          value_type;
    typedef value_type&                                                  reference;
    typedef const value_type&                                            const_reference;
    typedef typename std::array<_Tp, _N>::iterator                       iterator;
    typedef typename std::array<_Tp, _N>::const_iterator                 const_iterator;
    typedef value_type*                                                  pointer;
    typedef const value_type*                                            const_pointer;
    typedef std::size_t                                                  size_type;
    typedef ptrdiff_t                                                    difference_type;
    typedef std::reverse_iterator<iterator>                              reverse_iterator;
    typedef std::reverse_iterator<const_iterator>                        const_reverse_iterator;
    
    enum class DistType {
        EUC = 0,
        EUCSQ,
        MAN,
        HAV,
        HAVCOMP
    };
    
    
    template <typename... T>
    //PointND(T... ts) : coords{ts...} {}
    PointND(T... ts) : coords_{static_cast<value_type>(ts)...} {}
    
    PointND() = default;
    
    // std::size_t size() const;
    // Usage: for (std::size_t i = 0; i < myPointND.size(); ++i)
    // ------------------------------------------------------------------------
    // Returns N, the dimension of the point.
    constexpr std::size_t size() const;
    
    bool operator == (const PointND<_Tp, _N>& rhs) const;
    bool operator != (const PointND<_Tp, _N>& rhs) const;
    
    
    // _Tp& operator[](std::size_t index);
    // _Tp operator[](std::size_t index) const;
    // Usage: myPointND[3] = 137;
    // ------------------------------------------------------------------------
    // Queries or retrieves the value of the point at a particular point. The
    // index is assumed to be in-range.
    value_type& operator[](std::size_t index);
    value_type operator[](std::size_t index) const;
    
    const value_type* data() const;
    value_type* data();
    
    const std::array<_Tp, _N>& dataArray() const;
    std::array<_Tp, _N>& dataArray();
    
    // static _Tp eulDist(const PointND<_Tp, _N>& one, const PointND<_Tp, _N>& two);
    // Usage: _Tp d = Distance(one, two);
    // ----------------------------------------------------------------------------
    // Returns the Euclidean distance between two points.
    template <DistType, typename... VArgs>
    static value_type dist(const PointND<_Tp, _N>& pt1, const PointND<_Tp, _N>& pt2, VArgs... args);
   
    template <DistType, typename... VArgs>
    value_type dist(const PointND<_Tp, _N>&, VArgs... args) const;
    
    // iterator begin();
    // iterator end();
    // const_iterator begin() const;
    // const_iterator end() const;
    // Usage: for (PointND<3>::iterator itr = myPointND.begin(); itr != myPointND.end(); ++itr)
    // ------------------------------------------------------------------------
    // Returns iterators delineating the full range of elements in the PointND.
    iterator begin();
    iterator end();
    
    const_iterator cbegin() const;
    const_iterator cend() const;
    
private:
    // The point's actual coordinates are stored in an array.
    std::array<value_type, _N> coords_;
    
    template <std::size_t Dimension = _N, std::enable_if_t<Dimension == 2, bool> = true>
    static value_type eucDist(const PointND<_Tp, _N>& pt1, const PointND<_Tp, _N>& pt2);
    
    template <std::size_t Dimension = _N, std::enable_if_t<Dimension == 3, bool> = true>
    static value_type eucDist(const PointND<_Tp, _N>& pt1, const PointND<_Tp, _N>& pt2);
    
    template <std::size_t Dimension = _N, std::enable_if_t<Dimension != 2 && Dimension != 3, bool> = true>
    static value_type eucDist(const PointND<_Tp, _N>& pt1, const PointND<_Tp, _N>& pt2);

    static value_type eucSqDist(const PointND<_Tp, _N>& pt1, const PointND<_Tp, _N>& pt2);
    static value_type manDist(const PointND<_Tp, _N>& pt1, const PointND<_Tp, _N>& pt2);
    
    //template <typename... VArgs>
    static value_type havDist(const PointND<_Tp, _N>& pt1, const PointND<_Tp, _N>& pt2, _Tp radius);
    static value_type havDistCompCalcA(const PointND<_Tp, _N>& pt1, const PointND<_Tp, _N>& pt2);

};

/** PointND class implementation details */


template <typename _Tp, std::size_t _N>
constexpr std::size_t PointND<_Tp, _N>::size() const {
    return _N;
}

template <typename _Tp, std::size_t _N>
_Tp& PointND<_Tp, _N>::operator[] (std::size_t index) {
    return coords_[index];
}

template <typename _Tp, std::size_t _N>
_Tp PointND<_Tp, _N>::operator[] (std::size_t index) const {
    return coords_[index];
}

template <typename _Tp, std::size_t _N>
const _Tp* PointND<_Tp, _N>::data() const {
    return coords_.data();
}

template <typename _Tp, std::size_t _N>
_Tp* PointND<_Tp, _N>::data() {
    return coords_.data();
}

template <typename _Tp, std::size_t _N>
const std::array<_Tp, _N>& PointND<_Tp, _N>::dataArray() const {
    return coords_;
}

template <typename _Tp, std::size_t _N>
std::array<_Tp, _N>& PointND<_Tp, _N>::dataArray() {
    return coords_;
}

template <typename _Tp, std::size_t _N>
typename PointND<_Tp, _N>::iterator PointND<_Tp, _N>::begin() {
    return coords_.begin();
}

template <typename _Tp, std::size_t _N>
typename PointND<_Tp, _N>::const_iterator PointND<_Tp, _N>::cbegin() const {
    return coords_.cbegin();
}

template <typename _Tp, std::size_t _N>
typename PointND<_Tp, _N>::iterator PointND<_Tp, _N>::end() {
    return coords_.end();
}

template <typename _Tp, std::size_t _N>
typename PointND<_Tp, _N>::const_iterator PointND<_Tp, _N>::cend() const {
    return coords_.cend();
}

template <typename _Tp, std::size_t _N>
template <typename PointND<_Tp, _N>::DistType DT, typename... VArgs>
_Tp PointND<_Tp, _N>::dist(const PointND<_Tp, _N>& one, const PointND<_Tp, _N>& two, VArgs... args) {
    switch (DT) {
        case DistType::EUC:
            return eucDist(one, two);
        case DistType::EUCSQ:
            return eucSqDist(one, two);
        case DistType::MAN:
            return manDist(one, two);
        case DistType::HAV:
            // TODO type-safe
            return havDist(one, two, {args...});
        case DistType::HAVCOMP:
            return havDistCompCalcA(one, two);
        }
}

template <typename _Tp, std::size_t _N>
template <typename PointND<_Tp, _N>::DistType DT, typename... VArgs>
_Tp PointND<_Tp, _N>::dist(const PointND<_Tp, _N> &other, VArgs ...args) const {
    return PointND<_Tp, _N>::template dist<DT>(*this, other, args...);
}



template <typename _Tp, std::size_t _N>
template <std::size_t Dimension, std::enable_if_t<Dimension == 2, bool>>
_Tp PointND<_Tp, _N>::eucDist(const PointND<_Tp, _N>& one, const PointND<_Tp, _N>& two) {
    return std::hypot(one[0] - two[0], one[1] - two[1]);
}

template <typename _Tp, std::size_t _N>
template <std::size_t Dimension, std::enable_if_t<Dimension == 3, bool>>
_Tp PointND<_Tp, _N>::eucDist(const PointND<_Tp, _N>& one, const PointND<_Tp, _N>& two) {
    return std::hypot(one[0] - two[0], one[1] - two[1], one[2] - two[2]);
}

template <typename _Tp, std::size_t _N>
template <std::size_t Dimension, std::enable_if_t<Dimension != 2 && Dimension != 3, bool>>
_Tp PointND<_Tp, _N>::eucDist(const PointND<_Tp, _N>& one, const PointND<_Tp, _N>& two) {
    return sqrt(eucSqDist(one, two));
}

template <typename _Tp, std::size_t _N>
_Tp PointND<_Tp, _N>::eucSqDist(const PointND<_Tp, _N>& one, const PointND<_Tp, _N>& two) {
    return std::transform_reduce(one.cbegin(), one.cend(), two.cbegin(), 0.0, std::plus<_Tp>(),
                                 [](const _Tp a, const _Tp b){return (a-b)*(a-b);});
    /*
    _Tp diff = one[0] - two[0], result = diff*diff;
    for (std::size_t i = 1; i != _N; ++i) {
        diff = one[i] - two[i];
        result += diff*diff;
    }
    return result;  */
}


template <typename _Tp, std::size_t _N>
_Tp PointND<_Tp, _N>::havDist(const PointND<_Tp, _N>& one, const PointND<_Tp, _N>& two, _Tp radius) {
    return 2.0 * std::asin(std::sqrt(havDistCompCalcA(one, two))) * radius;
}

template <typename _Tp, std::size_t _N>
_Tp PointND<_Tp, _N>::havDistCompCalcA(const PointND<_Tp, _N>& one, const PointND<_Tp, _N>& two) {
    _Tp dLatOneHalf = 0.5*(one[0] - two[0]);
    _Tp dLonOneHalf = 0.5*(one[1] - two[1]);
    return sin(dLatOneHalf) * sin(dLatOneHalf) +
           sin(dLonOneHalf) * sin(dLonOneHalf) * std::cos(one[0]) * std::cos(two[0]);
}


template <typename _Tp, std::size_t _N>
_Tp PointND<_Tp, _N>::manDist(const PointND<_Tp, _N>& one, const PointND<_Tp, _N>& two) {
    _Tp result = std::fabs(one[0] - two[0]);
    for (std::size_t i = 1; i < _N; ++i)
        result += std::fabs(one[i] - two[i]);
    return result;
}

template <typename _Tp, std::size_t _N>
bool PointND<_Tp, _N>::operator == (const PointND<_Tp, _N>& rhs) const {
    return coords_ == rhs.coords_;
}

template <typename _Tp, std::size_t _N>
bool PointND<_Tp, _N>::operator != (const PointND<_Tp, _N>& rhs) const {
    return !(*this == rhs);
}


#endif // POINTND_INCLUDED
