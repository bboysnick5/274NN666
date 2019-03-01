/**
 * File: Point.h
 * -------------
 * A class representing a point in N-dimensional space. Unlike the other class
 * templates you've seen before, Point is parameterized over an integer rather
 * than a type. This allows the compiler to verify that the type is being used
 * correctly.
 */
#ifndef POINT_INCLUDED
#define POINT_INCLUDED

#include <cmath>
#include <stdlib.h>
#include <array>
#include <numeric>

template <typename _Tp, size_t _N>
class Point {
public:
    
    typedef Point                                 __self;
    typedef _Tp                                   value_type;
    typedef value_type&                           reference;
    typedef const value_type&                     const_reference;
    typedef value_type*                           iterator;
    typedef const value_type*                     const_iterator;
    typedef value_type*                           pointer;
    typedef const value_type*                     const_pointer;
    typedef size_t                                size_type;
    typedef ptrdiff_t                             difference_type;
    typedef std::reverse_iterator<iterator>       reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    
    enum class DistType {
        EUC = 0,
        EUCSQ,
        MAN,
        HAV
    };
    
    
    template <typename... T>
    //Point(T... ts) : coords{ts...} {}
    Point(T... ts) : coords{static_cast<value_type>(ts)...} {}
    
    Point() = default;
    
    // size_t size() const;
    // Usage: for (size_t i = 0; i < myPoint.size(); ++i)
    // ------------------------------------------------------------------------
    // Returns N, the dimension of the point.
    constexpr size_t size() const;
    
    bool operator == (const Point<_Tp, _N>& rhs) const;
    bool operator != (const Point<_Tp, _N>& rhs) const;
    
    
    // _Tp& operator[](size_t index);
    // _Tp operator[](size_t index) const;
    // Usage: myPoint[3] = 137;
    // ------------------------------------------------------------------------
    // Queries or retrieves the value of the point at a particular point. The
    // index is assumed to be in-range.
    value_type& operator[](size_t index);
    value_type operator[](size_t index) const;
    
    const value_type* data() const;
    value_type* data();
    
    const std::array<_Tp, _N>& dataArray() const;
    std::array<_Tp, _N>& dataArray();
    
    // static _Tp eulDist(const Point<_Tp, _N>& one, const Point<_Tp, _N>& two);
    // Usage: _Tp d = Distance(one, two);
    // ----------------------------------------------------------------------------
    // Returns the Euclidean distance between two points.
    template <DistType>
    static value_type dist(const Point<_Tp, _N>& pt1, const Point<_Tp, _N>& pt2);
   
    template <DistType>
    value_type dist(const Point<_Tp, _N>&) const;
    
    // iterator begin();
    // iterator end();
    // const_iterator begin() const;
    // const_iterator end() const;
    // Usage: for (Point<3>::iterator itr = myPoint.begin(); itr != myPoint.end(); ++itr)
    // ------------------------------------------------------------------------
    // Returns iterators delineating the full range of elements in the Point.
    iterator begin();
    iterator end();
    
    const_iterator cbegin() const;
    const_iterator cend() const;
    
private:
    // The point's actual coordinates are stored in an array.
    std::array<value_type, _N> coords;
    
    static value_type eucDist(const Point<_Tp, _N>& pt1, const Point<_Tp, _N>& pt2);
    static value_type eucSqDist(const Point<_Tp, _N>& pt1, const Point<_Tp, _N>& pt2);
    static value_type havDist(const Point<_Tp, _N>& pt1, const Point<_Tp, _N>& pt2);
    static value_type manDist(const Point<_Tp, _N>& pt1, const Point<_Tp, _N>& pt2);

};

/** Point class implementation details */


template <typename _Tp, size_t _N>
constexpr size_t Point<_Tp, _N>::size() const {
    return _N;
}

template <typename _Tp, size_t _N>
_Tp& Point<_Tp, _N>::operator[] (size_t index) {
    return coords[index];
}

template <typename _Tp, size_t _N>
_Tp Point<_Tp, _N>::operator[] (size_t index) const {
    return coords[index];
}

template <typename _Tp, size_t _N>
const _Tp* Point<_Tp, _N>::data() const {
    return coords.data();
}

template <typename _Tp, size_t _N>
_Tp* Point<_Tp, _N>::data() {
    return coords.data();
}

template <typename _Tp, size_t _N>
const std::array<_Tp, _N>& Point<_Tp, _N>::dataArray() const {
    return coords;
}

template <typename _Tp, size_t _N>
std::array<_Tp, _N>& Point<_Tp, _N>::dataArray() {
    return coords;
}

template <typename _Tp, size_t _N>
typename Point<_Tp, _N>::iterator Point<_Tp, _N>::begin() {
    return coords.begin();
}

template <typename _Tp, size_t _N>
typename Point<_Tp, _N>::const_iterator Point<_Tp, _N>::cbegin() const {
    return coords.cbegin();
}

template <typename _Tp, size_t _N>
typename Point<_Tp, _N>::iterator Point<_Tp, _N>::end() {
    return coords.end();
}

template <typename _Tp, size_t _N>
typename Point<_Tp, _N>::const_iterator Point<_Tp, _N>::cend() const {
    return coords.cend();
}

template <typename _Tp, size_t _N>
template <typename Point<_Tp, _N>::DistType DT>
_Tp Point<_Tp, _N>::dist(const Point<_Tp, _N>& one, const Point<_Tp, _N>& two) {
    switch (DT) {
        case DistType::EUC:
            return eucDist(one, two);
        case DistType::EUCSQ:
            return eucSqDist(one, two);
        case DistType::MAN:
            return manDist(one, two);
        case DistType::HAV:
            return havDist(one, two);
    }
}

template <typename _Tp, size_t _N>
template <typename Point<_Tp, _N>::DistType DT>
_Tp Point<_Tp, _N>::dist(const Point<_Tp, _N> &other) const {
    return Point<_Tp, _N>::template dist<DT>(*this, other);
}

template <typename _Tp, size_t _N>
_Tp Point<_Tp, _N>::eucDist(const Point<_Tp, _N>& one, const Point<_Tp, _N>& two) {
    return sqrt(eucSqDist(one, two));
}

template <typename _Tp, size_t _N>
_Tp Point<_Tp, _N>::eucSqDist(const Point<_Tp, _N>& one, const Point<_Tp, _N>& two) {
    //return std::transform_reduce(one.cbegin(), one.cend(), two.cbegin(), 0.0, std::plus<_Tp>(),
                 //                [](const _Tp a, const _Tp b){return (a-b)*(a-b);});
    
    _Tp diff = one[0] - two[0], result = diff*diff;
    for (size_t i = 1; i < _N; ++i) {
        diff = one[i] - two[i];
        result += diff*diff;
    }
    return result; 
}


template <typename _Tp, size_t _N>
_Tp Point<_Tp, _N>::havDist(const Point<_Tp, _N>& one, const Point<_Tp, _N>& two) {
    constexpr static _Tp EARTH_RADIUS = 6371.0;
    _Tp dLat = (one[1]-two[1])*M_PI/180.0;
    _Tp dLon = (one[0]-two[0])*M_PI/180.0;
    _Tp lat1 = one[1]*M_PI/180.0;
    _Tp lat2 = two[1]*M_PI/180.0;
    _Tp a = sin(dLat*0.5) * sin(dLat*0.5) +
               sin(dLon*0.5) * sin(dLon*0.5) * cos(lat1) * cos(lat2);
    _Tp c = 2.0 * atan2(sqrt(a), sqrt(1.0-a));
    return EARTH_RADIUS * c;
}

template <typename _Tp, size_t _N>
_Tp Point<_Tp, _N>::manDist(const Point<_Tp, _N>& one, const Point<_Tp, _N>& two) {
    _Tp result = std::fabs(one[0] - two[0]);
    for (size_t i = 1; i < _N; ++i)
        result += std::fabs(one[i] - two[i]);
    return result;
}

template <typename _Tp, size_t _N>
bool Point<_Tp, _N>::operator == (const Point<_Tp, _N>& rhs) const {
    return coords == rhs.coords;
}

template <typename _Tp, size_t _N>
bool Point<_Tp, _N>::operator != (const Point<_Tp, _N>& rhs) const {
    return !(*this == rhs);
}


#endif // POINT_INCLUDED
