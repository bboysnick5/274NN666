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
//#include <initializer_list>

template <size_t N>
class Point {
public:
    
    enum class DistType {
        EUC = 0,
        EUCSQ,
        MAN,
        HAV
    };
    
    template <typename... T>
    Point(T... ts) : coords{ts...} {}
    
    Point() = default;
    
    // Type: iterator
    // Type: const_iterator
    // ------------------------------------------------------------------------
    // Types representing iterators that can traverse and optionally modify the
    // elements of the Point.
    typedef double* iterator;
    typedef const double* const_iterator;
    
    // size_t size() const;
    // Usage: for (size_t i = 0; i < myPoint.size(); ++i)
    // ------------------------------------------------------------------------
    // Returns N, the dimension of the point.
    size_t size() const;
    
    // double& operator[](size_t index);
    // double operator[](size_t index) const;
    // Usage: myPoint[3] = 137;
    // ------------------------------------------------------------------------
    // Queries or retrieves the value of the point at a particular point. The
    // index is assumed to be in-range.
    double& operator[](size_t index);
    double operator[](size_t index) const;
    
    
    // static double eulDist(const Point<N>& one, const Point<N>& two);
    // Usage: double d = Distance(one, two);
    // ----------------------------------------------------------------------------
    // Returns the Euclidean distance between two points.
    template <DistType>
    static double dist(const Point<N>& pt1, const Point<N>& pt2);
    
    // iterator begin();
    // iterator end();
    // const_iterator begin() const;
    // const_iterator end() const;
    // Usage: for (Point<3>::iterator itr = myPoint.begin(); itr != myPoint.end(); ++itr)
    // ------------------------------------------------------------------------
    // Returns iterators delineating the full range of elements in the Point.
    iterator begin();
    iterator end();
    
    const_iterator begin() const;
    const_iterator end() const;
    
private:
    // The point's actual coordinates are stored in an array.
    double coords[N];
    
    static double eucDist(const Point<N>& pt1, const Point<N>& pt2);
    static double eucSqDist(const Point<N>& pt1, const Point<N>& pt2);
    static double havDist(const Point<N>& pt1, const Point<N>& pt2);
    static double manDist(const Point<N>& pt1, const Point<N>& pt2);

};


// bool operator==(const Point<N>& one, const Point<N>& two);
// bool operator!=(const Point<N>& one, const Point<N>& two);
// Usage: if (one == two)
// ----------------------------------------------------------------------------
// Returns whether two points are equal or not equal.
template <size_t N>
bool operator==(const Point<N>& one, const Point<N>& two);

template <size_t N>
bool operator!=(const Point<N>& one, const Point<N>& two);

/** Point class implementation details */

#include <algorithm>

template <size_t N>
size_t Point<N>::size() const {
    return N;
}

template <size_t N>
double& Point<N>::operator[] (size_t index) {
    return *(coords + index);
}

template <size_t N>
double Point<N>::operator[] (size_t index) const {
    return *(coords + index);
}

template <size_t N>
typename Point<N>::iterator Point<N>::begin() {
    return coords;
}

template <size_t N>
typename Point<N>::const_iterator Point<N>::begin() const {
    return coords;
}

template <size_t N>
typename Point<N>::iterator Point<N>::end() {
    return begin() + size();
}

template <size_t N>
typename Point<N>::const_iterator Point<N>::end() const {
    return begin() + size();
}

template <size_t N>
template <typename Point<N>::DistType DT>
double Point<N>::dist(const Point<N>& one, const Point<N>& two) {
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

template <size_t N>
double Point<N>::eucDist(const Point<N>& one, const Point<N>& two) {
    return sqrt(eucSqDist(one, two));
}

template <size_t N>
double Point<N>::eucSqDist(const Point<N>& one, const Point<N>& two) {
    double diff = *one.coords - *two.coords, result = diff*diff;
    for (size_t i = 1; i < N; ++i) {
        diff = *(one.coords+i) - *(two.coords+i);
        result += diff*diff;
    }
    return result;
}


template <size_t N>
double Point<N>::havDist(const Point<N>& one, const Point<N>& two) {
    constexpr static double EARTH_RADIUS = 6371;
    double dLat = (one[1]-two[1])*M_PI/180;
    double dLon = (one[0]-two[0])*M_PI/180;
    double lat1 = one[1]*M_PI/180;
    double lat2 = two[1]*M_PI/180;
    double a = sin(dLat/2) * sin(dLat/2) +
               sin(dLon/2) * sin(dLon/2) * cos(lat1) * cos(lat2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return EARTH_RADIUS * c;
}

template <size_t N>
double Point<N>::manDist(const Point<N>& one, const Point<N>& two) {
    double result = std::fabs(one[0] - two[0]);
    for (size_t i = 1; i < N; ++i)
        result += std::fabs(one[i] - two[i]);
    return result;
}

// Equality is implemented using the equal algorithm, which takes in two ranges
// and reports whether they contain equal values.
template <size_t N>
bool operator==(const Point<N>& one, const Point<N>& two) {
    return std::equal(one.begin(), one.end(), two.begin());
}

template <size_t N>
bool operator!=(const Point<N>& one, const Point<N>& two) {
    return !(one == two);
}


#endif // POINT_INCLUDED
