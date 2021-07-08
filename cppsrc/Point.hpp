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
#include <iterator>
#include <utility>
//#include <execution>

namespace ns {
    struct Point {};
    struct X {};
    struct Y {};
    struct Z {};
}

template <typename FPType, std::uint8_t N>
class PointND {
public:



    typedef PointND                                                        __self;
    typedef FPType                                                          value_type;
    typedef value_type&                                                  reference;
    typedef const value_type&                                            const_reference;
    typedef typename std::array<FPType, N>::iterator                       iterator;
    typedef typename std::array<FPType, N>::const_iterator                 const_iterator;
    typedef value_type*                                                  pointer;
    typedef const value_type*                                            const_pointer;
    typedef std::uint8_t                                                  size_type;
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
    
    // std::uint8_t size() const;
    // Usage: for (std::uint8_t i = 0; i < myPointND.size(); ++i)
    // ------------------------------------------------------------------------
    // Returns N, the dimension of the point.
    constexpr std::uint8_t size() const;
    
    bool operator == (const PointND<FPType, N>& rhs) const;
    bool operator != (const PointND<FPType, N>& rhs) const;
    
    
    // FPType& operator[](std::uint8_t dim);
    // FPType operator[](std::uint8_t dim) const;
    // Usage: myPointND[3] = 137;
    // ------------------------------------------------------------------------
    // Queries or retrieves the value of the point at a particular point. The
    // dim is assumed to be in-range.
    value_type& operator[](std::uint8_t dim);
    value_type operator[](std::uint8_t dim) const;
    
    const value_type* data() const;
    value_type* data();
    
    const std::array<FPType, N>& DataArray() const;
    std::array<FPType, N>& DataArray();
    
    // static FPType eulDist(const PointND<FPType, N>& one, const PointND<FPType, N>& two);
    // Usage: FPType d = Distance(one, two);
    // ----------------------------------------------------------------------------
    // Returns the Euclidean distance between two points.
    template <DistType, typename... VArgs>
    static value_type dist(const PointND<FPType, N>& pt1, const PointND<FPType, N>& pt2, VArgs... args);
   
    template <DistType, typename... VArgs>
    value_type dist(const PointND<FPType, N>&, VArgs... args) const;
    
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
    std::array<value_type, N> coords_;
    
    template <std::uint8_t Dim = N, std::enable_if_t<Dim == 2, bool> = true>
    static value_type EucDist(const PointND<FPType, N>& pt1, const PointND<FPType, N>& pt2);
    
    template <std::uint8_t Dim = N, std::enable_if_t<Dim == 3, bool> = true>
    static value_type EucDist(const PointND<FPType, N>& pt1, const PointND<FPType, N>& pt2);
    
    template <std::uint8_t Dim = N, std::enable_if_t<Dim != 2 && Dim != 3, bool> = true>
    static value_type EucDist(const PointND<FPType, N>& pt1, const PointND<FPType, N>& pt2);

    static value_type EucSqDist(const PointND<FPType, N>& pt1, const PointND<FPType, N>& pt2);
    static value_type ManDist(const PointND<FPType, N>& pt1, const PointND<FPType, N>& pt2);
    
    //template <typename... VArgs>
    static value_type HavDist(const PointND<FPType, N>& pt1, const PointND<FPType, N>& pt2, FPType radius);
    static value_type HavDistCompCalcA(const PointND<FPType, N>& pt1, const PointND<FPType, N>& pt2);

};

/** PointND class implementation details */


template <typename FPType, std::uint8_t N>
constexpr std::uint8_t PointND<FPType, N>::size() const {
    return N;
}

template <typename FPType, std::uint8_t N>
FPType& PointND<FPType, N>::operator[] (std::uint8_t dim) {
    return coords_[dim];
}

template <typename FPType, std::uint8_t N>
FPType PointND<FPType, N>::operator[] (std::uint8_t dim) const {
    return coords_[dim];
}

template <typename FPType, std::uint8_t N>
const FPType* PointND<FPType, N>::data() const {
    return coords_.data();
}

template <typename FPType, std::uint8_t N>
FPType* PointND<FPType, N>::data() {
    return coords_.data();
}

template <typename FPType, std::uint8_t N>
const std::array<FPType, N>& PointND<FPType, N>::DataArray() const {
    return coords_;
}

template <typename FPType, std::uint8_t N>
std::array<FPType, N>& PointND<FPType, N>::DataArray() {
    return coords_;
}

template <typename FPType, std::uint8_t N>
typename PointND<FPType, N>::iterator PointND<FPType, N>::begin() {
    return coords_.begin();
}

template <typename FPType, std::uint8_t N>
typename PointND<FPType, N>::const_iterator PointND<FPType, N>::cbegin() const {
    return coords_.cbegin();
}

template <typename FPType, std::uint8_t N>
typename PointND<FPType, N>::iterator PointND<FPType, N>::end() {
    return coords_.end();
}

template <typename FPType, std::uint8_t N>
typename PointND<FPType, N>::const_iterator PointND<FPType, N>::cend() const {
    return coords_.cend();
}

template <typename FPType, std::uint8_t N>
template <typename PointND<FPType, N>::DistType DT, typename... VArgs>
FPType PointND<FPType, N>::dist(const PointND<FPType, N>& one, const PointND<FPType, N>& two, VArgs... args) {
    switch (DT) {
        case DistType::EUC:
            return EucDist(one, two);
        case DistType::EUCSQ:
            return EucSqDist(one, two);
        case DistType::MAN:
            return ManDist(one, two);
        case DistType::HAV:
            // TODO type-safe
            return HavDist(one, two, {args...});
        case DistType::HAVCOMP:
            return HavDistCompCalcA(one, two);
        }
}

template <typename FPType, std::uint8_t N>
template <typename PointND<FPType, N>::DistType DT, typename... VArgs>
FPType PointND<FPType, N>::dist(const PointND<FPType, N> &other, VArgs ...args) const {
    return PointND<FPType, N>::template dist<DT>(*this, other, args...);
}



template <typename FPType, std::uint8_t N>
template <std::uint8_t Dim, std::enable_if_t<Dim == 2, bool>>
FPType PointND<FPType, N>::EucDist(const PointND<FPType, N>& one, const PointND<FPType, N>& two) {
    return std::hypot(one[0] - two[0], one[1] - two[1]);
}

template <typename FPType, std::uint8_t N>
template <std::uint8_t Dim, std::enable_if_t<Dim == 3, bool>>
FPType PointND<FPType, N>::EucDist(const PointND<FPType, N>& one, const PointND<FPType, N>& two) {
    return std::hypot(one[0] - two[0], one[1] - two[1], one[2] - two[2]);
}

template <typename FPType, std::uint8_t N>
template <std::uint8_t Dim, std::enable_if_t<Dim != 2 && Dim != 3, bool>>
FPType PointND<FPType, N>::EucDist(const PointND<FPType, N>& one, const PointND<FPType, N>& two) {
    return sqrt(EucSqDist(one, two));
}

template <typename FPType, std::uint8_t N>
FPType PointND<FPType, N>::EucSqDist(const PointND<FPType, N>& one, const PointND<FPType, N>& two) {
    return std::transform_reduce(one.cbegin(), one.cend(), two.cbegin(), FPType(0), std::plus<FPType>(),
                                 [](FPType a, FPType b){return (a-b)*(a-b);});
    /*
    FPType diff = one[0] - two[0], result = diff*diff;
    for (std::uint8_t i = 1; i != N; ++i) {
        diff = one[i] - two[i];
        result += diff*diff;
    }
    return result;  */
}


template <typename FPType, std::uint8_t N>
FPType PointND<FPType, N>::HavDist(const PointND<FPType, N>& one, const PointND<FPType, N>& two, FPType radius) {
    return FPType(2.0) * std::asin(std::sqrt(HavDistCompCalcA(one, two))) * radius;
}

template <typename FPType, std::uint8_t N>
FPType PointND<FPType, N>::HavDistCompCalcA(const PointND<FPType, N>& one, const PointND<FPType, N>& two) {
    FPType dLatOneHalf = FPType(0.5)*(one[0] - two[0]);
    FPType dLonOneHalf = FPType(0.5)*(one[1] - two[1]);
    return sin(dLatOneHalf) * sin(dLatOneHalf) +
           sin(dLonOneHalf) * sin(dLonOneHalf) * std::cos(one[0]) * std::cos(two[0]);
}


template <typename FPType, std::uint8_t N>
FPType PointND<FPType, N>::ManDist(const PointND<FPType, N>& one, const PointND<FPType, N>& two) {
    FPType result = std::fabs(one[0] - two[0]);
    for (std::uint8_t i = 1; i < N; ++i)
        result += std::fabs(one[i] - two[i]);
    return result;
}

template <typename FPType, std::uint8_t N>
bool PointND<FPType, N>::operator == (const PointND<FPType, N>& rhs) const {
    return coords_ == rhs.coords_;
}

template <typename FPType, std::uint8_t N>
bool PointND<FPType, N>::operator != (const PointND<FPType, N>& rhs) const {
    return !(*this == rhs);
}


#endif // POINTND_INCLUDED
