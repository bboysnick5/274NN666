//
//  SBLoc<FPType>.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/13/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef SBLoc_hpp
#define SBLoc_hpp

#include "Definition.hpp"
#include "Point.hpp"
#include <iostream>
#include <stdio.h>
#include <string>

template <typename FPType> struct SBLoc {

    using GeoPtType = PointND<FPType, 2>;
    using CartPtType = PointND<FPType, 3>;

    typename SBLoc<FPType>::GeoPtType geo_pt;
    std::string city;
    std::string addr;

    static constexpr FPType EARTH_RADIUS = FPType(6371.0);

    SBLoc() = default;

    bool operator<(const SBLoc<FPType> &rhs) const;

    bool operator==(const SBLoc<FPType> &rhs) const;

    bool operator!=(const SBLoc<FPType> &rhs) const;

    static FPType toRadians(FPType degree);

    static FPType toDegree(FPType radians);

    static FPType havDist(const SBLoc<FPType> &, const SBLoc<FPType> &);

    static FPType havDist(const GeoPtType &p1, const GeoPtType &p2);

    FPType havDist(const SBLoc<FPType> &other) const;

    FPType havDist(const GeoPtType &other) const;

    FPType havDistComp(const GeoPtType &other) const;

    static FPType deltaLatOnSameLngFromHavDist(FPType dist);

    static FPType lngFromSameLatHavDist(FPType dist, FPType lng1, FPType lat);

    CartPtType LocTo3dCartPt() const;

    static CartPtType GeoPtToCartPt(const GeoPtType &);

    static FPType CART3DDistFromLatDeltaLng(FPType lat1, FPType lat2,
                                            FPType deltaLng);

    static FPType CART3DDistSqFromLatDeltaLng(FPType lat1, FPType lat2,
                                              FPType deltaLng);
};

template <typename FPType>
inline bool SBLoc<FPType>::operator<(const SBLoc<FPType> &rhs) const {
    return geo_pt[0] < rhs.geo_pt[0];
}

template <typename FPType>
inline bool SBLoc<FPType>::operator==(const SBLoc<FPType> &rhs) const {
    if (geo_pt == rhs.geo_pt)
        return true;
    return std::fabs(geo_pt[0] - rhs.geo_pt[0]) <
               static_cast<FPType>(0.000001) &&
           std::fabs(geo_pt[1] - rhs.geo_pt[1]) < static_cast<FPType>(0.000001);
}

template <typename FPType>
inline bool SBLoc<FPType>::operator!=(const SBLoc<FPType> &rhs) const {
    return !(*this == rhs);
}

template <typename FPType>
inline FPType SBLoc<FPType>::toDegree(FPType radians) {
    return radians * 180.0 / def::kMathPi<FPType>;
}

template <typename FPType>
inline FPType SBLoc<FPType>::toRadians(FPType degree) {
    return degree * def::kMathPi<FPType> / 180.0;
}

template <typename FPType>
FPType SBLoc<FPType>::havDist(const SBLoc<FPType> &other) const {
    return havDist(geo_pt, other.geo_pt);
}

template <typename FPType>
FPType SBLoc<FPType>::havDist(const GeoPtType &otherGeoPt) const {
    return GeoPtType::template dist<def::DistType::kHav>(geo_pt, otherGeoPt,
                                                         EARTH_RADIUS);
}

template <typename FPType>
FPType SBLoc<FPType>::havDistComp(const GeoPtType &otherGeoPt) const {
    return GeoPtType::template dist<def::DistType::kHavComp>(geo_pt, otherGeoPt,
                                                             EARTH_RADIUS);
}

template <typename FPType>
inline FPType SBLoc<FPType>::havDist(const SBLoc<FPType> &l1,
                                     const SBLoc<FPType> &l2) {
    return GeoPtType::template dist<def::DistType::kHav>(l1.geo_pt, l2.geo_pt,
                                                         EARTH_RADIUS);
}

template <typename FPType>
inline FPType SBLoc<FPType>::havDist(const GeoPtType &p1, const GeoPtType &p2) {
    return GeoPtType::template dist<def::DistType::kHav>(p1, p2, EARTH_RADIUS);
}

template <typename FPType>
inline FPType SBLoc<FPType>::deltaLatOnSameLngFromHavDist(FPType dist) {
    FPType c = dist / EARTH_RADIUS;
    FPType sum = sin(c * 0.5);
    FPType dLat = std::asin(sum) * 2.0;
    return dLat;
}

template <typename FPType>
inline FPType SBLoc<FPType>::lngFromSameLatHavDist(FPType dist, FPType lng1,
                                                   FPType lat) {
    FPType c = dist / EARTH_RADIUS;
    FPType sum = sin(c * 0.5);
    FPType dLng = std::asin(sum / std::cos(lat)) * 2.0;
    return dLng + lng1;
}

template <typename FPType>
inline typename SBLoc<FPType>::CartPtType SBLoc<FPType>::LocTo3dCartPt() const {
    return GeoPtToCartPt(geo_pt);
}

template <typename FPType>
inline typename SBLoc<FPType>::CartPtType
SBLoc<FPType>::GeoPtToCartPt(const GeoPtType &geo_pt) {
    return CartPtType{std::cos(geo_pt[0]) * std::cos(geo_pt[1]),
                      std::cos(geo_pt[0]) * sin(geo_pt[1]), sin(geo_pt[0])};
}

template <typename FPType>
inline FPType SBLoc<FPType>::CART3DDistFromLatDeltaLng(FPType lat1, FPType lat2,
                                                       FPType deltaLng) {
    return sqrt(CART3DDistSqFromLatDeltaLng(lat1, lat2, deltaLng));
}

template <typename FPType>
inline FPType SBLoc<FPType>::CART3DDistSqFromLatDeltaLng(FPType lat1,
                                                         FPType lat2,
                                                         FPType deltaLng) {
    return 2.0 - 2.0 * (sin(lat1) * sin(lat2) +
                        std::cos(lat1) * std::cos(lat2) * std::cos(deltaLng));
}

namespace std {
template <typename FPType> struct hash<SBLoc<FPType>> {
    std::size_t operator()(const SBLoc<FPType> &l) const {
        return (static_cast<std::size_t>(
                    (l.geo_pt[0] + 0.5 * def::kMathPi<FPType>)*1000000.0)
                << 20) +
               static_cast<std::size_t>(
                   (l.geo_pt[1] + def::kMathPi<FPType>)*1000000.0);
    }
};

template <typename FPType> struct hash<const SBLoc<FPType> *> {
    std::size_t operator()(const SBLoc<FPType> *l) const {
        return std::hash<SBLoc<FPType>>()(*l);
    }
};
} // namespace std

template <typename FPType>
inline std::ostream &operator<<(std::ostream &os, const SBLoc<FPType> &l) {
    return os << "Lat: " << SBLoc<FPType>::toDegree(l.geo_pt[0])
              << ", Lng: " << SBLoc<FPType>::toDegree(l.geo_pt[1])
              << "\nCity: " << l.city << std::endl
              << "Addr: " << l.addr << std::endl;
}

template <typename FPType>
inline std::istream &operator>>(std::istream &is, SBLoc<FPType> &l) {
    std::getline(is, l.city, ',');
    is >> l.geo_pt[0];
    l.geo_pt[0] = SBLoc<FPType>::toRadians(l.geo_pt[0]);
    is.ignore(256, ',');
    is >> l.geo_pt[1];
    l.geo_pt[1] = SBLoc<FPType>::toRadians(l.geo_pt[1]);
    is.ignore(256, ',');
    std::getline(is, l.addr, '\r');
    return is;
}

#endif /* SBLoc_hpp */
