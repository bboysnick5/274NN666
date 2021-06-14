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
#include <stdio.h>
#include <string>
#include <iostream>

template <typename FPType>
struct SBLoc {
    
    PointND<FPType, 2> geoPt;
    std::string city;
    std::string addr;
    
    static constexpr FPType EARTH_RADIUS = FPType(6371.0);
    
    SBLoc() = default;
    SBLoc(const PointND<FPType, 3>&);
    
    bool operator<(const SBLoc<FPType> &other) const;
    
    bool operator==(const SBLoc<FPType> &other) const;
    
    bool operator!=(const SBLoc<FPType> &other) const;
    
    static FPType toRadians(FPType degree);
    
    static FPType toDegree(FPType radians);
    
    static FPType havDist(const SBLoc<FPType>&, const SBLoc<FPType>&);
        
    static FPType havDist(const PointND<FPType, 2>& p1, const PointND<FPType, 2>& p2);
    
    FPType havDist(const SBLoc<FPType>& other) const;
    
    FPType havDist(const PointND<FPType, 2>& other) const;
    
    FPType havDistComp(const PointND<FPType, 2>& other) const;

    static FPType deltaLatOnSameLngFromHavDist(FPType dist);
    
    static FPType lngFromSameLatHavDist(FPType dist, FPType lng1, FPType lat);
    
    PointND<FPType, 3> locToCart3DPt() const;
    
    static PointND<FPType, 3> geoPtToCart3DPt(const PointND<FPType, 2>&);
    
    static FPType EUC3DDistFromLatDeltaLng(FPType lat1, FPType lat2, FPType deltaLng);
    
    static FPType EUC3DDistSqFromLatDeltaLng(FPType lat1, FPType lat2, FPType deltaLng);

};


template <typename FPType>
inline SBLoc<FPType>::SBLoc(const PointND<FPType, 3> &pt) :
geoPt(std::asin(pt[2]), lng(std::asin(pt[1]/std::cos(std::asin(pt[2]))))), city(""), addr("") {}


template <typename FPType>
inline bool SBLoc<FPType>::operator<(const SBLoc<FPType> &other) const {
    return geoPt[0] < other.geoPt[0];
}

template <typename FPType>
inline bool SBLoc<FPType>::operator==(const SBLoc<FPType> &other) const {
    if (geoPt == other.geoPt)
        return true;
    return std::fabs(geoPt[0] - other.geoPt[0]) < static_cast<FPType>(0.000001) &&
           std::fabs(geoPt[1] - other.geoPt[1]) < static_cast<FPType>(0.000001);
}

template <typename FPType>
inline bool SBLoc<FPType>::operator!=(const SBLoc<FPType> &other) const {
    return !(*this == other);
}

template <typename FPType>
inline FPType SBLoc<FPType>::toDegree(FPType radians) {
    return radians*180.0/def::kMathPi<FPType>;
}

template <typename FPType>
inline FPType SBLoc<FPType>::toRadians(FPType degree) {
    return degree*def::kMathPi<FPType>/180.0;
}

template <typename FPType>
FPType SBLoc<FPType>::havDist(const SBLoc<FPType>& other) const {
    return havDist(geoPt, other.geoPt);
}

template <typename FPType>
FPType SBLoc<FPType>::havDist(const PointND<FPType, 2>& otherGeoPt) const {
    return PointND<FPType, 2>::template dist<PointND<FPType, 2>::DistType::HAV>(geoPt, otherGeoPt, EARTH_RADIUS);
}

template <typename FPType>
FPType SBLoc<FPType>::havDistComp(const PointND<FPType, 2>& otherGeoPt) const {
    return PointND<FPType, 2>::template dist<PointND<FPType, 2>::DistType::HAVCOMP>(geoPt, otherGeoPt, EARTH_RADIUS);
}

template <typename FPType>
inline FPType SBLoc<FPType>::havDist(const SBLoc<FPType>& l1, const SBLoc<FPType>& l2) {
    return PointND<FPType, 2>::template dist<PointND<FPType, 2>::DistType::HAV>(l1.geoPt, l2.geoPt, EARTH_RADIUS);
}

template <typename FPType>
inline FPType SBLoc<FPType>::havDist(const PointND<FPType, 2>& p1, const PointND<FPType, 2>& p2) {
    return PointND<FPType, 2>::template dist<PointND<FPType, 2>::DistType::HAV>(p1, p2, EARTH_RADIUS);
}

template <typename FPType>
inline FPType SBLoc<FPType>::deltaLatOnSameLngFromHavDist(FPType dist) {
    FPType c = dist/EARTH_RADIUS;
    FPType sum = sin(c*0.5);
    FPType dLat = std::asin(sum)*2.0;
    return dLat;
}

template <typename FPType>
inline FPType SBLoc<FPType>::lngFromSameLatHavDist(FPType dist, FPType lng1, FPType lat) {
    FPType c = dist/EARTH_RADIUS;
    FPType sum = sin(c*0.5);
    FPType dLng = std::asin(sum/std::cos(lat))*2.0;
    return dLng + lng1;
}

template <typename FPType>
inline PointND<FPType, 3> SBLoc<FPType>::locToCart3DPt() const {
    return geoPtToCart3DPt(geoPt);
}

template <typename FPType>
inline PointND<FPType, 3> SBLoc<FPType>::geoPtToCart3DPt(const PointND<FPType, 2>& geoPt) {
    return PointND<FPType, 3>{std::cos(geoPt[0])*std::cos(geoPt[1]), std::cos(geoPt[0])*sin(geoPt[1]), sin(geoPt[0])};
}

template <typename FPType>
inline FPType SBLoc<FPType>::EUC3DDistFromLatDeltaLng(FPType lat1, FPType lat2, FPType deltaLng) {
    return sqrt(EUC3DDistSqFromLatDeltaLng(lat1, lat2, deltaLng));
}

template <typename FPType>
inline FPType SBLoc<FPType>::EUC3DDistSqFromLatDeltaLng(FPType lat1, FPType lat2, FPType deltaLng) {
    return 2.0 - 2.0 * (sin(lat1)*sin(lat2) + std::cos(lat1)*std::cos(lat2)*std::cos(deltaLng));
}

namespace std {
    template <typename FPType>
    struct hash<SBLoc<FPType>> {
        std::size_t operator()(const SBLoc<FPType>& l) const {
            return (static_cast<std::size_t>((l.geoPt[0] + 0.5*def::kMathPi<FPType>) * 1000000.0) << 20) +
                    static_cast<std::size_t>((l.geoPt[1] + def::kMathPi<FPType>) * 1000000.0);
        }
    };
    
    template <typename FPType>
    struct hash<const SBLoc<FPType>*> {
        std::size_t operator()(const SBLoc<FPType>* l) const {
            return std::hash<SBLoc<FPType>>()(*l);
        }
    };
}

template <typename FPType>
inline std::ostream& operator<<(std::ostream &os, const SBLoc<FPType> &l) {
    return os << "Lat: " << SBLoc<FPType>::toDegree(l.geoPt[0]) << ", Lng: "
              << SBLoc<FPType>::toDegree(l.geoPt[1]) << "\nCity: "
              << l.city << std::endl << "Addr: " << l.addr << std::endl;
}

template <typename FPType>
inline std::istream& operator>>(std::istream &is, SBLoc<FPType> &l) {
    std::getline(is, l.city, ',');
    is >> l.geoPt[0];
    l.geoPt[0] = SBLoc<FPType>::toRadians(l.geoPt[0]);
    is.ignore(256, ',');
    is >> l.geoPt[1];
    l.geoPt[1] = SBLoc<FPType>::toRadians(l.geoPt[1]);
    is.ignore(256, ',');
    std::getline(is, l.addr, '\r');
    return is;
}


#endif /* SBLoc_hpp */
