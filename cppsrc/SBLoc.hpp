//
//  SBLoc<dist_type>.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/13/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef SBLoc_hpp
#define SBLoc_hpp

#include "Point.hpp"
#include <stdio.h>
#include <string>
#include <iostream>
#include <math.h>

template <typename dist_type>
struct SBLoc {
    
    dist_type lat;
    dist_type lng;
    std::string city;
    std::string addr;
    
    static constexpr dist_type EARTH_RADIUS = 6371.0;
    
    SBLoc<dist_type>() = default;
    SBLoc<dist_type>(const Point<dist_type, 3>&);
    SBLoc<dist_type>(dist_type lat, dist_type lng);
    
    Point<dist_type, 3> locToCart3DPt() const;
    
    bool operator==(const SBLoc<dist_type> &other) const;
    
    bool operator!=(const SBLoc<dist_type> &other) const;
    
    static dist_type toRadians(dist_type degree);
    
    static dist_type toDegree(dist_type radians);
    
    static dist_type havDist(dist_type lng1, dist_type lat1, dist_type lng2, dist_type lat2);
    
    static dist_type latFromHavDist(dist_type dist, dist_type lat1);
    
    static dist_type lngFromHavDist(dist_type dist, dist_type lng1, dist_type lat);
    
    static Point<dist_type, 3> latLngToCart3DPt(dist_type lat, dist_type lng);
    
    static dist_type xyzDistFromLngLat(dist_type lat1, dist_type lat2, dist_type lngDiff);
};

template <typename dist_type>
inline SBLoc<dist_type>::SBLoc(dist_type lat, dist_type lng) : lat(lat), lng(lng), city(""), addr("") {}

template <typename dist_type>
inline SBLoc<dist_type>::SBLoc(const Point<dist_type, 3> &pt) : lat(asin(pt[2])), lng(asin(pt[1]/cos(asin(pt[2])))), city(""), addr("") {}

template <typename dist_type>
inline bool SBLoc<dist_type>::operator==(const SBLoc<dist_type> &other) const {
    if (lat == other.lat && lng == other.lng)
        return true;
    return std::fabs(lat - other.lat) < static_cast<dist_type>(0.000001) &&
           std::fabs(lng - other.lng) < static_cast<dist_type>(0.000001);
}

template <typename dist_type>
inline bool SBLoc<dist_type>::operator!=(const SBLoc<dist_type> &other) const {
    return !(*this == other);
}

template <typename dist_type>
inline dist_type SBLoc<dist_type>::toDegree(dist_type radians) {
    return radians*180.0/M_PI;
}

template <typename dist_type>
inline dist_type SBLoc<dist_type>::toRadians(dist_type degree) {
    return degree*M_PI/180.0;
}

template <typename dist_type>
inline dist_type SBLoc<dist_type>::havDist(dist_type lng1, dist_type lat1, dist_type lng2, dist_type lat2) {
    dist_type dLat = lat2-lat1, dLon = lng2-lng1;
    dist_type a = sin(dLat/2.0) * sin(dLat/2.0) +
               sin(dLon/2.0) * sin(dLon/2.0) * cos(lat1) * cos(lat2);
    return 2.0 * atan2(sqrt(a), sqrt(1.0-a)) * EARTH_RADIUS;
}

template <typename dist_type>
inline dist_type SBLoc<dist_type>::latFromHavDist(dist_type dist, dist_type lat1) {
    dist_type c = dist/EARTH_RADIUS;
    dist_type sum = sin(c/2.0);
    dist_type dLat = asin(sum)*2.0;
    return dLat + lat1;
}

template <typename dist_type>
inline dist_type SBLoc<dist_type>::lngFromHavDist(dist_type dist, dist_type lng1, dist_type lat) {
    dist_type c = dist/EARTH_RADIUS;
    dist_type sum = sin(c/2);
    dist_type dLng = asin(sum/cos(lat))*2.0;
    return dLng + lng1;
}

template <typename dist_type>
inline Point<dist_type, 3> SBLoc<dist_type>::locToCart3DPt() const {
    return latLngToCart3DPt(lat, lng);
}

template <typename dist_type>
inline Point<dist_type, 3> SBLoc<dist_type>::latLngToCart3DPt(dist_type lat, dist_type lng) {
    return Point<dist_type, 3>{cos(lat)*cos(lng), cos(lat)*sin(lng), sin(lat)};
}

template <typename dist_type>
inline dist_type SBLoc<dist_type>::xyzDistFromLngLat(dist_type lat1, dist_type lat2, dist_type lngDiff) {
    return sqrt(-2 * (cos(lat1)*cos(lat2)*cos(lngDiff) +
                      sin(lat1)*sin(lat2) - 1.0));
}

namespace std {
    template <typename dist_type>
    struct hash<SBLoc<dist_type>> {
        size_t operator()(const SBLoc<dist_type>& l) const {
            return (static_cast<size_t>((l.lat + 0.5*M_PI) * 1000000) << 20) +
                    static_cast<size_t>((l.lng + M_PI) * 1000000);
        }
    };
    
    template <typename dist_type>
    struct hash<const SBLoc<dist_type>*> {
        size_t operator()(const SBLoc<dist_type>* l) const {
            return std::hash<SBLoc<dist_type>>()(*l);
        }
    };
}

template <typename dist_type>
inline std::ostream& operator<<(std::ostream &os, const SBLoc<dist_type> &l) {
    return os << "Lat: " << SBLoc<dist_type>::toDegree(l.lat) << ", Lng: "
              << SBLoc<dist_type>::toDegree(l.lng) << "\nCity: "
              << l.city << std::endl << "Addr: " << l.addr << std::endl;
}

template <typename dist_type>
inline std::istream& operator>>(std::istream &is, SBLoc<dist_type> &l) {
    std::getline(is, l.city, ',');
    is >> l.lat;
    l.lat = SBLoc<dist_type>::toRadians(l.lat);
    is.ignore(256, ',');
    is >> l.lng;
    l.lng = SBLoc<dist_type>::toRadians(l.lng);
    is.ignore(256, ',');
    std::getline(is, l.addr, '\r');
    return is;
}


#endif /* SBLoc_hpp */
