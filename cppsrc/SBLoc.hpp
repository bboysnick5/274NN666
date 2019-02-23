//
//  SBLoc.hpp
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

struct SBLoc {
    
    double lat;
    double lng;
    std::string city;
    std::string addr;
    
    static constexpr double EARTH_RADIUS = 6371;
    
    SBLoc() = default;
    SBLoc(const Point<3>&);
    SBLoc(double lat, double lng);
    
    Point<3> locToCart3DPt() const;
    
    bool operator==(const SBLoc &other) const;
    
    bool operator!=(const SBLoc &other) const;
    
    static double toRadians(double degree);
    
    static double toDegree(double radians);
    
    static double havDist(double lng1, double lat1, double lng2, double lat2);
    
    static double latFromHavDist(double dist, double lat1);
    
    static double lngFromHavDist(double dist, double lng1, double lat);
    
    static Point<3> latLngToCart3DPt(double lat, double lng);
    
    static double xyzDistFromLngLat(double lat1, double lat2, double lngDiff);
};

inline SBLoc::SBLoc(double lat, double lng) : lat(lat), lng(lng), city(""), addr("") {}

inline SBLoc::SBLoc(const Point<3> &pt) : lat(asin(pt[2])), lng(asin(pt[1]/cos(asin(pt[2])))), city(""), addr("") {}


inline bool SBLoc::operator==(const SBLoc &other) const {
    if (lat == other.lat && lng == other.lng)
        return true;
    return std::fabs(lat - other.lat) < 0.000001 &&
           std::fabs(lng - other.lng) < 0.000001;
}

inline bool SBLoc::operator!=(const SBLoc &other) const {
    return !(*this == other);
}

inline double SBLoc::toDegree(double radians) {
    return radians*180.0/M_PI;
}

inline double SBLoc::toRadians(double degree) {
    return degree*M_PI/180.0;
}

inline double SBLoc::havDist(double lng1, double lat1, double lng2, double lat2) {
    double dLat = lat2-lat1, dLon = lng2-lng1;
    double a = sin(dLat/2.0) * sin(dLat/2.0) +
               sin(dLon/2.0) * sin(dLon/2.0) * cos(lat1) * cos(lat2);
    return 2.0 * atan2(sqrt(a), sqrt(1-a)) * EARTH_RADIUS;
}

inline double SBLoc::latFromHavDist(double dist, double lat1) {
    double c = dist/EARTH_RADIUS;
    double sum = sin(c/2);
    double dLat = asin(sum)*2;
    return dLat + lat1;
}

inline double SBLoc::lngFromHavDist(double dist, double lng1, double lat) {
    double c = dist/EARTH_RADIUS;
    double sum = sin(c/2);
    double dLng = asin(sum/cos(lat))*2;
    return dLng + lng1;
}

inline Point<3> SBLoc::locToCart3DPt() const {
    return latLngToCart3DPt(lat, lng);
}

inline Point<3> SBLoc::latLngToCart3DPt(double lat, double lng) {
    return Point<3>{cos(lat)*cos(lng), cos(lat)*sin(lng), sin(lat)};
}

inline double SBLoc::xyzDistFromLngLat(double lat1, double lat2, double lngDiff) {
    return sqrt(-2 * (cos(lat1)*cos(lat2)*cos(lngDiff) +
                      sin(lat1)*sin(lat2) - 1.0));
}

namespace std {
    template <>
    struct hash<SBLoc> {
        size_t operator()(const SBLoc& l) const {
            return (static_cast<size_t>((l.lat + 0.5*M_PI) * 1000000) << 20) +
                    static_cast<size_t>((l.lng + M_PI) * 1000000);
        }
    };
    
    template <>
    struct hash<const SBLoc*> {
        size_t operator()(const SBLoc* l) const {
            return std::hash<SBLoc>()(*l);
        }
    };
}

inline std::ostream& operator<<(std::ostream &os, const SBLoc &l) {
    return os << "Lat: " << SBLoc::toDegree(l.lat) << ", Lng: "
              << SBLoc::toDegree(l.lng) << "\nCity: "
              << l.city << std::endl << "Addr: " << l.addr << std::endl;
}

inline std::istream& operator>>(std::istream &is, SBLoc &l) {
    std::getline(is, l.city, ',');
    is >> l.lat;
    l.lat = SBLoc::toRadians(l.lat);
    is.ignore(256, ',');
    is >> l.lng;
    l.lng = SBLoc::toRadians(l.lng);
    is.ignore(256, ',');
    std::getline(is, l.addr, '\r');
    return is;
}


#endif /* SBLoc_hpp */
