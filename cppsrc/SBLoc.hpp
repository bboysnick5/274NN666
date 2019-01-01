//
//  SBLoc.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/13/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef SBLoc_hpp
#define SBLoc_hpp

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
    SBLoc(double lat, double lng);
    
    Point<3, DistType::EUC> locToCart3DPt() const;
    
    bool operator==(const SBLoc &other) const;
    
    bool operator!=(const SBLoc &other) const;
    
    static double toRadians(double degree);
    
    static double havDist(double lng1, double lat1, double lng2, double lat2);
    
    static double latFromHavDist(double dist, double lat1);
    
    static double lngFromHavDist(double dist, double lng1, double lat);
    
    static Point<3, DistType::EUC> latLngToCart3DPt(double lng, double lat);
    
    static double xyzDistFromLngLat(double lat1, double lat2, double lngDiff);
};

inline SBLoc::SBLoc(double lat, double lng) : lat(lat), lng(lng) {}

inline bool SBLoc::operator==(const SBLoc &other) const {
    if (other.lng == lng && other.lat == lat)
        return true;
    return std::fabs(other.lng - lng) < 0.00001 &&
           std::fabs(other.lat - lat) < 0.00001;
}

inline bool SBLoc::operator!=(const SBLoc &other) const {
    return !(*this == other);
}

inline double SBLoc::toRadians(double degree) {
    return degree*M_PI/180.0;
}

inline double SBLoc::havDist(double lng1, double lat1, double lng2, double lat2) {
    double dLat = (lat2-lat1)*M_PI/180;
    double dLon = (lng2-lng1)*M_PI/180;
    lat1 = lat1*M_PI/180;
    lat2 = lat2*M_PI/180;
    double a = sin(dLat/2) * sin(dLat/2)
    + sin(dLon/2) * sin(dLon/2) * cos(lat1) * cos(lat2);
    double c = 2 * atan2(sqrt(a), sqrt(1-a));
    return EARTH_RADIUS * c;
}

inline double SBLoc::latFromHavDist(double dist, double lat1) {
    double c = dist/EARTH_RADIUS;
    double sum = sin(c/2);
    double dLat = asin(sum)*2;
    return dLat*180/M_PI + lat1;
}

inline double SBLoc::lngFromHavDist(double dist, double lng1, double lat) {
    double c = dist/EARTH_RADIUS;
    double sum = sin(c/2);
    double dLng = asin(sum/cos(lat/180*M_PI))*2;
    return dLng*180/M_PI + lng1;
}

inline Point<3, DistType::EUC> SBLoc::locToCart3DPt() const {
    return latLngToCart3DPt(lng, lat);
}

inline Point<3, DistType::EUC> SBLoc::latLngToCart3DPt(double lng, double lat) {
    lng = toRadians(lng);
    lat = toRadians(lat);
    return Point<3, DistType::EUC>{cos(lat)*cos(lng), cos(lat)*sin(lng),
                                   sin(lat)};
}

inline double SBLoc::xyzDistFromLngLat(double lat1, double lat2, double lngDiff) {
    lat1 *= M_PI/180.0;
    lat2 *= M_PI/180.0;
    return sqrt(-2 * (cos(lat1)*cos(lat2)*cos(lngDiff*M_PI/180.0) +
                      sin(lat1)*sin(lat2) - 1.0));
}

namespace std {
    template <>
    struct hash<SBLoc> {
        size_t operator()(const SBLoc& l) const {
            return (static_cast<size_t>((l.lat + 90.0) * 1000000) << 20) +
                   static_cast<size_t>((l.lng + 180.0) * 1000000);
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
    return os << "Lng: " << l.lng << ", Lat: " << l.lat << std::endl << "City: "
              << l.city << std::endl << "Addr: " << l.addr << std::endl;
}

inline std::istream& operator>>(std::istream &is, SBLoc &l) {
    std::getline(is, l.city, ',');
    is >> l.lat;
    is.ignore(256, ',');
    is >> l.lng;
    is.ignore(256, ',');
    std::getline(is, l.addr, '\r');
    return is;
}


#endif /* SBLoc_hpp */
