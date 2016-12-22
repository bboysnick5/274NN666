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
    
    double lng;
    double lat;
    std::string city;
    std::string addr;
    
    static constexpr double EARTH_RADIUS = 6371;
        
    inline bool operator==(const SBLoc &other) const {
        return other.lng == lng && other.lat == lat;
    }
    
    inline bool operator!=(const SBLoc &other) const {
        return other.lng != lng || other.lat != lat;
    }
    
    inline static double
    havDist(double lng1, double lat1, double lng2, double lat2) {
        double dLat = (lat2-lat1)*M_PI/180;
        double dLon = (lng2-lng1)*M_PI/180;
        lat1 = lat1*M_PI/180;
        lat2 = lat2*M_PI/180;
        double a = sin(dLat/2) * sin(dLat/2)
                   + sin(dLon/2) * sin(dLon/2) * cos(lat1) * cos(lat2);
        double c = 2 * atan2(sqrt(a), sqrt(1-a));
        return EARTH_RADIUS * c;
    }
    
    inline static Point<3> latLngToCart3DXYZ(double lng, double lat) {
        Point<3> p;
        p[0] = cos(lat*M_PI/180)*cos(lng*M_PI/180);
        p[1] = cos(lat*M_PI/180)*sin(lng*M_PI/180);
        p[2] = sin(lat*M_PI/180);
        return p;
    }
};

namespace std {
    template <>
    struct hash<SBLoc> {
        size_t operator()(const SBLoc& l) const {
            return static_cast<long long>(l.lat * 1000000 +
                   std::numeric_limits<double>::epsilon())* 100000000 +
                   static_cast<long long>(l.lng * 1000000 +
                   std::numeric_limits<double>::epsilon());
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