//
//  SBSolver.hpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/10/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#ifndef SBSolver_hpp
#define SBSolver_hpp

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream>
#include <math.h>
#include "Point.hpp"


struct SBLoc {
    double lng;
    double lat;
    std::string city;
    std::string addr;
    
    inline bool operator==(const SBLoc &other) {
        return other.lng == lng && other.lat == lat;
    }
    
    inline static double distance(double lng1, double lat1, double lng2, double lat2){
        double R = 6371; // radius of earth, in km
        double dLat = (lat2-lat1)*M_PI/180;
        double dLon = (lng2-lng1)*M_PI/180;
        lat1 = lat1*M_PI/180;
        lat2 = lat2*M_PI/180;
        
        double a = sin(dLat/2) * sin(dLat/2) +
        sin(dLon/2) * sin(dLon/2) * cos(lat1) * cos(lat2);
        double c = 2 * atan2(sqrt(a), sqrt(1-a));
        double d = R * c;
        
        return d;
    }
    
};


inline std::istream& operator>>(std::istream &is, SBLoc &l) {
    std::getline(is, l.city, ',');
    is >> l.lat;
    is.ignore(256, ',');
    is >> l.lng;
    is.ignore(256, ',');
    std::getline(is, l.addr, '\r');
    return is;
}


class SBSolver {
public:
    
    inline Point<3> transLatLngToXYZPt(double lng, double lat) {
        Point<3> p;
        p[0] = cos(lat*M_PI/180)*cos(lng*M_PI/180);
        p[1] = cos(lat*M_PI/180)*sin(lng*M_PI/180);
        p[2] = sin(lat*M_PI/180);
        return p;
    }
    
    virtual void build(const std::vector<SBLoc> &sbData) = 0;
    virtual SBLoc findNearest(double lng, double lat) = 0;
    
};




#endif /* SBSolver_hpp */
