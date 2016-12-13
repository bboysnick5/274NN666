//
//  main.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/9/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <chrono>
#include <memory>
#include <typeinfo>
#include <random>
#include "KDTSBSolver.hpp"
#include "BFSBSolver.hpp"
#include "BKDTSBSolver.hpp"


bool accuracyTest(SBSolver *testSolver, SBSolver *refSolver) {
    double testTotal = 0.0, refTotal = 0.0;
    int numTrials = 1000;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (int i=0; i < numTrials; i++){
        double x = -125.0 + 73.0*dist(mt);
        double y = 24.0 + 25.0*dist(mt);
        auto testLoc = testSolver->findNearest(x, y);
        auto refLoc = refSolver->findNearest(x, y);
        testTotal += SBLoc::distance(testLoc.lng, testLoc.lat, x, y);
        refTotal += SBLoc::distance(refLoc.lng, refLoc.lat, x, y);
    }
    
    double error = testTotal/refTotal;
    std::cout << "Error percentage is: " << 100.0*(error-1.0) << std::endl;;
    double epsilon = 0.2;
    return error <= 1 + epsilon && error >= 1 - epsilon;
}

void timeBuild(const std::vector<SBLoc> &sbData, SBSolver *solver) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    solver->build(sbData);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    std::cout << typeid(*solver).name() << " Build time: "
        << elapsedSeconds.count() << std::endl;
}

void timeNN(SBSolver *solver) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsedSeconds;
    int numTrials = 1000;
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    do {
        start = std::chrono::system_clock::now();
        for(int i=0; i<numTrials; i++){
            double x = -125.0 + 73.0*dist(mt);
            double y = 24.0 + 25.0*dist(mt);
            auto tmp = solver->findNearest(x, y);
        }
        end = std::chrono::system_clock::now();
        elapsedSeconds = end - start;
        numTrials *= 10;
    } while (elapsedSeconds.count()*1000.0 < 5000 && numTrials < 1000000);
    std::cout << "Time: " << (elapsedSeconds.count()*1000.0)/numTrials
        << " ms per search, " << numTrials << " trials" << std::endl;
}


int main(int argc, const char * argv[]) {
    
    std::ifstream sbFile(argv[1]);
    std::string line;

    sbFile.ignore(256, '\r');
    sbFile.ignore(256, '\r');
    
    std::vector<SBLoc> sbData;
    std::copy(std::istream_iterator<SBLoc>(sbFile),
              std::istream_iterator<SBLoc>(), std::back_inserter(sbData));
    std::random_shuffle(sbData.begin(), sbData.end());
    
    std::vector<shared_ptr<SBSolver>> solvers{
        std::make_shared<BFSBSolver>(),
        std::make_shared<KDTSBSolver>(),
        std::make_shared<BKDTSBSolver>()
    };
    
    for (int i = 0; i < solvers.size(); ++i) {
        //sbData = std::vector<SBLoc>(sbData.begin(), sbData.begin()+sbData.size()/2);
        timeBuild(sbData, solvers[i].get());
        if (i != 0 &&!accuracyTest(solvers[i].get(), solvers[0].get()))
            continue;
        timeNN(solvers[i].get());
        std::cout << std::endl;
    }
    
    
    
    
    return 0;
}

 
 
 
 
 
