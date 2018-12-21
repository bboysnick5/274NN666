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
#include "GridSBSolver.hpp"
#include "BKDTGridSBSolver.hpp"


std::vector<SBLoc> generateTestLocs(size_t numTrials) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<SBLoc> testLocs;
    std::generate_n(std::back_inserter(testLocs), numTrials,
                    [&]()->SBLoc{ return SBLoc(24.0 + 25.0*dist(mt), -125.0 + 73.0*dist(mt));});
    return testLocs;
}

bool accuracyTest(const std::vector<SBLoc> &testLocs,
                  const std::vector<SBLoc> &testResults,
                  const std::vector<SBLoc> &refResults) {
    double testTotal = 0.0, refTotal = 0.0;
    size_t errorCount = 0;
    for (size_t i=0; i < refResults.size(); i++){
        const auto &testLoc = testLocs[i], &refResult = refResults[i],
                   &testResult = testResults[i];
        double refDist = SBLoc::havDist(refResult.lng, refResult.lat, testLoc.lng, testLoc.lat);
        refTotal += refDist;
        if (testResult != refResult) {
            errorCount++;
            std::cout << "Test Point: " << testLoc << std::endl
                      << "Return point: " << testResult << "Ref point: " << refResult
                      << std::endl;
            testTotal += SBLoc::havDist(testResult.lng, testResult.lat, testLoc.lng, testLoc.lat);
        } else {
            testTotal += refDist;
        }
    }
    
    double error = testTotal/refTotal;
    std::cout << "A total test of " << testResults.size() << " locations\n"
              << "Error percentage in num diff is: "
              << (double)errorCount*100/refResults.size() << "%\n"
              << "Error percentage in hav dist is: " << 100.0*(error-1.0) << "%\n";
    double epsilon = 0.2;
    return error <= 1 + epsilon && error >= 1 - epsilon;
}

void timeBuild(const std::shared_ptr<std::vector<SBLoc>> &sbData,
               SBSolver *solver) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    solver->build(sbData);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    std::cout << typeid(*solver).name() << " Build time: "
        << elapsedSeconds.count() << std::endl;
}

void timeNN(SBSolver *solver, const std::vector<SBLoc> &testLocs,
           std::vector<SBLoc> &resultLocs, size_t maxResults) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsedSeconds;
    size_t numTrials = 2;
    do {
        numTrials *= 5;
        start = std::chrono::system_clock::now();
        if (resultLocs.size() < maxResults) {
            std::transform(testLocs.begin() + resultLocs.size(),
                           testLocs.begin() + resultLocs.size() + numTrials,
                           std::back_inserter(resultLocs), [&](const SBLoc &l){
                               return *solver->findNearest(l.lng, l.lat);});
        } else {
            std::for_each(testLocs.begin() + maxResults,
                          testLocs.begin() + maxResults + numTrials,
                          [&](const SBLoc &l){
                              solver->findNearest(l.lng, l.lat);});
        }
        end = std::chrono::system_clock::now();
        elapsedSeconds = end - start;
    } while (elapsedSeconds.count()*1000.0 < 8000 && numTrials * 5 < testLocs.size());
    std::cout << "Time: " << (elapsedSeconds.count()*1000.0)/numTrials
              << " ms per search, " << numTrials << " trials" << std::endl;
}

int main(int argc, const char * argv[]) {
    
    std::ifstream sbFile(argv[1]);
    sbFile.ignore(256, '\r');
    sbFile.ignore(256, '\r');
    
    auto sbData = std::make_shared<std::vector<SBLoc>>();
    std::copy(std::istream_iterator<SBLoc>(sbFile),
              std::istream_iterator<SBLoc>(), std::back_inserter(*sbData));
    std::stable_sort(sbData->begin(), sbData->end(), [](const auto &l1,
        const auto& l2){return l1.lng*10000000+l1.lat<l2.lng*10000000+l2.lat;});
    sbData->erase(sbData->begin(),
                  std::unique(sbData->rbegin(), sbData->rend()).base());
    std::random_shuffle(sbData->begin(), sbData->end());
    
    double aveLocPerCell = argc == 2 ? 1 : std::stod(argv[2]);
    
    size_t MAX_TRIALS = 0xFFFFF;
    std::vector<SBLoc> testResults, refResults, testLocs = generateTestLocs(MAX_TRIALS);
    
    std::vector<shared_ptr<SBSolver>> solvers{
        std::make_shared<BFSBSolver>(),
        std::make_shared<KDTSBSolver>(),
        std::make_shared<BKDTSBSolver>(),
        //std::make_shared<GridSBSolver>(),
        std::make_shared<BKDTGridSBSolver>(aveLocPerCell),
    };
    
    for (size_t i = 0; i < solvers.size(); ++i) {
        timeBuild(sbData, solvers[i].get());
        timeNN(solvers[i].get(), testLocs, testResults,
               refResults.size() == 0 ? MAX_TRIALS : refResults.size());
        if (i == 0) {
            refResults = std::move(testResults);
        } else {
            accuracyTest(testLocs, testResults, refResults);
        }
        std::cout << std::endl;
        testResults.clear();
    }
    
    return 0;
}

 
 
 
 
 
