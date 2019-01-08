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
#include "UpgradeBKDTGridSBSolver.hpp"
#include "UniCellBKDTGridSBSolver.hpp"


std::vector<SBLoc> generateTestLocs(size_t numTrials) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<SBLoc> testLocs;
    std::generate_n(std::back_inserter(testLocs), numTrials,
                    [&]()->SBLoc{return SBLoc(SBLoc::toRadians(-90.0 + 180.0*dist(mt)),
                                              SBLoc::toRadians(-180.0 + 360.0*dist(mt)));});
    return testLocs;
}

bool accuracyTest(const std::vector<SBLoc> &testLocs,
                  const std::vector<SBLoc> &testResults,
                  const std::vector<SBLoc> &refResults) {
    double testTotal = 0.0, refTotal = 0.0;
    size_t errorCount = 0, maxTests = std::min(testResults.size(), refResults.size());
    for (size_t i=0; i < maxTests; i++){
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
    std::cout << "A total test of " << maxTests << " locations\n"
              << "Error percentage in num diff is: "
              << (double)errorCount*100/maxTests << "%\n"
              << "Error percentage in hav dist is: " << 100.0*(error-1.0) << "%\n";
    double epsilon = 0.2;
    return error <= 1 + epsilon && error >= 1 - epsilon;
}

void timeBuild(const std::shared_ptr<std::vector<SBLoc>> &locData,
               SBSolver *solver) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    solver->build(locData);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    std::cout << typeid(*solver).name() << " Build time: "
        << elapsedSeconds.count() << std::endl;
}

void timeNN(SBSolver *solver, const std::vector<SBLoc> &testLocs,
           std::vector<SBLoc> &resultLocs, size_t maxResults) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsedSeconds;
    size_t numTrials = 3;
    resultLocs.reserve(maxResults);
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
    } while (elapsedSeconds.count()*1000.0 < 4000 &&
             numTrials * 5 < testLocs.size());
    std::cout << "Time: " << (elapsedSeconds.count()*1000.0)/numTrials
              << " ms per search, " << numTrials << " trials" << std::endl;
}

void writeResults(const char* argv[], const std::vector<SBLoc> &testLocs,
                  SBSolver *solver) {
    std::ofstream outRefResults(argv[2], std::ios::app);
    std::transform(testLocs.begin(), testLocs.end(),
                   std::ostream_iterator<std::string>(outRefResults),
                   [&](const SBLoc &l){
                       auto resultLoc = solver->findNearest(l.lng, l.lat);
                       return to_string(l.lat) + " " + to_string(l.lng) + " " +
                              to_string(resultLoc->lat) + " " + to_string(resultLoc->lng) +
                              resultLoc->city + "," + resultLoc->addr + "\n";});
}


int main(int argc, const char * argv[]) {
    std::ifstream infileLocs(argv[1]), inRefResults(argv[2]);
    double aveLocPerCell = argc == 3 ? 1 : std::stod(argv[3]);
    size_t numOfLocsToWriteToFile = argc < 5 ? false : std::stoi(argv[4]);
    
    infileLocs.ignore(256, '\r');
    infileLocs.ignore(256, '\r');
    
    auto locData = std::make_shared<std::vector<SBLoc>>();
    std::copy(std::istream_iterator<SBLoc>(infileLocs),
              std::istream_iterator<SBLoc>(), std::back_inserter(*locData));
    std::stable_sort(locData->begin(), locData->end(), [](const auto &l1,
        const auto& l2){return l1.lng*10000000+l1.lat<l2.lng*10000000+l2.lat;});
    locData->erase(locData->begin(),
                  std::unique(locData->rbegin(), locData->rend()).base());
    std::random_shuffle(locData->begin(), locData->end());
    
    
    /*
    
    auto bfSolver = std::make_shared<BFSBSolver>();
    timeBuild(locData, bfSolver.get());
    auto refResult = bfSolver->findNearest(SBLoc::toRadians(-100.603), SBLoc::toRadians(56.6682));
    auto bkdtSolver = std::make_shared<BKDTSBSolver>();
    timeBuild(locData, bkdtSolver.get());
    for (int i = 0; i < 1; ++i) {
    auto result = bkdtSolver->findNearest(SBLoc::toRadians(-100.603), SBLoc::toRadians(56.6682));
    std::cout << "Result: " << *result << std::endl
              << "Ref result: " << *refResult << std::endl;
    }
    return 0;
    
    */
    
    
    
    
    if (numOfLocsToWriteToFile) {
        //auto bfSolver = std::make_shared<BFSBSolver>();
        auto bfSolver = std::make_shared<BKDTSBSolver>();
        timeBuild(locData, bfSolver.get());
        writeResults(argv, generateTestLocs(numOfLocsToWriteToFile), bfSolver.get());
        return 0;
    }
    
    
    size_t MAX_TRIALS = 0xFFFFFF;
    std::vector<SBLoc> testResults, refResults, testLocs;
    
    /*
    std::string line;
    while (std::getline(inRefResults, line)) {
        std::stringstream ss(line);
        SBLoc test, ref;
        std::string cityAddr;
        ss >> test.lat >> test.lng >> ref.lat >> ref.lng;
        std::getline(inRefResults, ref.city, ',');
        std::getline(inRefResults, ref.addr);
        testLocs.push_back(test);
        refResults.push_back(ref);
    } */
    auto restTestLocs = generateTestLocs(MAX_TRIALS-testLocs.size());
    testLocs.insert(testLocs.begin(),{SBLoc(M_PI/2, M_PI),SBLoc(-M_PI/2, -M_PI),
        SBLoc(-M_PI/2, M_PI), SBLoc(M_PI/2, -M_PI)});
    testLocs.insert(testLocs.end(), std::make_move_iterator(restTestLocs.begin()),
                    std::make_move_iterator(restTestLocs.end()));
    
    std::vector<shared_ptr<SBSolver>> solvers{
        //std::make_shared<BFSBSolver>(),
        //std::make_shared<KDTSBSolver>(),
        std::make_shared<BKDTSBSolver>(),
        //std::make_shared<GridSBSolver>(),
        //std::make_shared<BKDTGridSBSolver>(aveLocPerCell),
        std::make_shared<UniCellBKDTGridSBSolver>(aveLocPerCell),
        //std::make_shared<UpgradeBKDTGridSBSolver>(0.7*aveLocPerCell),
    };
    
    
    
    refResults.clear();
    
    
    std::for_each(solvers.begin(), solvers.end(),
                  [&](const auto &solver){timeBuild(locData, solver.get());});
    for (size_t i = 0; i < solvers.size(); ++i) {
        //timeBuild(locData, solvers[i].get());
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

 
 
 
 
 
