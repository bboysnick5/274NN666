//
//  main.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/9/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//


#include "KDTSBSolver.hpp"
#include "BFSBSolver.hpp"
#include "BKDTSBSolver.hpp"
#include "GridSBSolver.hpp"
#include "BKDTGridSBSolver.hpp"
#include "UniLatLngBKDTGridSBSolver.hpp"
#include "UniCellBKDTGridSBSolver.hpp"

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

std::vector<std::pair<double, double>> generateTestLocs(size_t numTrials) {
    std::random_device rd;
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<std::pair<double, double>> testLocs;
    testLocs.reserve(numTrials+4);
    testLocs.insert(testLocs.begin(),{{M_PI/2, M_PI},{-M_PI/2, -M_PI},
        {-M_PI/2, M_PI}, {M_PI/2, -M_PI}});
    
    std::generate_n(std::back_inserter(testLocs), numTrials,
                    [&]()->std::pair<double, double>{
                        return {SBLoc::toRadians(-90.0 + 180.0*dist(mt)),
                                SBLoc::toRadians(-180.0 + 360.0*dist(mt))};});
    return testLocs;
}

bool accuracyTest(const std::vector<std::pair<double, double>> &testLocs,
                  const std::vector<const SBLoc*> &testResults,
                  const std::vector<const SBLoc*> &refResults) {
    double testTotal = 0.0, refTotal = 0.0;
    size_t errorCount = 0, maxTests = std::min(testResults.size(), refResults.size());
    for (size_t i = 0; i < maxTests; ++i) {
        const auto &refResult = refResults[i],
                   &testResult = testResults[i];
        const auto[testLat, testLng] = testLocs[i];
        double refDist = SBLoc::havDist(refResult->lng, refResult->lat, testLng, testLat);
        refTotal += refDist;
        if (testResult != refResult) {
            errorCount++;
            std::cout << "Test Point lat: " << SBLoc::toDegree(testLat)
                      << ", lng: " << SBLoc::toDegree(testLng)
                      << "\nReturn point: " << *testResult
                      << "Ref point: " << *refResult << "\n";
            testTotal += SBLoc::havDist(testResult->lng, testResult->lat, testLng, testLat);
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
              << elapsedSeconds.count() << "\n";
    solver->printSolverInfo();
    std::cout << "\n";
}

void timeNN(SBSolver *solver, const std::vector<std::pair<double, double>> &testLocs,
            std::vector<const SBLoc*> &resultLocs, size_t maxResults) {
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
                           std::back_inserter(resultLocs), [&](const auto &p){
                               return solver->findNearest(p.second, p.first);});
        } else {
            std::for_each(testLocs.begin() + maxResults,
                          testLocs.begin() + maxResults + numTrials,
                          [&](const auto &p){
                              solver->findNearest(p.second, p.first);});
        }
        end = std::chrono::system_clock::now();
        elapsedSeconds = end - start;
    } while (elapsedSeconds.count()*1000.0 < 4000 &&
             numTrials * 5 < testLocs.size());
    resultLocs.shrink_to_fit();
    std::cout << "Time: " << (elapsedSeconds.count()*1000.0)/numTrials
              << " ms per search, " << numTrials << " trials\n";
}

void writeResults(const char* argv[],
                  const std::vector<std::pair<double, double>> &testLocs,
                  SBSolver *solver) {
    std::ofstream outRefResults(argv[2], std::ios::app);
    std::transform(testLocs.begin(), testLocs.end(),
                   std::ostream_iterator<std::string>(outRefResults),
                   [&](const auto &p){
                       const auto resultLoc = solver->findNearest(p.second, p.first);
                       return to_string(p.first) + " " + to_string(p.second) + " " +
                       to_string(resultLoc->lat) + " " + to_string(resultLoc->lng) +
                              resultLoc->city + "," + resultLoc->addr + "\n";});
}


int main(int argc, const char * argv[]) {
    
    
    //std::cout << "Size of: " << sizeof(KDTree<3, int*>) << std::endl;
    //return 0;
    
    size_t MAX_TRIALS = 0xFFFFFF;
    std::vector<const SBLoc*> testResults, refResults;
    std::vector<std::pair<double, double>> testLocs = generateTestLocs(MAX_TRIALS) ;
    
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
    
    
    std::ifstream infileLocs(argv[1]), inRefResults(argv[2]);
    double aveLocPerCell = argc == 3 ? 0.8 : std::stod(argv[3]);
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
    locData->shrink_to_fit();
    infileLocs.close();
    
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(locData->begin(), locData->end(), g);
    
    
    
    
    if (numOfLocsToWriteToFile) {
        //auto solver = std::make_shared<BFSBSolver>();
        auto solver = std::make_shared<BKDTSBSolver<KDTree>>();
        timeBuild(locData, solver.get());
        writeResults(argv, generateTestLocs(numOfLocsToWriteToFile), solver.get());
        return 0;
    }
    


    
    std::vector<shared_ptr<SBSolver>> solvers{
        //std::make_shared<BFSBSolver>(),
        //std::make_shared<KDTSBSolver<KDTree>>(),
        //std::make_shared<BKDTSBSolver<KDTree>>(),
        //std::make_shared<BKDTSBSolver<KDTreeCusMem>>(),
        //std::make_shared<GridSBSolver>(),
        //std::make_shared<BKDTGridSBSolver>(aveLocPerCell),
        //std::make_shared<UniLatLngBKDTGridSBSolver<KDTree>>(0.85*aveLocPerCell),
        //std::make_shared<UniLatLngBKDTGridSBSolver<KDTreeCusMem>>(0.85*aveLocPerCell),
        //std::make_shared<UniCellBKDTGridSBSolver<KDTree>>(aveLocPerCell),
        std::make_shared<UniCellBKDTGridSBSolver<KDTreeCusMem>>(aveLocPerCell),
    };
    
    
    std::for_each(solvers.begin(), solvers.end(),
                  [&](const auto &solver) {timeBuild(locData, solver.get());});
    for (size_t i = 0; i < solvers.size(); ++i) {
        //timeBuild(locData, solvers[i].get());
        timeNN(solvers[i].get(), testLocs, testResults,
               refResults.size() == 0 ? MAX_TRIALS : refResults.size());
        if (i == 0) {
            refResults = std::move(testResults);
        } else {
            accuracyTest(testLocs, testResults, refResults);
        }
        testResults.clear();
        std::cout << "\n";
    }
    
    return 0;
}

 
 
 
 
 
