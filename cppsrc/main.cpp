//
//  main.cpp
//  274F16NearestSB
//
//  Created by Yunlong Liu on 12/9/16.
//  Copyright © 2016 Yunlong Liu. All rights reserved.
//

#include "BFEUCPtSolver.hpp"
#include "BFSBSolver.hpp"
#include "KDTSBSolver.hpp"
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
#include <regex>

std::vector<std::pair<double, double>> generateTestLocs(size_t numTrials, std::mt19937& mt) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<std::pair<double, double>> testLocs;
    testLocs.reserve(numTrials+4);
    testLocs.insert(testLocs.cbegin(),{{M_PI/2, M_PI},{-M_PI/2, -M_PI},
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
        if (*testResult != *refResult) {
            errorCount++;
            std::cout << "Test Point lng: " << SBLoc::toDegree(testLng)
                      << ", lat: " << SBLoc::toDegree(testLat)
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
               const std::shared_ptr<SBSolver> &solver) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    solver->build(locData);
    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsedSeconds = end - start;
    auto &solverRef = *solver.get();
    std::cout << std::regex_replace(typeid(solverRef).name(),
                                    std::regex("[A-Z]?[0-9]+|.$"), "")
              << "\nBuild time: " << elapsedSeconds.count() << "\n";
    solver->printSolverInfo();
    std::cout << "\n";
}

void timeNN(const std::shared_ptr<SBSolver> &solver,
            const std::vector<std::pair<double, double>> &testLocs,
            std::vector<const SBLoc*> &resultLocs, size_t maxResults) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsedSeconds;
    size_t numTrials = 3;
    resultLocs.reserve(maxResults);
    do {
        numTrials *= 5;
        start = std::chrono::system_clock::now();
        if (resultLocs.size() < maxResults) {
            std::transform(testLocs.cbegin() + resultLocs.size(),
                           testLocs.cbegin() + resultLocs.size() + numTrials,
                           std::back_inserter(resultLocs), [&](const auto &p){
                               return solver->findNearest(p.second, p.first);});
        } else {
            std::for_each(testLocs.cbegin() + maxResults,
                          testLocs.cbegin() + maxResults + numTrials,
                          [&](const auto &p){
                              solver->findNearest(p.second, p.first);});
        }
        end = std::chrono::system_clock::now();
        elapsedSeconds = end - start;
    } while (elapsedSeconds.count()*1000.0 < 4000 &&
             numTrials * 5 < testLocs.size());
    resultLocs.shrink_to_fit();
    auto &solverRef = *solver.get();
    std::cout << std::regex_replace(typeid(solverRef).name(),
                                    std::regex("[A-Z]?[0-9]+|.$"), "")
              << "\nSearch Time: " << (elapsedSeconds.count()*1000.0)/numTrials
              << " ms per search, " << numTrials << " trials\n";
}

void writeResults(const char* argv[],
                  const std::vector<std::pair<double, double>> &testLocs,
                  const std::shared_ptr<SBSolver>& solver) {
    std::ofstream outRefResults(argv[2], std::ios::app);
    std::transform(testLocs.cbegin(), testLocs.cend(),
                   std::ostream_iterator<std::string>(outRefResults),
                   [&](const auto &p){
                       const auto resultLoc = solver->findNearest(p.second, p.first);
                       return std::to_string(p.first) + " " + std::to_string(p.second) + " " +
                       std::to_string(resultLoc->lat) + " " + std::to_string(resultLoc->lng) +
                              resultLoc->city + "," + resultLoc->addr + "\n";});
    outRefResults.close();
}


int main(int argc, const char * argv[]) {

    //for (int i = 0; i < argc; ++i)
   //     printf(argv[i]);
    //return 0;
    //std::cout << std::is_pod<Point<3>>::value;
    
    //std::cout << std::is_const<typename std::remove_pointer<typename std::iterator_traits<std::vector<int>::const_iterator>::pointer>::type>::value << "\n";
    
   //std::cout << "Size of: " << sizeof(std::variant<KDT<KDTreeCusMem>, std::vector<std::pair<Point<3>, const SBLoc*>>, const SBLoc*>) << std::endl;
   // return 0;
    
    std::random_device rd;
    std::mt19937 mt(rd());
    
    size_t MAX_TRIALS = 0xFFFFFF;
    std::vector<const SBLoc*> testResults, refResults;
    std::vector<std::pair<double, double>> testLocs = generateTestLocs(MAX_TRIALS, mt);
    
    std::ifstream infileLocs(argv[1]), inRefResults(argv[2]);
    double aveLocPerCell = argc < 4 ? 0.8 : std::stod(argv[3]);
    size_t maxCacheCellVecSize = argc < 5 ? 1200 : std::stoi(argv[4]);
    size_t numOfLocsToWriteToFile = argc < 6 ? false : std::stoi(argv[5]);
    
    infileLocs.ignore(256, '\r');
    infileLocs.ignore(256, '\r');
    
    
    auto locData = std::make_shared<std::vector<SBLoc>>();
    std::copy(std::istream_iterator<SBLoc>(infileLocs),
              std::istream_iterator<SBLoc>(), std::back_inserter(*locData));
    infileLocs.close();
    locData->shrink_to_fit();
    std::shuffle(locData->begin(), locData->end(), mt);
    
    /*
    std::stable_sort(locData->begin(), locData->end(), [](const auto &l1,
        const auto& l2){return l1.lng*10000000+l1.lat<l2.lng*10000000+l2.lat;});
    locData->erase(locData->begin(),
                   std::unique(locData->rbegin(), locData->rend()).base());
    locData->shrink_to_fit();
    infileLocs.close();
    
    std::random_device rd;
    std::mt19937 g(rd());
    std::shuffle(locData->begin(), locData->end(), g);
    
    std::for_each(locData->begin(), locData->end(), []( auto& l){
        l.city = std::regex_replace(l.city, std::regex("\n|\r"), "");
        l.addr = std::regex_replace(l.addr, std::regex("\n|\r"), "");
    });
    
    std::ofstream outFile(argv[2]);
    std::transform(locData->begin(), locData->end(),
                   std::ostream_iterator<std::string>(outFile),
                   [&](const auto &l){
                       return l.city + "," +to_string(SBLoc::toDegree( l.lat)) + "," + to_string(SBLoc::toDegree( l.lng)) + "," + l.addr + "\r";});
    
    return 0;
    
     */
    if (numOfLocsToWriteToFile) {
        auto solver = std::make_shared<BFSBSolver>();
        //auto solver = std::make_shared<BKDTSBSolver<KDTree>>();
        timeBuild(locData, solver);
        writeResults(argv, generateTestLocs(numOfLocsToWriteToFile, mt), solver);
        return 0;
    }
    


    
    std::vector<std::shared_ptr<SBSolver>> solvers{
        //std::make_shared<BFSBSolver>(),
        //std::make_shared<BFEUCPtSBSolver>(),
        //std::make_shared<KDTSBSolver<KDTree>>(),
        //std::make_shared<BKDTSBSolver<KDTree>>(),
        //std::make_shared<BKDTSBSolver<KDTreeCusMem>>(),
        std::make_shared<BKDTSBSolver<KDTreeExpandLongest>>(),
        //std::make_shared<GridSBSolver>(),
        //std::make_shared<BKDTGridSBSolver>(aveLocPerCell),
        //std::make_shared<UniLatLngBKDTGridSBSolver<KDTree>>(0.85*aveLocPerCell, maxCacheCellVecSize),
        //std::make_shared<UniLatLngBKDTGridSBSolver<KDTreeCusMem>>(0.85*aveLocPerCell, maxCacheCellVecSize),
        //std::make_shared<UniCellBKDTGridSBSolver<KDTree>>(aveLocPerCell, maxCacheCellVecSize),
        // std::make_shared<UniCellBKDTGridSBSolver<KDTreeCusMem>>(aveLocPerCell, maxCacheCellVecSize),
        std::make_shared<UniCellBKDTGridSBSolver<KDTreeExpandLongest>>(aveLocPerCell, maxCacheCellVecSize),
    };
    
    
    std::for_each(solvers.cbegin(), solvers.cend(),
                  [&](const auto &solver) {timeBuild(locData, solver);});
    for (size_t i = 0; i < solvers.size(); ++i) {
        //timeBuild(locData, solvers[i].get());
        timeNN(solvers[i], testLocs, testResults,
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

 
 
 
 
 
