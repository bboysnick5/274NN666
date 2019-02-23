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
#include "UnionUniLatLngBKDTGridSBSolver.hpp"
#include "UnionUniCellBKDTGridSBSolver.hpp"




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
#include <thread>
#include <cstdlib>



std::vector<Point<2>> generateTestLocs(size_t numTrials, std::mt19937_64& mt) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<Point<2>> testLocs;
    testLocs.reserve(numTrials+4);
    testLocs.insert(testLocs.cbegin(),{{M_PI/2, M_PI},{-M_PI/2, -M_PI},
                    {-M_PI/2, M_PI}, {M_PI/2, -M_PI}});
    
    std::generate_n(std::back_inserter(testLocs), numTrials, [&]()->Point<2>{
        return {SBLoc::toRadians(-90.0 + 180.0*dist(mt)),
                SBLoc::toRadians(-180.0 + 360.0*dist(mt))};});
    return testLocs;
}

void accuracyTest(const std::vector<Point<2>> &testLocs,
                  const std::vector<const SBLoc*> &testResults,
                  const std::vector<const SBLoc*> &refResults) {
    double testTotal = 0.0, refTotal = 0.0;
    size_t errorCount = 0, maxTests = std::min(testResults.size(), refResults.size());
    for (size_t i = 0; i < maxTests; ++i) {
        const auto &refResult = refResults[i],
                   &testResult = testResults[i];
        const auto &[testLat, testLng] = testLocs[i].dataArray();
        double refDist = SBLoc::havDist(refResult->lng, refResult->lat, testLng, testLat);
        refTotal += refDist;
        if (*testResult == *refResult) {
            testTotal += refDist;
        } else {
            testTotal += SBLoc::havDist(testResult->lng, testResult->lat, testLng, testLat);
            ++errorCount;
            std::cout << "Test Point lng: " << SBLoc::toDegree(testLng)
            << ", lat: " << SBLoc::toDegree(testLat)
            << "\nReturn point: " << *testResult
            << "Ref point: " << *refResult << "\n";
        }
    }
    
    double error = testTotal/refTotal;
    std::cout << "A total test of " << maxTests << " locations\n"
              << "Error percentage in num diff is: "
              << errorCount*100.0/maxTests << "%\n"
              << "Error percentage in hav dist is: " << 100.0*(error-1.0) << "%\n";
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
            const std::vector<Point<2>> &testLocs,
            std::vector<const SBLoc*> &resultLocs, size_t maxResults,
            unsigned int seed, bool testAccuracy = true) {
    if (!testAccuracy)
        maxResults = 0;
    std::mt19937_64 mt(seed);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<double> elapsedSeconds;
    size_t numTrials = 64;
    resultLocs.reserve(maxResults);
    do {
        numTrials *= 4;
        if (testAccuracy && resultLocs.size() < maxResults) {
            start = std::chrono::system_clock::now();
            std::transform(testLocs.cbegin() + resultLocs.size(),
                           testLocs.cbegin() + resultLocs.size() + numTrials,
                           std::back_inserter(resultLocs), [&](const auto &p){
                               return solver->findNearest(p[0], p[1]);});
            end = std::chrono::system_clock::now();
        } else {
            //std::shuffle(testLocs.begin(), testLocs.end(), mt);
            start = std::chrono::system_clock::now();
            std::for_each_n(testLocs.cbegin() + maxResults, numTrials,
                            [&](const auto &p){solver->findNearest(p[0], p[1]);});
            end = std::chrono::system_clock::now();
        }
        elapsedSeconds = end - start;
    } while (elapsedSeconds.count()*1000.0 < 4000.0 && numTrials * 5 < testLocs.size());
    resultLocs.shrink_to_fit();
    auto &polySolverRef = *solver.get();
    std::cout << std::regex_replace(typeid(polySolverRef).name(),
                                    std::regex("[A-Z]?[0-9]+|.$"), "")
              << "\nSearch Time: " << (elapsedSeconds.count()*1000.0)/numTrials
              << " ms per search, " << numTrials << " trials\n";
}

void writeResults(const char* argv[],
                  const std::vector<Point<2>> &testLocs,
                  const std::shared_ptr<SBSolver>& solver) {
    std::ofstream outRefResults(argv[2], std::ios::app);
    std::transform(testLocs.cbegin(), testLocs.cend(),
                   std::ostream_iterator<std::string>(outRefResults),
                   [&](const auto &p){
                       const auto resultLoc = solver->findNearest(p[0], p[1]);
                       return std::to_string(p[0]) + " " + std::to_string(p[1]) + " " +
                       std::to_string(resultLoc->lat) + " " + std::to_string(resultLoc->lng) +
                       resultLoc->city + "," + resultLoc->addr + "\n";});
    outRefResults.close();
}



int main(int argc, const char * argv[]) {
    
    
    /*
    struct Cell {
        
        Cell(size_t, const std::vector<std::pair<Point<3>, const SBLoc*>>&);
        ~Cell();
        size_t size() const;
        
    private:
        size_t _size;
        union {
            const SBLoc *cacheLoc;
            std::pair<Point<3>, const SBLoc*>* cacheLocs;
            //KDT<KDTType>* cacheTree;
        };
    };
    
    std::cout << std::is_standard_layout<Cell>::value;
    return 0;
     */
    
    //int *x = new int[20];
    //std::cout << sizeof(std::pair<double, std::pair<const Point<3>*, const SBLoc*>>);
    //delete[] x;
    //return 0; */
    
    /*
    std::ifstream infileLocs(argv[1]);
    std::ofstream outfileLocs(argv[2]);
    std::string line;
    while (std::getline(infileLocs, line)) {
        std::istringstream iss(line, '\t');
        std::string id, token, desc;
        std::getline(iss, id, '\t');
        std::getline(iss, token, '\t');
        desc = token;
        std::getline(iss, token, '\t');
        desc += " " + token;
        std::getline(iss, token, '\t');
        desc += " " + token;
        double lat, lng;
        std::getline(iss, token, '\t');
        lat = std::stod(token);
        std::getline(iss, token, '\t');
        lng = std::stod(token);
        outfileLocs << id + "," + std::to_string(lat) + ","
                       + std::to_string(lng) + "," + desc + "\r";
    }
    return 0;
    
    */
    
    std::random_device rd;
    std::seed_seq seq{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937_64 mt(seq);
    //mt.seed(686868);
    
    size_t MAX_TRIALS = 0xFFFFFF;
    std::vector<const SBLoc*> testResults, refResults;
    std::vector<Point<2>> testLocs = generateTestLocs(MAX_TRIALS, mt);
    
    std::ifstream infileLocs(argv[1]), inRefResults(argv[2]);
    double aveLocPerCell = argc < 4 ? 0.4 : std::stod(argv[3]);
    size_t MAX_CACHE_CELL_VEC_SIZE = argc < 5 ? 1200 : std::stoi(argv[4]);
    bool testAccuracy = argc < 6 ? true : std::tolower(argv[5][0]) == 'y';
    size_t numOfLocsToWriteToFile = argc < 7 ? false : std::stoi(argv[6]);
    
    //infileLocs.ignore(256, '\r');
    //infileLocs.ignore(256, '\r');
    
    
    auto locData = std::make_shared<std::vector<SBLoc>>();
    locData->reserve(0xFFFFFF);
    locData->assign(std::istream_iterator<SBLoc>(infileLocs),
                    std::istream_iterator<SBLoc>());
    locData->shrink_to_fit();
    infileLocs.close();
    std::shuffle(locData->begin(), locData->end(), mt);
    
    
    /*
    std::stable_sort(locData->begin(), locData->end(), [](const auto &l1,
        const auto& l2){return l1.lng*100000000+l1.lat<l2.lng*100000000+l2.lat;});
    locData->erase(locData->begin(),
                   std::unique(locData->rbegin(), locData->rend()).base());
    locData->shrink_to_fit();
    infileLocs.close();
    
    //std::random_device rd;
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
                       return l.city + "," +std::to_string(SBLoc::toDegree(l.lat)) + "," + std::to_string(SBLoc::toDegree(l.lng)) + "," + l.addr + "\r";});
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
        std::make_shared<BKDTSBSolver<KDTreeExpandLongestVec>>(),
        //std::make_shared<GridSBSolver>(),
        //std::make_shared<BKDTGridSBSolver>(aveLocPerCell),
        //std::make_shared<UniLatLngBKDTGridSBSolver<KDTree>>(0.85*aveLocPerCell, maxCacheCellVecSize),
        //std::make_shared<UniLatLngBKDTGridSBSolver<KDTreeCusMem>>(0.85*aveLocPerCell, maxCacheCellVecSize),
        //std::make_shared<UniLatLngBKDTGridSBSolver<KDTreeExpandLongest>>(0.85*aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        //std::make_shared<UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        //std::make_shared<UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        //std::make_shared<UniCellBKDTGridSBSolver<KDTree>>(aveLocPerCell, maxCacheCellVecSize),
        // std::make_shared<UniCellBKDTGridSBSolver<KDTreeCusMem>>(aveLocPerCell, maxCacheCellVecSize),
        std::make_shared<UniCellBKDTGridSBSolver<KDTreeExpandLongest>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        std::make_shared<UniCellBKDTGridSBSolver<KDTreeExpandLongestVec>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        std::make_shared<UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        std::make_shared<UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
    };
    
    
   // std::for_each(solvers.cbegin(), solvers.cend(),
         //         [&](const auto &solver) {timeBuild(locData, solver);});
    //testAccuracy = false;
    unsigned int timeNNSeed = rd();
    for (size_t i = 0; i < solvers.size(); ++i) {
        using namespace std::chrono_literals;
        std::this_thread::sleep_for(100ms);
        timeBuild(locData, solvers[i]);
        std::this_thread::sleep_for(100ms);
        timeNN(solvers[i], testLocs, testResults,
               refResults.size() == 0 ? MAX_TRIALS : refResults.size(), timeNNSeed, testAccuracy);
        solvers[i].reset();
        
        if (testAccuracy) {
            if (i == 0) {
                refResults = std::move(testResults);
            } else {
                accuracyTest(testLocs, testResults, refResults);
            }
            testResults.clear();
        }
        std::cout << "\n";
    }
    
    return 0;
}

 
 
 
 
 
