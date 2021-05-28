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


template <typename dist_type>
std::vector<Point<dist_type, 2>> generateTestLocs(size_t numTests, std::mt19937_64& mt) {
    std::uniform_real_distribution<dist_type> dist(0.0, 1.0);
    std::vector<Point<dist_type, 2>> testLocs;
    testLocs.reserve(numTests + 4);
    testLocs.insert(testLocs.end(),
                    {{Def::PI<dist_type>*0.5, Def::PI<dist_type>},
                     {-Def::PI<dist_type>*0.5, -Def::PI<dist_type>},
                     {-Def::PI<dist_type>*0.5, Def::PI<dist_type>},
                     {Def::PI<dist_type>*0.5, -Def::PI<dist_type>}});
    std::generate_n(std::back_inserter(testLocs), numTests, [&]()->Point<dist_type, 2> {
        return {SBLoc<dist_type>::toRadians(-90.0 + 180.0*dist(mt)),
                SBLoc<dist_type>::toRadians(-180.0 + 360.0*dist(mt))};});
    return testLocs;
}

template <typename dist_type>
void accuracyTest(const std::vector<Point<dist_type, 2>> &testLocs,
                  const std::vector<const SBLoc<dist_type>*> &testResults,
                  const std::vector<const SBLoc<dist_type>*> &refResults) {
    dist_type testErrTotal = 0.0, refTotal = 0.0;
    size_t errorCount = 0, maxTests = std::min(testResults.size(), refResults.size());
    for (size_t i = 0; i < maxTests; ++i) {
        const auto &refResult = refResults[i],
                   &testResult = testResults[i];
        const auto &[testLat, testLng] = testLocs[i].dataArray();
        dist_type testDist, refDist;
        if (*testResult != *refResult
            && (std::fabs((testDist = testResult->havDistComp({testLat, testLng}))
                          - (refDist = refResult->havDistComp({testLat, testLng})))
                > static_cast<dist_type>(0.000001))) {
            testErrTotal += testResult->havDist({testLat, testLng});
            refTotal += refResult->havDist({testLat, testLng});
            ++errorCount;
            std::cout << "Test Point lng: " << SBLoc<dist_type>::toDegree(testLng)
            << ", lat: " << SBLoc<dist_type>::toDegree(testLat)
            << "\nReturn point: " << *testResult
            << "Ref point: " << *refResult << "\n";
        }
    }
    
    std::cout << "A total test of " << maxTests << " locations\n"
              << "Num of errors: " << errorCount << std::endl
              << "Error percentage in num diff is: "
              << errorCount*100.0/maxTests << "%\n"
              << "Error in hav dist is: " << std::fabs(testErrTotal - refTotal) << "\n";
}

template <typename dist_type>
void timeBuild(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &locData,
               SBSolver<dist_type> &solver) {
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();
    solver.build(locData);
    end = std::chrono::system_clock::now();
    std::chrono::duration<dist_type> elapsedSeconds = end - start;
    //auto &solverRef = *solver.get();
    std::cout << std::regex_replace(typeid(solver).name(),
                                    std::regex("[A-Z]?[0-9]+|.$"), "")
              << "\nBuild time: " << elapsedSeconds.count() << "\n";
    solver.printSolverInfo();
    std::cout << "\n";
}

template <typename dist_type>
void timeNN(const SBSolver<dist_type> &solver,
            const std::vector<Point<dist_type, 2>> &testLocs,
            std::vector<const SBLoc<dist_type>*> &resultLocs, size_t maxResults,
            unsigned int seed, bool testAccuracy = true) {
    if (!testAccuracy)
        maxResults = 0;
    std::mt19937_64 mt(seed);
    std::chrono::time_point<std::chrono::system_clock> start, end;
    std::chrono::duration<dist_type> elapsedSeconds{};
    size_t numTrials = 4, totalTrials = 0;
    resultLocs.reserve(maxResults);
    do {
        numTrials *= 4;
        totalTrials += numTrials;
        if (testAccuracy && resultLocs.size() < maxResults) {
            start = std::chrono::system_clock::now();
            std::transform(testLocs.cbegin() + resultLocs.size(),
                           testLocs.cbegin() + resultLocs.size() + numTrials,
                           std::back_inserter(resultLocs), [&](const auto &p){
                               return solver.findNearest(p);});
            end = std::chrono::system_clock::now();
        } else {
            //std::shuffle(testLocs.begin(), testLocs.end(), mt);
            start = std::chrono::system_clock::now();
            std::for_each_n(testLocs.cbegin() + maxResults, numTrials,
                            [&](const auto &p){solver.findNearest(p);});
            end = std::chrono::system_clock::now();
        }
        elapsedSeconds += end - start;
    } while (elapsedSeconds.count()*1000.0 < 4000.0 && numTrials * 5 < testLocs.size());
    resultLocs.shrink_to_fit();
    //auto &polySolverRef = *solver.get();
    std::cout << std::regex_replace(typeid(solver).name(),
                                    std::regex("[A-Z]?[0-9]+|.$"), "")
              << "\nSearch Time: " << (elapsedSeconds.count()*1000.0)/totalTrials
              << " ms per search, " << totalTrials << " trials\n";
}

template <typename dist_type>
void writeResults(const char* argv[],
                  const std::vector<Point<dist_type, 2>> &testLocs,
                  const SBSolver<dist_type> *solver) {
    std::ofstream outRefResults(argv[2], std::ios::app);
    std::transform(testLocs.cbegin(), testLocs.cend(),
                   std::ostream_iterator<std::string>(outRefResults),
                   [&](const auto &p){
                       const auto resultLoc = solver->findNearest(p);
                       return std::to_string(p[0]) + " " + std::to_string(p[1]) + " " +
                       std::to_string(resultLoc->geoPt[0]) + " " + std::to_string(resultLoc->geoPt[1]) +
                       resultLoc->city + "," + resultLoc->addr + "\n";});
    outRefResults.close();
}


int main(int argc, const char * argv[]) {
    
    //std::bitset<64> bit(~1ull);
    //std::cout << bit.to_string() << std::endl;
    
    using dist_type = double;
    
    std::random_device rd;
    std::seed_seq seq{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937_64 mt(seq);
    //mt.seed(686868);
    
    size_t MAX_TRIALS = (1 << 24);
    std::vector<const SBLoc<dist_type>*> testResults, refResults;
    std::vector<Point<dist_type, 2>> testLocs = generateTestLocs<dist_type>(MAX_TRIALS, mt);
    
    //std::cout << argv[1] << std::endl;
    
    std::ifstream infileLocs(argv[1]), inRefResults(argv[2]);
    dist_type aveLocPerCell = argc < 4 ? 0.4 : std::stod(argv[3]);
    size_t MAX_CACHE_CELL_VEC_SIZE = argc < 5 ? 1200 : std::stoi(argv[4]);
    bool testAccuracy = argc < 6 ? false : std::tolower(argv[5][0]) == 'y';
    size_t numOfLocsToWriteToFile = argc < 7 ? false : std::stoi(argv[6]);
    
    testAccuracy = true;
    //MAX_CACHE_CELL_VEC_SIZE = (1 << 16);
    //aveLocPerCell = 10;
    
    //infileLocs.ignore(256, '\r');
    //infileLocs.ignore(256, '\r');
    
    
    auto locData = std::make_shared<std::vector<SBLoc<dist_type>>>();
    locData->reserve(1 << 24);
    locData->assign(std::istream_iterator<SBLoc<dist_type>>(infileLocs),
                    std::istream_iterator<SBLoc<dist_type>>());
    locData->shrink_to_fit();
    std::shuffle(locData->begin(), locData->end(), mt);
    infileLocs.close();
    
    if (numOfLocsToWriteToFile) {
        auto solver = std::make_unique<BFSBSolver<dist_type>>();
        //auto solver = std::make_unique<BKDTSBSolver<dist_type><KDTree>>();
        timeBuild<dist_type>(locData, *solver);
        writeResults<dist_type>(argv, generateTestLocs<dist_type>(numOfLocsToWriteToFile, mt), solver.get());
        return 0;
    }
    

    std::unique_ptr<SBSolver<dist_type>> solvers[] = {
        //std::make_unique<BFSBSolver<dist_type>>(),
        //std::make_unique<BFEUCPtSBSolver<dist_type>>(),
        //std::make_unique<KDTSBSolver<KDTree,dist_type>>(),
        //std::make_unique<BKDTSBSolver<KDTree, dist_type>>(),
        //std::make_unique<BKDTSBSolver<dist_type><KDTreeCusMem>>(),
        //std::make_unique<BKDTSBSolver<KDTreeExpandLongest, dist_type>>(),
        std::make_unique<BKDTSBSolver<KDTreeExpandLongestVec, dist_type>>(),
        //std::make_unique<GridSBSolver<dist_type>>(),
        //std::make_unique<BKDTGridSBSolver<dist_type>>(aveLocPerCell),
        //std::make_unique<UniLatLngBKDTGridSBSolver<KDTree, dist_type>>(0.85*aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        //std::make_unique<UniLatLngBKDTGridSBSolver<KDTreeCusMem, dist_type>>(0.85*aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        //std::make_unique<UniLatLngBKDTGridSBSolver<KDTreeExpandLongest, dist_type>>(0.85*aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        //std::make_unique<UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec,dist_type>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        //std::make_unique<UnionUniLatLngBKDTGridSBSolver<dist_type><KDTreeExpandLongestVec>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        //std::make_unique<UniCellBKDTGridSBSolver<dist_type><KDTree>>(aveLocPerCell, maxCacheCellVecSize),
        // std::make_unique<UniCellBKDTGridSBSolver<dist_type><KDTreeCusMem>>(aveLocPerCell, maxCacheCellVecSize),
        //std::make_unique<UniCellBKDTGridSBSolver<KDTreeExpandLongest, dist_type>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        //std::make_unique<UniCellBKDTGridSBSolver<KDTreeExpandLongestVec, dist_type>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        std::make_unique<UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, dist_type, Def::Threading_Policy::SINGLE>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        std::make_unique<UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, dist_type, Def::Threading_Policy::SINGLE>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        std::make_unique<UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, dist_type, Def::Threading_Policy::MULTI_OMP>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),
        std::make_unique<UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, dist_type, Def::Threading_Policy::MULTI_OMP>>(aveLocPerCell, MAX_CACHE_CELL_VEC_SIZE),

    };
    
    //std::vector<std::unique_ptr<SBSolver<dist_type>>> solvers(std::make_move_iterator(std::begin(move_only_solvers)), std::make_move_iterator(std::end(move_only_solvers)));
    
    
    
   // std::for_each(solvers.cbegin(), solvers.cend(),
         //         [&](const auto &solver) {timeBuild(locData, solver);});
   // testAccuracy = false;
    unsigned int timeNNSeed = rd();
    for (size_t i = 0; i < std::size(solvers); ++i) {
        using namespace std::chrono_literals;
        std::this_thread::sleep_for(300ms);
        timeBuild(locData, *solvers[i]);
        std::this_thread::sleep_for(300ms);
        timeNN(*solvers[i], testLocs, testResults,
               refResults.size() == 0 ? MAX_TRIALS : refResults.size(), timeNNSeed, testAccuracy);
        solvers[i].reset();
        
        if (testAccuracy) {
            if (i == 0) {
                refResults = std::move(testResults);
            } else {
                accuracyTest<dist_type>(testLocs, testResults, refResults);
            }
            testResults.clear();
        }
        std::cout << "\n";
    }
    
    return 0;
}

 
 
 
 
 
