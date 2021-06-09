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
#include <iterator>
#include <memory>
#include <typeinfo>
#include <random>
#include <regex>
#include <thread>
#include <cstdlib>
#include <cstdint>
#include <map>
#include <cassert>
#include <execution>



template <typename dist_type>
std::vector<Point<dist_type, 2>> GenerateTestLatLngPts(std::size_t num_tests, std::mt19937_64& mt) {
    std::uniform_real_distribution<dist_type> dist(0.0, 1.0);
    std::vector<Point<dist_type, 2>> test_lat_lng_pts;
    test_lat_lng_pts.reserve(num_tests + 4);
    test_lat_lng_pts.insert(test_lat_lng_pts.end(),
                    {{def::kMathPi<dist_type>*0.5, def::kMathPi<dist_type>},
                     {-def::kMathPi<dist_type>*0.5, -def::kMathPi<dist_type>},
                     {-def::kMathPi<dist_type>*0.5, def::kMathPi<dist_type>},
                     {def::kMathPi<dist_type>*0.5, -def::kMathPi<dist_type>}});
    std::generate_n(std::back_inserter(test_lat_lng_pts), num_tests, [&]()->Point<dist_type, 2> {
        return {-0.5*def::kMathPi<dist_type> + def::kMathPi<dist_type>*dist(mt),
                -def::kMathPi<dist_type> + 2.0*def::kMathPi<dist_type>*dist(mt)};});
    return test_lat_lng_pts;
}

template <typename dist_type>
void accuracyTestFromRefFile() {
    
}


template <typename dist_type>
void AccuracyTestFromRefSolver(const std::vector<Point<dist_type, 2>> &test_lat_lng_pts,
                               std::vector<const SBLoc<dist_type>*> &ref_locs,
                               const SBSolver<dist_type>& test_solver,
                               const std::chrono::duration<dist_type> &duration) {
    dist_type test_dist_err_total = 0.0, ref_dist_for_err_pts_total = 0.0;
    std::size_t errot_count = 0;
    auto start_time = std::chrono::steady_clock::now();
    if (ref_locs.empty()) {
        ref_locs.reserve(1ull<<24);
        for (int locIdx = 0; locIdx < test_lat_lng_pts.size() && (std::chrono::steady_clock::now() - start_time < duration); ++locIdx) {
            ref_locs.emplace_back(test_solver.FindNearestLoc(test_lat_lng_pts[locIdx]));
        }
        ref_locs.shrink_to_fit();
        std::cout << "A total number of " << ref_locs.size()
                  << " reference nearest locations were generated by ref solver" << std::endl;
    } else {
        for (int locIdx = 0; locIdx < ref_locs.size() && std::chrono::steady_clock::now() - start_time < duration; ++locIdx) {
            std::chrono::time_point<std::chrono::steady_clock> this_search_start_time = std::chrono::steady_clock::now();
            const SBLoc<dist_type>* testLoc = test_solver.FindNearestLoc(test_lat_lng_pts[locIdx]);
            std::chrono::duration<dist_type, std::micro> this_search_duration = std::chrono::steady_clock::now() - this_search_start_time;
            dist_type test_dist, ref_dist;
            auto test_lat_lng_pt = test_lat_lng_pts[locIdx];
            auto refLoc = ref_locs[locIdx];
            if (testLoc != refLoc
                && (std::fabs((test_dist = testLoc->havDist(test_lat_lng_pt)) - (ref_dist = refLoc->havDist(test_lat_lng_pt))))
                > static_cast<dist_type>(0.000001)) {
                test_dist_err_total += testLoc->havDist(test_lat_lng_pt);
                ref_dist_for_err_pts_total += refLoc->havDist(test_lat_lng_pt);
                ++errot_count;
                auto[test_lat, test_lng] = test_lat_lng_pt.dataArray();
                std::cout << "Test solver this one search time is: " << this_search_duration.count() << "μs" << std::endl
                          << (test_dist < ref_dist ? "Test" : "Ref") << " dist is closer" << std::endl
                          << "Test Point: Lat: " << SBLoc<dist_type>::toDegree(test_lat)
                          << ", Lng: " << SBLoc<dist_type>::toDegree(test_lng) << std::endl
                          << "Test result point: " << *testLoc << "Ref point: " << *refLoc << std::endl;
            }
        }
        std::cout << "A total test of " << ref_locs.size() << " locations" << std::endl
                  << "Num of errors: " << errot_count << std::endl
                  << "Error percentage is: "
                  << errot_count*100.0/ref_locs.size() << "%" << std::endl
                  << "Error in total hav dist diff is: " << std::fabs(test_dist_err_total - ref_dist_for_err_pts_total) << std::endl;
    }
}

template <typename dist_type>
void TimeBuild(const std::shared_ptr<std::vector<SBLoc<dist_type>>> &loc_data_vec,
               SBSolver<dist_type> &solver) {
    std::chrono::time_point<std::chrono::steady_clock> start = std::chrono::steady_clock::now();
    solver.Build(loc_data_vec);
    std::chrono::duration<dist_type> elapsed_time_in_secs = std::chrono::steady_clock::now() - start;
    std::cout << std::regex_replace(typeid(solver).name(),
                                    std::regex("[A-Z]?[0-9]+|.$"), "")
              << std::endl << "Build time: " << elapsed_time_in_secs.count() << std::endl;
    solver.PrintSolverInfo();
}

template <typename dist_type>
void TimeNNSearch(const SBSolver<dist_type> &solver, std::vector<Point<dist_type, 2>> &test_lat_lng_pts,
                  uint_fast64_t seed, std::chrono::duration<dist_type> search_duration_in_secs) {
    std::mt19937_64 mt(seed);
    std::shuffle(test_lat_lng_pts.begin(), test_lat_lng_pts.end(), mt);
    std::vector<std::chrono::duration<dist_type, std::micro>> per_search_time_vec_in_micro_secs;
    std::chrono::duration<dist_type, std::micro> total_elapsed_time_in_micro_secs, max_one_search_time_in_micro_secs{}, min_one_search_time_in_micro_secs(std::numeric_limits<dist_type>::max());
    per_search_time_vec_in_micro_secs.reserve(test_lat_lng_pts.size());
    std::chrono::time_point<std::chrono::steady_clock> startTime = std::chrono::steady_clock::now();
    for (std::size_t loc_idx = 0;
         loc_idx < test_lat_lng_pts.size() && (total_elapsed_time_in_micro_secs = std::chrono::steady_clock::now() - startTime) < search_duration_in_secs;
         ++loc_idx) {
        std::chrono::time_point<std::chrono::steady_clock> this_search_start_time = std::chrono::steady_clock::now();
        solver.FindNearestLoc(test_lat_lng_pts[loc_idx]);
        per_search_time_vec_in_micro_secs.emplace_back(std::chrono::steady_clock::now() - this_search_start_time);
    }
    
    std::for_each(per_search_time_vec_in_micro_secs.cbegin(), per_search_time_vec_in_micro_secs.cend(), [&](auto this_search_time_in_micro_secs) mutable {
        max_one_search_time_in_micro_secs = std::max(max_one_search_time_in_micro_secs, this_search_time_in_micro_secs);
        min_one_search_time_in_micro_secs = std::min(min_one_search_time_in_micro_secs, this_search_time_in_micro_secs);
    });
    
    
    std::cout << std::endl << "Ave per Search Time: "
              << total_elapsed_time_in_micro_secs.count()/per_search_time_vec_in_micro_secs.size()
              << "us from " << per_search_time_vec_in_micro_secs.size() << " searches" << std::endl
              << "Max sinlge search time: " << max_one_search_time_in_micro_secs.count() << std::endl
              << "Min single search time: " << min_one_search_time_in_micro_secs.count() << std::endl << std::endl;
    
}

template <typename dist_type>
void WriteResults(const char* argv[],
                  const std::vector<Point<dist_type, 2>> &test_lat_lng_pts,
                  const SBSolver<dist_type> *solver) {
    std::ofstream outrefLocPtrPerSearchTimePairVecs(argv[2], std::ios::app);
    std::transform(test_lat_lng_pts.cbegin(), test_lat_lng_pts.cend(),
                   std::ostream_iterator<std::string>(outrefLocPtrPerSearchTimePairVecs),
                   [&](const auto &p){
                       const auto resultLoc = solver->FindNearestLoc(p);
                       return std::to_string(p[0]) + " " + std::to_string(p[1]) + " " +
                       std::to_string(resultLoc->geoPt[0]) + " " + std::to_string(resultLoc->geoPt[1]) +
                       resultLoc->city + "," + resultLoc->addr + std::endl;});
    outrefLocPtrPerSearchTimePairVecs.close();
}


template <typename dist_type>
void MainContent(int argc, const char * argv[]) {
    
    std::ifstream infilestream_locs(argv[1]);
    dist_type ave_actual_locs_per_cell = argc < 4 ? def::kDefalutAveActualLocsPerCell : std::stod(argv[3]);
    std::size_t max_cached_cell_vec_size = argc < 5 ? def::kDefaultMaxCacheCellVecSize : std::stoi(argv[4]);
    bool to_test_accuracy = argc < 6 ? def::kDefaultToTestAccuracy : std::tolower(argv[5][0]) == 'y';
    std::chrono::duration<dist_type> accuracy_test_time_in_secs(def::kDefaultAccuracyTestDurationInSecs);
    if (to_test_accuracy) {
        accuracy_test_time_in_secs = std::chrono::duration<dist_type>(std::stoi(std::regex_replace(argv[5], std::regex("[^0-9]*([0-9]+).*"), std::string("$1"))));
    }
    std::chrono::duration<dist_type> search_benchmark_duration_in_secs(def::kDefaultSearchBenchDurationInSecs);
    if (argc >= 6) {
        accuracy_test_time_in_secs = std::chrono::duration<dist_type>(std::stoi(std::regex_replace(argv[6],  std::regex("[^0-9]*([0-9]+).*"), std::string("$1"))));
    }
    std::ifstream inRefLatLngPtLocPairVec;
    if (argc >= 7)
        inRefLatLngPtLocPairVec.open(argv[7]);
    std::ofstream outRefLatLngPtLocPairVec;
    if (argc >= 8)
        outRefLatLngPtLocPairVec.open(argv[8]);
    
        
    std::random_device rd;
    std::seed_seq seq{rd(), rd(), rd(), rd(), rd(), rd(), rd(), rd()};
    std::mt19937_64 mt(seq);
    //mt.seed(686868);

    
    to_test_accuracy = true;
    //maxCacheCellVecSize = (1 << 16ull);
    //maxCacheCellVecSize = (1 << 9ull);
    //aveActualLocsPerCell = 0.2;
    
    auto locData = std::make_shared<std::vector<SBLoc<dist_type>>>();
    locData->reserve(1 << 24);
    locData->assign(std::istream_iterator<SBLoc<dist_type>>(infilestream_locs),
                    std::istream_iterator<SBLoc<dist_type>>());
    infilestream_locs.close();
    assert(locData->size() != 0);
    //std::sort(locData->begin(), locData->end());
    //locData->erase(std::unique(locData->begin(), locData->end()), locData->end());
    locData->shrink_to_fit();
    std::shuffle(locData->begin(), locData->end(), mt);
    
    
    /*
    if (numOfLocsToWriteToFile) {
        auto solver = std::make_unique<BFSBSolver<dist_type>>();
        //auto solver = std::make_unique<BKDTSBSolver<dist_type><KDTree>>();
        timeBuild<dist_type>(locData, *solver);
        writeResults<dist_type>(argv, generatetestSearchLatLngPts<dist_type>(numOfLocsToWriteToFile, mt), solver.get());
        return 0;
    } */
    
    

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
        std::make_unique<UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, dist_type, def::ThreadingPolicy::kSingle>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
        std::make_unique<UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, dist_type, def::ThreadingPolicy::kSingle>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
        std::make_unique<UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, dist_type, def::ThreadingPolicy::kMultiOmp>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
        std::make_unique<UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, dist_type, def::ThreadingPolicy::kMultiOmp>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
    };
    
    uint_fast64_t seed = rd();
    std::vector<Point<dist_type, 2>> search_bench_test_lat_lng_pts = GenerateTestLatLngPts<dist_type>(def::kMaxTestLocs, mt);
    std::vector<const SBLoc<dist_type>*> ref_locs;
    std::vector<Point<dist_type, 2>> accuracy_test_lat_lng_pts;
    if (to_test_accuracy)
        accuracy_test_lat_lng_pts = GenerateTestLatLngPts<dist_type>(def::kMaxTestLocs, mt);

    for (std::size_t i = 0; i < std::size(solvers); ++i) {
        TimeBuild(locData, *solvers[i]);
        TimeNNSearch(*solvers[i], search_bench_test_lat_lng_pts, seed, search_benchmark_duration_in_secs);
        if (to_test_accuracy) {
            AccuracyTestFromRefSolver(accuracy_test_lat_lng_pts, ref_locs, *solvers[i], accuracy_test_time_in_secs);
        }
        solvers[i].reset();
        std::cout << std::endl;
    }
}



int main(int argc, const char * argv[]) {
    if (argc >= 3 && std::strcmp(argv[2], "float") == 0){
        MainContent<float>(argc, argv);
    } else {
        //mainContent<float>(argc, argv);
        MainContent<def::kDefaultDistType>(argc, argv);
    }
    return 0;
}

 
 
 
 
 
