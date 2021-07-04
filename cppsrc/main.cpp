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

//#include <benchmark/benchmark.h>
#include <absl/random/random.h>

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
#include <utility>
#include <map>
#include <cassert>
#include <concepts>
#include <span>
//#include <ranges>



template <typename FPType>
std::vector<PointND<FPType, 2>> GenerateTestLatLngPts(std::size_t num_tests, absl::BitGen &bitgen) {
    std::vector<PointND<FPType, 2>> test_lat_lng_pts;
    test_lat_lng_pts.reserve(num_tests + 4);
    test_lat_lng_pts.insert(test_lat_lng_pts.end(),
                    {{def::kMathPi<FPType>*0.5, def::kMathPi<FPType>},
                     {-def::kMathPi<FPType>*0.5, -def::kMathPi<FPType>},
                     {-def::kMathPi<FPType>*0.5, def::kMathPi<FPType>},
                     {def::kMathPi<FPType>*0.5, -def::kMathPi<FPType>}});
    std::generate_n(std::back_inserter(test_lat_lng_pts), num_tests, [&bitgen]()->PointND<FPType, 2> {
        return {-0.5*def::kMathPi<FPType> + def::kMathPi<FPType>*absl::Uniform(absl::IntervalClosed, bitgen, 0.0, 1.0),
                -def::kMathPi<FPType> + def::kMathPi<FPType>*absl::Uniform(absl::IntervalClosed, bitgen, 0.0, 2.0)};});
    return test_lat_lng_pts;
}

template <typename FPType>
void accuracyTestFromRefFile() {
}


template <typename FPType>
void AccuracyTestFromRefSolver(const std::vector<PointND<FPType, 2>> &test_lat_lng_pts,
                               std::vector<const SBLoc<FPType>*> &ref_locs,
                               const SBSolver<FPType>& test_solver,
                               const std::chrono::duration<FPType> &duration) {
    FPType test_dist_err_total = 0.0, ref_dist_for_err_pts_total = 0.0;
    std::size_t errot_count = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    if (ref_locs.empty()) {
        ref_locs.reserve(1ull<<24);
        for (std::size_t loc_idx = 0; loc_idx < test_lat_lng_pts.size() && (std::chrono::high_resolution_clock::now() - start_time < duration); ++loc_idx) {
            ref_locs.emplace_back(test_solver.FindNearestLoc(test_lat_lng_pts[loc_idx]));
        }
        ref_locs.shrink_to_fit();
        std::cout << "A total number of " << ref_locs.size()
                  << " reference nearest locations were generated by ref solver" << std::endl;
    } else {
        std::size_t loc_idx = 0;
        for (; loc_idx < ref_locs.size() && std::chrono::high_resolution_clock::now() - start_time < duration; ++loc_idx) {
            std::chrono::time_point<std::chrono::high_resolution_clock> this_search_start_time = std::chrono::high_resolution_clock::now();
            const SBLoc<FPType>* testLoc = test_solver.FindNearestLoc(test_lat_lng_pts[loc_idx]);
            std::chrono::duration<FPType, std::micro> this_search_duration = std::chrono::high_resolution_clock::now() - this_search_start_time;
            FPType test_dist, ref_dist;
            auto test_lat_lng_pt = test_lat_lng_pts[loc_idx];
            auto refLoc = ref_locs[loc_idx];
            if (testLoc != refLoc
                && (std::fabs((test_dist = testLoc->havDist(test_lat_lng_pt))
                              - (ref_dist = refLoc->havDist(test_lat_lng_pt))))
                    > std::numeric_limits<FPType>::epsilon()) {
                test_dist_err_total += testLoc->havDist(test_lat_lng_pt);
                ref_dist_for_err_pts_total += refLoc->havDist(test_lat_lng_pt);
                ++errot_count;
                auto[test_lat, test_lng] = test_lat_lng_pt.dataArray();
                std::cout << "Test solver this one search time is: " << this_search_duration.count() << "us" << std::endl
                          << (test_dist < ref_dist ? "Test" : "Ref") << " dist is closer" << std::endl
                          << "Test PointND: Lat: " << SBLoc<FPType>::toDegree(test_lat)
                          << ", Lng: " << SBLoc<FPType>::toDegree(test_lng) << std::endl
                          << "Test result point: " << *testLoc << "Ref point: " << *refLoc << std::endl;
            }
        }
        std::cout << loc_idx << " locations test" << std::endl
                  << "Num of errors: " << errot_count << std::endl
                  << "Error percentage is: "
                  << errot_count*100.0/loc_idx << "%" << std::endl
                  << "Error in total hav dist diff is: " << std::fabs(test_dist_err_total - ref_dist_for_err_pts_total) << std::endl;
    }
}


template <typename FPType>
void TimeBuild(std::span<const SBLoc<FPType>> loc_data_span, std::unique_ptr<SBSolver<FPType>>& solver, std::uint8_t num_builds = 1) {
    std::chrono::duration<FPType> total_elapsed_time{0};
    for (std::uint8_t ui = num_builds; ui > 0; --ui) {
        //solver = std::make_unique<BKDTSBSolver<KDTreeExpandLongestVec, FPType>>();
        std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
        solver->Build(loc_data_span);
        std::chrono::duration<FPType> elapsed_time_in_secs = std::chrono::high_resolution_clock::now() - start;
        total_elapsed_time += elapsed_time_in_secs;
        std::cout << std::regex_replace(typeid(*solver).name(), std::regex("[A-Z]?[0-9]+|.$"), "")
            << std::endl << "Build time: " << elapsed_time_in_secs.count() << std::endl;
        solver->PrintSolverInfo();
    }
    std::cout << "-----------------------------\n";
    std::cout << std::regex_replace(typeid(*solver).name(), std::regex("[A-Z]?[0-9]+|.$"), "")
        << std::endl << "Ave build time: " << total_elapsed_time.count()/num_builds << std::endl;
}

template <typename FPType>
volatile SBLoc<FPType>* volatile nn_result{}; // ensures a side effect

template <typename FPType>
void TimeNNSearch(const SBSolver<FPType> &solver, std::vector<PointND<FPType, 2>> &test_lat_lng_pts,
                  std::chrono::duration<FPType> search_duration_in_secs, absl::BitGen &bitgen) {
    std::shuffle(test_lat_lng_pts.begin(), test_lat_lng_pts.end(), bitgen);
    std::vector<std::chrono::duration<FPType, std::micro>> per_search_time_vec_in_micro_secs;
    std::chrono::duration<FPType, std::micro> total_elapsed_time_in_micro_secs, max_one_search_time_in_micro_secs{}, min_one_search_time_in_micro_secs(std::numeric_limits<FPType>::max());
    per_search_time_vec_in_micro_secs.reserve(test_lat_lng_pts.size());
    std::chrono::time_point<std::chrono::high_resolution_clock> startTime = std::chrono::high_resolution_clock::now();
    for (std::size_t loc_idx = 0;
         loc_idx < test_lat_lng_pts.size() && (total_elapsed_time_in_micro_secs = std::chrono::high_resolution_clock::now() - startTime) < search_duration_in_secs;
         ++loc_idx) {
        std::chrono::time_point<std::chrono::high_resolution_clock> this_search_start_time = std::chrono::high_resolution_clock::now();
        nn_result<FPType> = const_cast<SBLoc<FPType>*>(solver.FindNearestLoc(test_lat_lng_pts[loc_idx]));
        per_search_time_vec_in_micro_secs.emplace_back(std::chrono::high_resolution_clock::now() - this_search_start_time);
    }
    
    std::for_each(per_search_time_vec_in_micro_secs.crbegin(), per_search_time_vec_in_micro_secs.crend(), [&](auto this_search_time_in_micro_secs) mutable {
        max_one_search_time_in_micro_secs = std::max(max_one_search_time_in_micro_secs, this_search_time_in_micro_secs);
        min_one_search_time_in_micro_secs = std::min(min_one_search_time_in_micro_secs, this_search_time_in_micro_secs);
    });
    
    
    std::cout << std::endl << "Ave per Search Time: "
              << total_elapsed_time_in_micro_secs.count()/per_search_time_vec_in_micro_secs.size()
              << "us from " << per_search_time_vec_in_micro_secs.size() << " searches" << std::endl
              << "Max sinlge search time: " << max_one_search_time_in_micro_secs.count() << std::endl
              << "Min single search time: " << min_one_search_time_in_micro_secs.count() << std::endl << std::endl;
    
}

template <typename FPType>
void WriteResults(const char* argv[],
                  const std::vector<PointND<FPType, 2>> &test_lat_lng_pts,
                  const SBSolver<FPType> *solver) {
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

template<typename FPType>
std::vector<SBLoc<FPType>> ConstructLocDataVec(const char* loc_file_path, absl::BitGen &bitgen) {
    std::ifstream infilestream_locs(loc_file_path);
    std::vector<SBLoc<FPType>> loc_data_vec;
    //loc_data_vec.reserve(1 << 24); 
    loc_data_vec.reserve(11569609); 
    loc_data_vec.assign(std::istream_iterator<SBLoc<FPType>>(infilestream_locs),
        std::istream_iterator<SBLoc<FPType>>());
    assert(loc_data_vec.size() != 0);
    //std::sort(locData->begin(), locData->end());
    //locData->erase(std::unique(locData->begin(), locData->end()), locData->end());
    loc_data_vec.shrink_to_fit();
    std::shuffle(loc_data_vec.rbegin(), loc_data_vec.rend(), bitgen);
    return loc_data_vec;
}

template <typename FPType>
void MainContent(int argc, const char * argv[]) {
    FPType ave_actual_locs_per_cell = argc < 4 ? def::kDefalutAveActualLocsPerCell : std::stod(argv[3]);
    std::size_t max_cached_cell_vec_size = argc < 5 ? def::kDefaultMaxCacheCellVecSize : std::stoi(argv[4]);
    bool to_test_accuracy = argc < 6 ? def::kDefaultToTestAccuracy : std::tolower(argv[5][0]) == 'y';
    bool to_test_search_time = argc < 7 ? def::kDefaultToTestSearchTime : std::tolower(argv[6][0]) == 'y';
    std::chrono::duration<FPType> accuracy_test_time_in_secs(def::kDefaultAccuracyTestDurationInSecs);
    std::chrono::duration<FPType> search_benchmark_duration_in_secs(def::kDefaultSearchBenchDurationInSecs);
    if (to_test_accuracy && argc >= 5) {
        accuracy_test_time_in_secs = std::chrono::duration<FPType>(std::stoi(std::regex_replace(argv[5], std::regex("[^0-9]*([0-9]+).*"), std::string("$1"))));
    }
    if (to_test_search_time && argc >= 6) {
        search_benchmark_duration_in_secs = std::chrono::duration<FPType>(std::stoi(std::regex_replace(argv[6],  std::regex("[^0-9]*([0-9]+).*"), std::string("$1"))));
    }
    std::ifstream inRefLatLngPtLocPairVec;
    if (argc >= 7)
        inRefLatLngPtLocPairVec.open(argv[7]);
    std::ofstream outRefLatLngPtLocPairVec;
    if (argc >= 8)
        outRefLatLngPtLocPairVec.open(argv[8]);
    
    
    to_test_accuracy = true;
    to_test_search_time = true;
    //maxCacheCellVecSize = (1 << 16ull);
    //maxCacheCellVecSize = (1 << 9ull);
    //ave_actual_locs_per_cell = 10;
    
    absl::BitGen bitgen;
    const std::vector<SBLoc<FPType>> loc_data_vec = ConstructLocDataVec<FPType>(argv[1], bitgen);
    
    /*
    if (numOfLocsToWriteToFile) {
        auto solver = std::make_unique<BFSBSolver<FPType>>();
        //auto solver = std::make_unique<BKDTSBSolver<FPType><KDTree>>();
        timeBuild<FPType>(locData, *solver);
        writeResults<FPType>(argv, generatetestSearchLatLngPts<FPType>(numOfLocsToWriteToFile, mt), solver.get());
        return 0;
    } */
    
    

    std::unique_ptr<SBSolver<FPType>> solvers[] = {
        //std::make_unique<BKDTSBSolver<KDTreeExpandLongestVec, FPType>>(),
        //std::make_unique<BFSBSolver<FPType>>(),
        //std::make_unique<BFEUCPtSBSolver<FPType>>(),
        //std::make_unique<KDTSBSolver<KDTree,FPType>>(),
        //std::make_unique<BKDTSBSolver<KDTree, FPType>>(),
        //std::make_unique<BKDTSBSolver<FPType><KDTreeCusMem>>(),
        std::make_unique<BKDTSBSolver<KDTreeExpandLongest, FPType>>(),
        //std::make_unique<UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongest, FPType, def::ThreadingPolicy::kSingle>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
        //std::make_unique<UnionUniCellBKDTGridSBSolver<KDTreeExpandLongest, FPType, def::ThreadingPolicy::kSingle>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
        std::make_unique<BKDTSBSolver<KDTreeExpandLongestVec, FPType>>(),
        //std::make_unique<GridSBSolver<FPType>>(),
        //std::make_unique<BKDTGridSBSolver<FPType>>(aveLocPerCell),
        //std::make_unique<UniLatLngBKDTGridSBSolver<KDTree, FPType>>(0.85*aveLocPerCell, kMaxCacheCellVecSize_),
        //std::make_unique<UniLatLngBKDTGridSBSolver<KDTreeCusMem, FPType>>(0.85*aveLocPerCell, kMaxCacheCellVecSize_),
        //std::make_unique<UniLatLngBKDTGridSBSolver<KDTreeExpandLongest, FPType>>(0.85*aveLocPerCell, kMaxCacheCellVecSize_),
        //std::make_unique<UniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec,FPType>>(aveLocPerCell, kMaxCacheCellVecSize_),
        //std::make_unique<UnionUniLatLngBKDTGridSBSolver<FPType><KDTreeExpandLongestVec>>(aveLocPerCell, kMaxCacheCellVecSize_),
        //std::make_unique<UniCellBKDTGridSBSolver<FPType><KDTree>>(aveLocPerCell, maxCacheCellVecSize),
        // std::make_unique<UniCellBKDTGridSBSolver<FPType><KDTreeCusMem>>(aveLocPerCell, maxCacheCellVecSize),
        //std::make_unique<UniCellBKDTGridSBSolver<KDTreeExpandLongest, FPType>>(aveLocPerCell, kMaxCacheCellVecSize_),
        //std::make_unique<UniCellBKDTGridSBSolver<KDTreeExpandLongestVec, FPType>>(aveLocPerCell, kMaxCacheCellVecSize_),
        //std::make_unique<UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, FPType, def::ThreadingPolicy::kSingle>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
        //std::make_unique<UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, FPType, def::ThreadingPolicy::kSingle>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
        //std::make_unique<UnionUniLatLngBKDTGridSBSolver<KDTreeExpandLongestVec, FPType, def::ThreadingPolicy::kMultiOmp>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
        //std::make_unique<UnionUniCellBKDTGridSBSolver<KDTreeExpandLongestVec, FPType, def::ThreadingPolicy::kMultiOmp>>(ave_actual_locs_per_cell, max_cached_cell_vec_size),
    };
    
    std::vector<PointND<FPType, 2>> search_bench_test_lat_lng_pts;
    if (to_test_search_time)
        search_bench_test_lat_lng_pts = GenerateTestLatLngPts<FPType>(def::kMaxTestLocs, bitgen);
    std::vector<const SBLoc<FPType>*> ref_locs;
    std::vector<PointND<FPType, 2>> accuracy_test_lat_lng_pts;
    if (to_test_accuracy)
        accuracy_test_lat_lng_pts = GenerateTestLatLngPts<FPType>(def::kMaxTestLocs, bitgen);

    for (std::size_t i = 0; i < std::size(solvers); ++i) {
        TimeBuild(std::forward<std::span<const SBLoc<FPType>>>(loc_data_vec), solvers[i]);
        if (to_test_search_time)
            TimeNNSearch(*solvers[i], search_bench_test_lat_lng_pts, search_benchmark_duration_in_secs, bitgen);
        if (to_test_accuracy)
            AccuracyTestFromRefSolver(accuracy_test_lat_lng_pts, ref_locs, *solvers[i], accuracy_test_time_in_secs);
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


/*
#include <thread>

int main(int argc, char** argv)
{
    std::chrono::time_point<std::chrono::high_resolution_clock> start = std::chrono::high_resolution_clock::now();
    //for (volatile int i = 0; i < 10; i++)
        std::thread([]() {}).detach();
    std::chrono::duration<double> elapsed_time_in_secs = std::chrono::high_resolution_clock::now() - start;
    std::cout << "time: " << elapsed_time_in_secs.count() << std::endl;
    return 0;
}*/
 
 
 
 
