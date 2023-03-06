/*
 * @Author: Nick Liu
 * @Date: 2022-08-02 21:32:38
 * @LastEditTime: 2022-09-07 15:58:32
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/SolverFactory.hpp
 */

#include "BFCartPtSolver.hpp"
#include "BFLocalStorageSBSolver.hpp"
#include "BFSBSolver.hpp"
#include "BKDTGridSBSolver.hpp"
#include "BKDTSBSolver.hpp"
#include "Definition.hpp"
#include "GridSBSolver.hpp"
#include "KDTSBSolver.hpp"
#include "SBSolver.hpp"
#include "UniCellBKDTGridSBSolver.hpp"
#include "UniLatLngBKDTGridSBSolver.hpp"
#include "UnionUniCellBKDTGridSBSolver.hpp"
#include "UnionUniLatLngBKDTGridSBSolver.hpp"

// #include <benchmark/benchmark.h>
#include <absl/random/random.h>

#include <algorithm>
#include <boost/hana.hpp>
#include <cassert>
#include <chrono>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <random>
#include <regex>
#include <span>
#include <sstream>
#include <string>
#include <thread>
#include <type_traits>
#include <typeinfo>
#include <utility>
#include <vector>

enum class SolverName {
    kSolverBf,
    kSolverBfCartPt,
    kSolverBfLocalStorage,
    kSolverKdt,
    kSolverBkdt,
    kSolverUniCell,
    kSolverUnionUniCell,
    kSolverUniLatLng,
    kSolverUnionUniLatLng,
    kSolverGrid,
};

template <typename T>
struct Deleter {
    void operator()(T* t) { delete t; }
};

template <std::floating_point FPType>
using SolverRegistry = std::variant<
    std::unique_ptr<SBSolver<kConfigStSoaNoSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigStSoaNoSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpSoaNoSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpSoaNoSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosNoSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigStAosNoSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpAosNoSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpAosNoSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosoaNoSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigStAosoaNoSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpAosoaNoSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpAosoaNoSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStSoaSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigStSoaSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpSoaSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpSoaSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigStAosSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpAosSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpAosSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosoaSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigStAosoaSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpAosoaSimdNoTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpAosoaSimdNoTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStSoaNoSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigStSoaNoSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpSoaNoSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpSoaNoSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosNoSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigStAosNoSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpAosNoSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpAosNoSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosoaNoSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigStAosoaNoSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpAosoaNoSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpAosoaNoSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStSoaSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigStSoaSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpSoaSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpSoaSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigStAosSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpAosSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpAosSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosoaSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigStAosoaSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpAosoaSimdKDTree<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpAosoaSimdKDTree<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStSoaNoSimdKDTreeCusMem<FPType>>,
                    Deleter<SBSolver<kConfigStSoaNoSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpSoaNoSimdKDTreeCusMem<FPType>>,
        Deleter<SBSolver<kConfigMtOmpSoaNoSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosNoSimdKDTreeCusMem<FPType>>,
                    Deleter<SBSolver<kConfigStAosNoSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosNoSimdKDTreeCusMem<FPType>>,
        Deleter<SBSolver<kConfigMtOmpAosNoSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStAosoaNoSimdKDTreeCusMem<FPType>>,
        Deleter<SBSolver<kConfigStAosoaNoSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosoaNoSimdKDTreeCusMem<FPType>>,
        Deleter<SBSolver<kConfigMtOmpAosoaNoSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStSoaSimdKDTreeCusMem<FPType>>,
                    Deleter<SBSolver<kConfigStSoaSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpSoaSimdKDTreeCusMem<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpSoaSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosSimdKDTreeCusMem<FPType>>,
                    Deleter<SBSolver<kConfigStAosSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigMtOmpAosSimdKDTreeCusMem<FPType>>,
                    Deleter<SBSolver<kConfigMtOmpAosSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<SBSolver<kConfigStAosoaSimdKDTreeCusMem<FPType>>,
                    Deleter<SBSolver<kConfigStAosoaSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosoaSimdKDTreeCusMem<FPType>>,
        Deleter<SBSolver<kConfigMtOmpAosoaSimdKDTreeCusMem<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStSoaNoSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigStSoaNoSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStAosNoSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigStAosNoSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStAosoaNoSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigStAosoaNoSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigMtOmpAosoaNoSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStSoaSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigStSoaSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStAosSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigStAosSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigMtOmpAosSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStAosoaSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigStAosoaSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongest<FPType>>,
        Deleter<SBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongest<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStSoaNoSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigStSoaNoSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStAosNoSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigStAosNoSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigMtOmpAosNoSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStAosoaNoSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigStAosoaNoSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<
            SBSolver<kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStSoaSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigStSoaSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigMtOmpSoaSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStAosSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigStAosSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigMtOmpAosSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigStAosoaSimdKDTreeExpandLongestVec<FPType>>,
        Deleter<SBSolver<kConfigStAosoaSimdKDTreeExpandLongestVec<FPType>>>>,
    std::unique_ptr<
        SBSolver<kConfigMtOmpAosoaSimdKDTreeExpandLongestVec<FPType>>>>;

template <std::floating_point FPType>
class SolverFactory {
   public:
    constexpr SolverFactory(auto ave_actual_locs_per_cell,
                            auto max_cached_cell_vec_size)
        : ave_actual_locs_per_cell_(ave_actual_locs_per_cell),
          max_cached_cell_vec_size_(max_cached_cell_vec_size) {}

    template <SolverConfig Config>
    constexpr std::unique_ptr<SBSolver<Config>> Create(SolverName) const;

   private:
    FPType ave_actual_locs_per_cell_;
    std::size_t max_cached_cell_vec_size_;
};

template <std::floating_point FPType>
template <SolverConfig Config>
constexpr std::unique_ptr<SBSolver<Config>> SolverFactory<FPType>::Create(
    SolverName name) const {
    switch (name) {
        case SolverName::kSolverBf:
            return std::make_unique<BFSBSolver<Config>>();
        case SolverName::kSolverBfCartPt:
            return std::make_unique<BFCartPtSBSolver<Config>>();
        case SolverName::kSolverKdt:
            return std::make_unique<KDTSBSolver<Config>>();
        case SolverName::kSolverBkdt:
            return std::make_unique<BKDTSBSolver<Config>>();
        case SolverName::kSolverBfLocalStorage:
            return std::make_unique<BFLocalStorageSBSolver<Config>>();
        case SolverName::kSolverUniCell:
            return std::make_unique<UniCellBKDTGridSBSolver<Config>>(
                ave_actual_locs_per_cell_, max_cached_cell_vec_size_);
        case SolverName::kSolverUnionUniCell:
            return std::make_unique<UnionUniCellBKDTGridSBSolver<Config>>(
                ave_actual_locs_per_cell_, max_cached_cell_vec_size_);
        case SolverName::kSolverUniLatLng:
            return std::make_unique<UniLatLngBKDTGridSBSolver<Config>>(
                ave_actual_locs_per_cell_, max_cached_cell_vec_size_);
        case SolverName::kSolverUnionUniLatLng:
            return std::make_unique<UnionUniLatLngBKDTGridSBSolver<Config>>(
                ave_actual_locs_per_cell_, max_cached_cell_vec_size_);
        case SolverName::kSolverGrid:
            return std::make_unique<GridSBSolver<Config>>();
    }
}