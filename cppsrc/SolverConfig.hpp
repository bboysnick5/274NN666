
#include "Definition.hpp"
#include "KDTree.hpp"
#include "KDTreeCusMem.hpp"
#include "KDTreeExpandLongest.hpp"
#include "KDTreeExpandLongestVec.hpp"
#include "SBLoc.hpp"

template <std::floating_point FPTypeT, def::ParPolicy ParPolicy,
          template <class FP, std::uint8_t N, class Elem, def::DistType DT>
          class KDTTypeT = KDTree>
struct SolverConfig {
    using KDTType =
        KDTTypeT<FPTypeT, 3, const SBLoc<FPTypeT>*, def::DistType::kEuc>;
    using FPType = FPTypeT;
    def::ParPolicy par_policy = ParPolicy;
};

// concept SolverConfigConcept

template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaNoSimdNoTree =
    SolverConfig<FPType, def::kParPolicyStSoaNoSimd, NoTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaNoSimdKDTree =
    SolverConfig<FPType, def::kParPolicyStSoaNoSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaNoSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyStSoaNoSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaNoSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyStSoaNoSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaNoSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyStSoaNoSimd, KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaNoSimdNoTree =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaNoSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaNoSimdKDTree =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaNoSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaNoSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaNoSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaNoSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaNoSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaNoSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaNoSimd,
                 KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosNoSimdNoTree =
    SolverConfig<FPType, def::kParPolicyStAosNoSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosNoSimdKDTree =
    SolverConfig<FPType, def::kParPolicyStAosNoSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosNoSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyStAosNoSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosNoSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyStAosNoSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosNoSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyStAosNoSimd, KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosNoSimdNoTree =
    SolverConfig<FPType, def::kParPolicyMtOmpAosNoSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosNoSimdKDTree =
    SolverConfig<FPType, def::kParPolicyMtOmpAosNoSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosNoSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyMtOmpAosNoSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosNoSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyMtOmpAosNoSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosNoSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyMtOmpAosNoSimd,
                 KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaNoSimdNoTree =
    SolverConfig<FPType, def::kParPolicyStAosoaNoSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaNoSimdKDTree =
    SolverConfig<FPType, def::kParPolicyStAosoaNoSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaNoSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyStAosoaNoSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaNoSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyStAosoaNoSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaNoSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyStAosoaNoSimd,
                 KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaNoSimdNoTree =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaNoSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaNoSimdKDTree =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaNoSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaNoSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaNoSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaNoSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaNoSimd,
                 KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaNoSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaNoSimd,
                 KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaSimdNoTree =
    SolverConfig<FPType, def::kParPolicyStSoaSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaSimdKDTree =
    SolverConfig<FPType, def::kParPolicyStSoaSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyStSoaSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyStSoaSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStSoaSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyStSoaSimd, KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaSimdNoTree =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaSimdKDTree =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpSoaSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyMtOmpSoaSimd, KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosSimdNoTree =
    SolverConfig<FPType, def::kParPolicyStAosSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosSimdKDTree =
    SolverConfig<FPType, def::kParPolicyStAosSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyStAosSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyStAosSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyStAosSimd, KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosSimdNoTree =
    SolverConfig<FPType, def::kParPolicyMtOmpAosSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosSimdKDTree =
    SolverConfig<FPType, def::kParPolicyMtOmpAosSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyMtOmpAosSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyMtOmpAosSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyMtOmpAosSimd, KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaSimdNoTree =
    SolverConfig<FPType, def::kParPolicyStAosoaSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaSimdKDTree =
    SolverConfig<FPType, def::kParPolicyStAosoaSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyStAosoaSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyStAosoaSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigStAosoaSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyStAosoaSimd, KDTreeExpandLongestVec>{};

template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaSimdNoTree =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaSimd>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaSimdKDTree =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaSimd, KDTree>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaSimdKDTreeCusMem =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaSimd, KDTreeCusMem>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaSimdKDTreeExpandLongest =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaSimd, KDTreeExpandLongest>{};
template <std::floating_point FPType>
static constexpr SolverConfig kConfigMtOmpAosoaSimdKDTreeExpandLongestVec =
    SolverConfig<FPType, def::kParPolicyMtOmpAosoaSimd,
                 KDTreeExpandLongestVec>{};
