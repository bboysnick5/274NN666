/*
 * @Author: Nick Liu
 * @Date: 2021-05-25 19:55:00
 * @LastEditTime: 2022-08-11 11:36:15
 * @LastEditors: Nick Liu
 * @Description:
 * @FilePath: /274F201666/cppsrc/Definition.hpp
 */

#ifndef Definition_hpp
#define Definition_hpp

#include <stdio.h>

#include <concepts>
#include <iterator>
#include <numbers>
#include <tuple>
#include <type_traits>
#include <utility>

namespace def {

using kDefaultDistType = double;
inline constexpr std::size_t kDefaultMaxCacheCellVecSize = 896;
// to better fit into 32kb L1 with each vec element being 3*sizeof(double) +
// pointer

inline constexpr double kDefalutAveActualLocsPerCell = 0.4;

inline constexpr bool kDefaultToTestAccuracy = false;
inline constexpr bool kDefaultToTestSearchTime = true;
inline constexpr std::size_t kDefaultAccuracyTestDurationInSecs = 10;
inline constexpr std::size_t kDefaultSearchBenchDurationInSecs = 10;
inline constexpr std::size_t kMaxTestLocs = 1 << 24;

template <class ConstIt>
concept const_iterator = std::is_const_v<typename std::remove_pointer_t<
    typename std::iterator_traits<ConstIt>::pointer>>;

template <class Non_Const_It>
concept non_const_iterator = !
std::is_const_v<typename std::remove_pointer_t<
    typename std::iterator_traits<Non_Const_It>::pointer>>;

template <typename T>
constexpr auto &kMathPi = std::numbers::pi_v<T>;

enum class ThreadingPolicy {
    kSingle,
    kMultiOmp,
    kMultiHand,
};

enum class MemoryLayout {
    kAos,
    kSoA,
    kAosoa,
};

enum class ExplicitSimd { kEnable, kDisable };
struct ParPolicy {
    ThreadingPolicy thread_policy;
    MemoryLayout mem_layout;
    ExplicitSimd exp_simd;
};

constexpr ParPolicy kParPolicyStSoaNoSimd = {def::ThreadingPolicy::kSingle,
                                             def::MemoryLayout::kSoA,
                                             def::ExplicitSimd::kDisable};
constexpr ParPolicy kParPolicyStSoaSimd = {def::ThreadingPolicy::kSingle,
                                           def::MemoryLayout::kSoA,
                                           def::ExplicitSimd::kEnable};
constexpr ParPolicy kParPolicyStAosSimd = {def::ThreadingPolicy::kSingle,
                                           def::MemoryLayout::kAos,
                                           def::ExplicitSimd::kEnable};
;
constexpr ParPolicy kParPolicyStAosoaSimd = {def::ThreadingPolicy::kSingle,
                                             def::MemoryLayout::kAosoa,
                                             def::ExplicitSimd::kEnable};
constexpr ParPolicy kParPolicyMtOmpSoaNoSimd = {def::ThreadingPolicy::kMultiOmp,
                                                def::MemoryLayout::kSoA,
                                                def::ExplicitSimd::kDisable};
constexpr ParPolicy kParPolicyMtOmpSoaSimd = {def::ThreadingPolicy::kMultiOmp,
                                              def::MemoryLayout::kSoA,
                                              def::ExplicitSimd::kEnable};
constexpr ParPolicy kParPolicyMtOmpAosNoSimd = {def::ThreadingPolicy::kMultiOmp,
                                                def::MemoryLayout::kAos,
                                                def::ExplicitSimd::kDisable};
constexpr ParPolicy kParPolicyMtOmpAosSimd = {def::ThreadingPolicy::kMultiOmp,
                                              def::MemoryLayout::kAos,
                                              def::ExplicitSimd::kEnable};
constexpr ParPolicy kParPolicyMtOmpAosoaNoSimd = {
    def::ThreadingPolicy::kMultiOmp, def::MemoryLayout::kAosoa,
    def::ExplicitSimd::kDisable};
constexpr ParPolicy kParPolicyMtOmpAosoaSimd = {def::ThreadingPolicy::kMultiOmp,
                                                def::MemoryLayout::kAosoa,
                                                def::ExplicitSimd::kEnable};
constexpr ParPolicy kParPolicyMtHandSoaNoSimd = {
    def::ThreadingPolicy::kMultiHand, def::MemoryLayout::kSoA,
    def::ExplicitSimd::kDisable};

constexpr ParPolicy kParPolicyStAosNoSimd = {def::ThreadingPolicy::kSingle,
                                             def::MemoryLayout::kAos,
                                             def::ExplicitSimd::kDisable};

constexpr ParPolicy kParPolicyStAosoaNoSimd = {def::ThreadingPolicy::kSingle,
                                               def::MemoryLayout::kAosoa,
                                               def::ExplicitSimd::kDisable};

template <ThreadingPolicy = ThreadingPolicy::kSingle>
struct ThreadingPolicyTag {};

template <>
struct ThreadingPolicyTag<ThreadingPolicy::kMultiOmp> {};

template <>
struct ThreadingPolicyTag<ThreadingPolicy::kMultiHand> {};

enum class DistType { kEuc = 0, kEucSq, kMan, kHav, kHavComp };

template <DistType = DistType::kEucSq>
struct DistTypeTag {};

template <>
struct DistTypeTag<DistType::kEuc> {};

template <>
struct DistTypeTag<DistType::kMan> {};

template <>
struct DistTypeTag<DistType::kHav> {};

template <>
struct DistTypeTag<DistType::kHavComp> {};

template <class Enum, Enum Value>
    requires std::is_enum_v<Enum>
using EnumConstant = std::integral_constant<std::underlying_type_t<Enum>,
                                            std::to_underlying(Value)>;

template <def::ParPolicy ParPolicy>
struct ParPolicyConstant {
    static constexpr def::ParPolicy value = ParPolicy;
};

}  // namespace def

namespace ns {}

namespace debug {
inline static bool mis_match = false;
}

#endif /* Definition_hpp */
