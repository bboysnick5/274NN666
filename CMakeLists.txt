set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -DNDEBUG")
cmake_minimum_required(VERSION 3.23)
project(274NN666 VERSION 0.1.0)

set(OpenMP_CXX_FLAGS "-fopenmp")
set(OpenMP_CXX_LIB_NAMES "libomp")

set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

include_directories(
        ${PROJECT_SOURCE_DIR}/cppsrc
        SYSTEM /usr/local/include
)

link_directories(/usr/local/lib/)

find_package(absl REQUIRED)
find_package(OpenMP)
find_package(Vc)

file(GLOB all_SRCS
        "${PROJECT_SOURCE_DIR}/cppsrc/*.h"
        "${PROJECT_SOURCE_DIR}/cppsrc/*.hpp"
        "${PROJECT_SOURCE_DIR}/cppsrc/*.cpp"
        "${PROJECT_SOURCE_DIR}/cppsrc/*.c"
)

add_executable(274NN666 ${all_SRCS})
target_link_libraries(274NN666 PUBLIC Vc::Vc)
target_link_libraries(274NN666 PUBLIC absl::random_random)
target_compile_features(274NN666 PRIVATE cxx_std_23)

if(OpenMP_CXX_FOUND)
        target_link_libraries(274NN666 PUBLIC OpenMP::OpenMP_CXX)
endif()

project(274NN666)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

message(${CMAKE_CXX_FLAGS_RELEASE})