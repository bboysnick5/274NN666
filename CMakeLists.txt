set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -march=native -DNDEBUG")
cmake_minimum_required(VERSION 3.12)
project(274NN666 VERSION 0.1.0)

include(CTest)
enable_testing()

set(OpenMP_CXX_FLAGS "-fopenmp")
set(OpenMP_CXX_LIB_NAMES "libomp")


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(
        ${PROJECT_SOURCE_DIR}/cppsrc
        SYSTEM /usr/local/include
)


link_directories( /usr/local/lib/)

find_package(absl REQUIRED)
find_package(OpenMP)


file(GLOB all_SRCS
        "${PROJECT_SOURCE_DIR}/cppsrc/*.h"
        "${PROJECT_SOURCE_DIR}/cppsrc/*.hpp"
        "${PROJECT_SOURCE_DIR}/cppsrc/*.cpp"
        "${PROJECT_SOURCE_DIR}/cppsrc/*.c"
        )

add_executable(274NN666 ${all_SRCS})
target_link_libraries(274NN666 PUBLIC absl::random_random)
if(OpenMP_CXX_FOUND)
    target_link_libraries(274NN666 PUBLIC OpenMP::OpenMP_CXX)
endif()

project(274NN666)



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)

message(${CMAKE_CXX_FLAGS_RELEASE})