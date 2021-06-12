cmake_minimum_required(VERSION 3.0.0)
project(274NN666 VERSION 0.1.0)

include(CTest)
enable_testing()

set(CMAKE_CXX_FLAGS_RELEASE "-Ofast -ffast-math")
set(CMAKE_CXX_FLAGS "-fopenmp")


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(
        ${PROJECT_SOURCE_DIR}/cppsrc
        SYSTEM /usr/local/include
)


file(GLOB all_SRCS
        "${PROJECT_SOURCE_DIR}/cppsrc/*.h"
        "${PROJECT_SOURCE_DIR}/cppsrc/*.hpp"
        "${PROJECT_SOURCE_DIR}/cppsrc/*.cpp"
        "${PROJECT_SOURCE_DIR}/cppsrc/*.c"
        )

add_executable(274NN666 ${all_SRCS})

project(274NN666)



set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
