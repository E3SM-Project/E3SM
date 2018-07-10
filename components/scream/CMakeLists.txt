cmake_minimum_required(VERSION 2.8.11)
project(SCREAM CXX)

include(CTest)

set(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Select build type.")

set(CATCH_INCL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/catch2/include)

enable_testing()

add_subdirectory(coupler)
add_subdirectory(p3)
add_subdirectory(rrtmgp)
add_subdirectory(shoc)
add_subdirectory(tests)
