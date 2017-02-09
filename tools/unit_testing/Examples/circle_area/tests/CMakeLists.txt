cmake_minimum_required(VERSION 2.8)
project(circle_area_tests Fortran C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_DIRECTORY})
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIRECTORY}/externals/pio2/cmake)
include(CIME_utils)

add_subdirectory(../src circle)

add_subdirectory(CTest)

add_subdirectory(pFUnit)
