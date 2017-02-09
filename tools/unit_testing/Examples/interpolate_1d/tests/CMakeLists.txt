cmake_minimum_required(VERSION 2.8)
project(interpolate_1d_tests Fortran C)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_DIRECTORY})
include(CIME_utils)

add_subdirectory(../src interpolate_1d)

add_subdirectory(CTest)

add_subdirectory(pFUnit)
