cmake_minimum_required(VERSION 2.8)

list(APPEND CMAKE_MODULE_PATH ${CIME_CMAKE_MODULE_DIRECTORY})
include(CIME_initial_setup)

project(interpolate_1d_tests Fortran C)

include(CIME_utils)

add_subdirectory(../src interpolate_1d)

add_subdirectory(CTest)

add_subdirectory(pFUnit)
