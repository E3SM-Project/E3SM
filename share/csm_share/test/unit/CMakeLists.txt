cmake_minimum_required(VERSION 2.8)
project(csm_share_tests Fortran C)

list(APPEND CMAKE_MODULE_PATH ${CESM_CMAKE_MODULE_DIRECTORY})
include(CESM_utils)

set(SHAREROOT ../..)
add_subdirectory(${SHAREROOT}/shr csm_share)
include_directories(${SHAREROOT}/include)
include_directories(${SHAREROOT}/shr)

add_subdirectory(mock)

add_subdirectory(shr_assert_test)

add_subdirectory(shr_spfn_test)

add_subdirectory(shr_infnan_test)

add_subdirectory(shr_strconvert_test)

add_subdirectory(shr_log_test)

add_subdirectory(dynamic_vector)

add_subdirectory(shr_vmath_test)
