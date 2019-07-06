#===============================================================================
#
# Common CMakeLists.txt: a framework for building all CIME components and more
#
# This is a port of cime/scripts/Tools/Makefile. As more components are ported to
# CMake, the directory level of this file will rise to the top-level directory.
#
# We will prefer space-separated strings over lists
#
#===============================================================================

# bmpersch
# - Read modern cmake docs, use modern features
# - Use find_package for trilinos and other TPLS

cmake_minimum_required(VERSION 3.9)
cmake_policy(SET CMP0057 NEW)

project(E3SM C CXX Fortran)

include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_util.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/build_mpas_model.cmake)

set(BUILDCONF ${CASEROOT}/Buildconf)

# Do any MPAS builds first since they manage their own flags
build_mpas_models()

include(${CMAKE_CURRENT_SOURCE_DIR}/common_setup.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/build_model.cmake)

# For each component that has registered itself to be built, read in Filepath and
# CCSM_cppdefs and build it. We identify the components by looking for Buildconf/*/Filepath
# files.
file(GLOB FILEPATHS "${BUILDCONF}/*/Filepath")
foreach(ITEM IN LISTS FILEPATHS)
  get_filename_component(FILEPATH_DIR ${ITEM} DIRECTORY)
  get_filename_component(MODELCONF ${FILEPATH_DIR} NAME)
  string(REPLACE "conf" "" MODEL "${MODELCONF}")

  if (NOT MODEL STREQUAL "cpl")
    message("Found component model ${MODEL}")
    list(APPEND MODELS ${MODEL})

    build_model(${MODEL} ${FILEPATH_DIR} "")
  endif()
endforeach()

message("Configuring exe with cpl and models ${MODELS}")
build_model("cpl" "${BUILDCONF}/cplconf" "${MODELS}")

