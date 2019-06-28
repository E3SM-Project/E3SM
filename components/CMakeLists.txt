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

include(${CMAKE_CURRENT_SOURCE_DIR}/common_setup.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/build_model.cmake)

# For each component that has registered itself to be built, read in Filepath and
# CCSM_cppdefs and build it. We identify the components by looking for Buildconf/*/Filepath
# files.

set(BUILDCONF ${CMAKE_BINARY_DIR}/../../Buildconf)
file(GLOB FILEPATHS "${BUILDCONF}/*/Filepath")
foreach(ITEM IN LISTS FILEPATHS)
  get_filename_component(FILEPATH_DIR ${ITEM} DIRECTORY)
  get_filename_component(MODELCONF ${FILEPATH_DIR} NAME)
  string(REPLACE "conf" "" MODEL "${MODELCONF}")

  message("Found component model ${MODEL}")
  list(APPEND MODELS ${MODEL})

  build_model(${MODEL} ${MODELCONF} ${ITEM})
endforeach()

set(EXEC_SE "${CIME_MODEL}.exe")
add_executable(${EXEC_SE})
set(ALL_LIBS "${ULIBDEP} ${MCTLIBS} ${PIOLIB} ${GPTLLIB} ${SLIBS} ${MLIBS} ${F90_LDFLAGS}")
separate_arguments(ALL_LIBS_LIST UNIX_COMMAND "${ALL_LIBS}")
foreach(ITEM IN LISTS ALL_LIBS_LIST)
  target_link_libraries(${EXEC_SE} ${ITEM})
endforeach()

