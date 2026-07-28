cmake_minimum_required(VERSION 3.16)
project(emulators LANGUAGES CXX Fortran)

# Path to externals (relative to this file)
get_filename_component(E3SM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(EXTERNALS_ROOT "${E3SM_ROOT}/externals")

option(BUILD_EMULATOR_TESTS "Build emulator unit tests" OFF)

# MPI is what turns InferenceContext from a serial stub into a real
# description of the component's ranks.  Without it everything still builds
# and the unit tests still run; the context just reports rank 0 of 1.
set(EMULATOR_ENABLE_MPI AUTO CACHE STRING
  "Use MPI in the inference context: AUTO, ON or OFF")

if(EMULATOR_ENABLE_MPI STREQUAL "AUTO")
  find_package(MPI QUIET COMPONENTS CXX)
  set(EMULATOR_HAVE_MPI ${MPI_CXX_FOUND})
elseif(EMULATOR_ENABLE_MPI)
  find_package(MPI REQUIRED COMPONENTS CXX)
  set(EMULATOR_HAVE_MPI TRUE)
else()
  set(EMULATOR_HAVE_MPI FALSE)
endif()
message(STATUS "Emulator inference: MPI=${EMULATOR_HAVE_MPI}")

# When built standalone (i.e. this is the top-level project), default to
# standalone mode so we don't try to compile Fortran sources (the project
# only declares CXX).  When included from the full E3SM build Fortran
# is already enabled by the parent project.
# if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
#   set(STANDALONE_MODE ON CACHE BOOL "Standalone build (CXX only, no Fortran)" FORCE)
# endif()

if(BUILD_EMULATOR_TESTS)
  # Catch2, as the v2 single header the tests include.  ekat used to vendor
  # this as a submodule and now pulls it with FetchContent, so the old path
  # is empty on a fresh clone; prefer a local copy wherever one exists and
  # otherwise fetch the one header we actually use.  Set CATCH2_INCLUDE_DIR
  # to point at a directory containing catch2/catch.hpp on a machine with no
  # network access from the build node.
  if(NOT CATCH2_INCLUDE_DIR)
    if(EXISTS "${EXTERNALS_ROOT}/ekat/extern/Catch2/single_include/catch2/catch.hpp")
      set(CATCH2_INCLUDE_DIR
        "${EXTERNALS_ROOT}/ekat/extern/Catch2/single_include")
    else()
      set(CATCH2_INCLUDE_DIR "${CMAKE_BINARY_DIR}/catch2")
      if(NOT EXISTS "${CATCH2_INCLUDE_DIR}/catch2/catch.hpp")
        message(STATUS "Emulator tests: downloading Catch2 v2.13.10")
        file(DOWNLOAD
          "https://raw.githubusercontent.com/catchorg/Catch2/v2.13.10/single_include/catch2/catch.hpp"
          "${CATCH2_INCLUDE_DIR}/catch2/catch.hpp"
          EXPECTED_HASH
            SHA256=3725c0f0a75f376a5005dde31ead0feb8f7da7507644c201b814443de8355170
          TLS_VERIFY ON
          STATUS EMULATOR_CATCH2_DOWNLOAD)
        list(GET EMULATOR_CATCH2_DOWNLOAD 0 EMULATOR_CATCH2_FAILED)
        if(EMULATOR_CATCH2_FAILED)
          list(GET EMULATOR_CATCH2_DOWNLOAD 1 EMULATOR_CATCH2_ERROR)
          message(FATAL_ERROR "Could not download Catch2: "
            "${EMULATOR_CATCH2_ERROR}.  Point CATCH2_INCLUDE_DIR at a "
            "directory containing catch2/catch.hpp instead.")
        endif()
      endif()
    endif()
  endif()

  if(NOT EXISTS "${CATCH2_INCLUDE_DIR}/catch2/catch.hpp")
    message(FATAL_ERROR
      "No catch2/catch.hpp under CATCH2_INCLUDE_DIR (${CATCH2_INCLUDE_DIR}).")
  endif()

  # Enable testing at the top level
  enable_testing()
endif()

# Add subdirectories
add_subdirectory(common)
add_subdirectory(emulatoratm)
add_subdirectory(driver)
