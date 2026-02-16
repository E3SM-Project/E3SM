cmake_minimum_required(VERSION 3.16)
project(emulators LANGUAGES CXX)

# Path to externals (relative to this file)
get_filename_component(E3SM_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(EXTERNALS_ROOT "${E3SM_ROOT}/externals")

option(BUILD_EMULATOR_TESTS "Build emulator unit tests" OFF)

# When built standalone (i.e. this is the top-level project), default to
# standalone mode so we don't try to compile Fortran sources (the project
# only declares CXX).  When included from the full E3SM build Fortran
# is already enabled by the parent project.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(STANDALONE_MODE ON CACHE BOOL "Standalone build (CXX only, no Fortran)" FORCE)
endif()

if(BUILD_EMULATOR_TESTS)
  # Catch2 from ekat
  set(CATCH2_INCLUDE_DIR "${EXTERNALS_ROOT}/ekat/extern/Catch2/single_include")

  # Verify Catch2 exists
  if(NOT EXISTS "${CATCH2_INCLUDE_DIR}/catch2/catch.hpp")
    message(FATAL_ERROR "Catch2 not found at ${CATCH2_INCLUDE_DIR}. "
      "Make sure externals/ekat is available.")
  endif()

  # Enable testing at the top level
  enable_testing()
endif()

# Add subdirectories
add_subdirectory(common)
add_subdirectory(eatm)
