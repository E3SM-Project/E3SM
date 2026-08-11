cmake_minimum_required(VERSION 3.18)
project(emulators LANGUAGES CXX Fortran)

include(FetchContent)
option(BUILD_EMULATOR_TESTS "Build emulator unit tests" OFF)

# ------------------------------------------------------------------------------
# EKAT
# ------------------------------------------------------------------------------

# Build only EKAT Core + YAML parser.
#
# Core is always enabled. All other optional EKAT packages default to OFF,
# so we only need to explicitly enable the YAML parser.
set(EKAT_ENABLE_ALL_PACKAGES OFF CACHE BOOL "" FORCE)
set(EKAT_ENABLE_YAML_PARSER ON CACHE BOOL "" FORCE)

set(EKAT_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
set(EKAT_ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
set(EKAT_ENABLE_MPI OFF CACHE BOOL "" FORCE) # unexpected ekat is not parser mpi.h correctly?

FetchContent_Declare(
    ekat
    GIT_REPOSITORY https://github.com/E3SM-Project/EKAT.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(ekat)

if(BUILD_EMULATOR_TESTS)
  # Enable testing at the top level
  enable_testing()
endif()



# Add subdirectories
add_subdirectory(common)
add_subdirectory(emulatoratm)
add_subdirectory(driver)
