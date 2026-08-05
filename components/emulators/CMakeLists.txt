cmake_minimum_required(VERSION 3.18)
project(emulators LANGUAGES C CXX Fortran)

include(FetchContent)
set(
    FETCHCONTENT_BASE_DIR
    "${CMAKE_SOURCE_DIR}/.cmake-deps"
    CACHE PATH "FetchContent dependency cache"
    FORCE
)
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

if (STANDALONE_MODE)
  find_package(MPI COMPONENTS C)
  if (NOT MPI_FOUND)
    message(FATAL_ERROR "MPI is required for Standalone Mode, but it couldn't be found!")
  endif()
  link_libraries(MPI::MPI_C)

  # Fetch MOAB so we don't have to have it preinstalled
  add_library(MOAB STATIC IMPORTED GLOBAL)
  set_target_properties(MOAB PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libMOAB.a)
  include(ExternalProject)
  set(MOAB_OPTIONS -G ${CMAKE_GENERATOR})
  set(MOAB_OPTIONS ${MOAB_OPTIONS} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR})
  set(MOAB_OPTIONS ${MOAB_OPTIONS} -DBUILD_SHARED_LIBS=OFF)
  set(MOAB_OPTIONS ${MOAB_OPTIONS} -DENABLE_MPI=ON)
  set(MOAB_OPTIONS ${MOAB_OPTIONS} -DENABLE_FORTRAN=OFF)
  set(MOAB_OPTIONS ${MOAB_OPTIONS} -DENABLE_PYMOAB=OFF)
  set(MOAB_OPTIONS ${MOAB_OPTIONS} -DENABLE_TESTING=OFF)
  ExternalProject_Add(MOAB_proj
    URL https://bitbucket.org/fathomteam/moab/get/3fa8ae9e8025e0c1349792dee1f39b82ad8c1a33.tar.gz # 5.6.0
    URL_MD5 640aced62afe329596ca5953df57e79c
    DOWNLOAD_NAME moab.tar.gz
    DOWNLOAD_DIR ${CMAKE_CURRENT_SOURCE_DIR}/downloads
    DOWNLOAD_EXTRACT_TIMESTAMP OFF
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/moab
    CMAKE_ARGS ${MOAB_OPTIONS}
    LOG_CONFIGURE TRUE
    LOG_BUILD TRUE
    LOG_INSTALL TRUE
    BUILD_BYPRODUCTS ${PROJECT_BINARY_DIR}/lib/libMOAB.a
  )
  add_dependencies(MOAB MOAB_proj)

  include_directories(${PROJECT_BINARY_DIR}/include)
  link_directories(${PROJECT_BINARY_DIR}/lib)
endif()

# Add subdirectories
add_subdirectory(common)
add_subdirectory(emulatoratm)
add_subdirectory(driver)
