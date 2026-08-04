cmake_minimum_required(VERSION 3.18)
project(emulators LANGUAGES CXX Fortran)

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
  # Enable testing at the top level
  enable_testing()
endif()

if (STANDALONE_MODE)
  # Fetch MOAB so we don't have to have it preinstalled
  add_library(MOAB STATIC IMPORTED GLOBAL)
  set_target_properties(MOAB PROPERTIES IMPORTED_LOCATION ${PROJECT_BINARY_DIR}/lib/libMOAB.a)
  include(ExternalProject)
  set(MOAB_OPTIONS -G ${CMAKE_GENERATOR})
  set(MOAB_OPTIONS ${MOAB_OPTIONS} -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR})
  set(MOAB_OPTIONS ${MOAB_OPTIONS} -DBUILD_SHARED_LIBS=OFF)
  set(MOAB_OPTIONS ${MOAB_OPTIONS} -DENABLE_FORTRAN=OFF)
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
