cmake_minimum_required(VERSION 3.18)
project(emulators LANGUAGES C CXX Fortran)

include(FetchContent)
set(
    FETCHCONTENT_BASE_DIR
    "${CMAKE_CURRENT_SOURCE_DIR}/downloads"
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
set(EKAT_ENABLE_MPI OFF CACHE BOOL "" FORCE) # unexpected: ekat is not parsing mpi.h correctly?

FetchContent_Declare(
    ekat
    GIT_REPOSITORY https://github.com/E3SM-Project/EKAT.git
    GIT_TAG master
    GIT_SHALLOW TRUE
)

FetchContent_MakeAvailable(ekat)

# ------------------------------------------------------------------------------
# LibTorch
# ------------------------------------------------------------------------------

option(EMULATOR_ENABLE_LIBTORCH "Build the LibTorch inference backend" OFF)
set(EMULATOR_LIBTORCH_VERSION "2.5.1" CACHE STRING "LibTorch version to download")
set(EMULATOR_LIBTORCH_VARIANT "auto" CACHE STRING
    "LibTorch build to download: auto, cpu, cu118, cu121, cu124, rocm6.2")
set(EMULATOR_LIBTORCH_URL "" CACHE STRING
    "LibTorch archive to download if no install is found (overrides the variant)")

if(EMULATOR_ENABLE_LIBTORCH)
  # Use an existing install (CMAKE_PREFIX_PATH / Torch_DIR) if there is one,
  # otherwise fetch the prebuilt CPU distribution into the download cache.
  # Prefer the torch the Python interpreter uses: two LibTorch versions in one
  # process (e.g. with the embedded Python backend) fail with symbol errors.
  set(EMULATOR_TORCH_PYTHON "python3" CACHE STRING
    "Python interpreter with torch, used to locate LibTorch and write the test fixtures")
  if(NOT Torch_DIR AND NOT CMAKE_PREFIX_PATH)
    execute_process(COMMAND ${EMULATOR_TORCH_PYTHON} -c
        "import torch.utils; print(torch.utils.cmake_prefix_path)"
      RESULT_VARIABLE _torch_py_failed
      OUTPUT_VARIABLE _torch_py_prefix OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
    if(NOT _torch_py_failed)
      message(STATUS "Using LibTorch from ${EMULATOR_TORCH_PYTHON}: ${_torch_py_prefix}")
      list(APPEND CMAKE_PREFIX_PATH "${_torch_py_prefix}")
    endif()
  endif()
  find_package(Torch QUIET)
  if(NOT Torch_FOUND)
    set(_torch_url "${EMULATOR_LIBTORCH_URL}")
    if(NOT _torch_url)
      set(_torch_variant "${EMULATOR_LIBTORCH_VARIANT}")
      if(_torch_variant STREQUAL "auto")
        find_program(_torch_nvcc nvcc)
        if(_torch_nvcc OR DEFINED ENV{CUDA_HOME} OR DEFINED ENV{CUDA_PATH})
          set(_torch_variant cu124)
        else()
          set(_torch_variant cpu)
        endif()
      endif()
      if(_torch_variant STREQUAL "cpu")
        set(_torch_tag "cpu")
      else()
        set(_torch_tag "${_torch_variant}")
      endif()
      set(_torch_url "https://download.pytorch.org/libtorch/${_torch_variant}/libtorch-cxx11-abi-shared-with-deps-${EMULATOR_LIBTORCH_VERSION}%2B${_torch_tag}.zip")
    endif()
    message(STATUS "LibTorch not found; fetching ${_torch_url}")
    FetchContent_Declare(
        libtorch
        URL ${_torch_url}
        DOWNLOAD_EXTRACT_TIMESTAMP OFF
    )
    FetchContent_GetProperties(libtorch)
    if(NOT libtorch_POPULATED)
      FetchContent_Populate(libtorch)
    endif()
    list(APPEND CMAKE_PREFIX_PATH "${libtorch_SOURCE_DIR}")
    find_package(Torch REQUIRED)
  endif()
endif()

if(BUILD_EMULATOR_TESTS)
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
# add_subdirectory(emulatoratm)
# add_subdirectory(driver)
