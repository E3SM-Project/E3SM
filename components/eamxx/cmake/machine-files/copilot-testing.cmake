include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# No resource manager in CI/container environments
set (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES True CACHE BOOL "")

# Fortran compatibility flags for newer gfortran
set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch" CACHE STRING "Fortran compiler flags" FORCE)

# Input data directory (set by setup-copilot-env.sh or agent)
if (DEFINED ENV{SCREAM_INPUT_ROOT})
  set(SCREAM_INPUT_ROOT "$ENV{SCREAM_INPUT_ROOT}" CACHE PATH "")
endif()

# BLAS/LAPACK - prefer env vars, fall back to system paths
if (DEFINED ENV{BLAS_ROOT})
  # Try lib first, then lib64 (common on HPC/module installs)
  set(_blas_lib_dir "lib")
  if (NOT EXISTS "$ENV{BLAS_ROOT}/${_blas_lib_dir}/libblas.so"
      AND EXISTS "$ENV{BLAS_ROOT}/lib64/libblas.so")
    set(_blas_lib_dir "lib64")
  endif()
  set(BLAS_LIBRARIES "$ENV{BLAS_ROOT}/${_blas_lib_dir}/libblas.so" CACHE STRING "")
  set(LAPACK_LIBRARIES "$ENV{BLAS_ROOT}/${_blas_lib_dir}/liblapack.so" CACHE STRING "")
elseif (EXISTS "/usr/lib/x86_64-linux-gnu/libblas.so")
  set(BLAS_LIBRARIES "/usr/lib/x86_64-linux-gnu/libblas.so" CACHE STRING "")
  set(LAPACK_LIBRARIES "/usr/lib/x86_64-linux-gnu/liblapack.so" CACHE STRING "")
endif()

# Help Scorpio find PnetCDF when pnetcdf-config is not available (e.g., apt install)
if (NOT DEFINED PnetCDF_C_PATH AND EXISTS "/usr/include/pnetcdf.h")
  set(PnetCDF_C_PATH "/usr" CACHE PATH "")
endif()

# MPI launcher settings
set(EKAT_MPIRUN_EXE "mpirun" CACHE STRING "")
set(EKAT_MPI_NP_FLAG "-n" CACHE STRING "")

# Allow running MPI as root (common in containers/CI)
set(EKAT_MPI_EXTRA_ARGS "--allow-run-as-root --oversubscribe" CACHE STRING "Extra args for mpirun")

# Disable deprecated Kokkos code warnings
option(Kokkos_ENABLE_DEPRECATED_CODE_4 "" OFF)
