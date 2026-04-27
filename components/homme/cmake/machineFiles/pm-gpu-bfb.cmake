# CMake initial cache file for pm-gpu for BFB testing
#
# This machine file works with either Intel or gnu
# (selected by which modules are loaded)
#
# Perlmutter generic MPI enabled compiler wrappers:
SET (CMAKE_Fortran_COMPILER ftn   CACHE FILEPATH "")
SET (CMAKE_C_COMPILER       cc    CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER     CC    CACHE FILEPATH "")

# For some reason, cmake does not yet know CMAKE_Fortran_COMPILER_ID here,
# so use a manual approach until can fix more elegantly
#project(TestProject LANGUAGES Fortran)
#enable_language(Fortran)
#MESSAGE(STATUS "Fortran Compiler: ${CMAKE_Fortran_COMPILER}")
#MESSAGE(STATUS "Fortran Compiler ID: ${CMAKE_Fortran_COMPILER_ID}")

# Manually probe the compiler before project()
execute_process(
    COMMAND ftn --version
    OUTPUT_VARIABLE FTN_VERSION_RAW
    ERROR_VARIABLE  FTN_VERSION_RAW
    RESULT_VARIABLE FTN_PROBE_RESULT
)

# Extract the first line for parsing
string(REGEX MATCH "^[^\n]*" FTN_VERSION_SUMMARY "${FTN_VERSION_RAW}")

# 1. Detection Logic
if(FTN_VERSION_SUMMARY MATCHES "ifx|ifort|Intel")
    set(MY_COMPILER_TYPE "Intel" CACHE STRING "Manual Compiler Identification")
    message(STATUS "Detected Intel Compiler (oneAPI/Classic)")

elseif(FTN_VERSION_SUMMARY MATCHES "gfortran|GNU")
    set(MY_COMPILER_TYPE "GNU" CACHE STRING "Manual Compiler Identification")
    message(STATUS "Detected GNU Compiler (gfortran)")

elseif(FTN_VERSION_SUMMARY MATCHES "nvfortran|NVIDIA|pgf90|PGI")
    set(MY_COMPILER_TYPE "NVHPC" CACHE STRING "Manual Compiler Identification")
    message(STATUS "Detected NVIDIA/NVHPC Compiler")
else()
    set(MY_COMPILER_TYPE "Unknown" CACHE STRING "Manual Compiler Identification")
    message(WARNING "Could not identify compiler from: ${FTN_VERSION_SUMMARY}")
endif()

# Set kokkos arch, to get correct avx flags
SET (Kokkos_ARCH_ZEN3 ON CACHE BOOL "")

SET (WITH_PNETCDF FALSE CACHE BOOL "")

EXECUTE_PROCESS(COMMAND nf-config --prefix
  RESULT_VARIABLE NFCONFIG_RESULT
  OUTPUT_VARIABLE NFCONFIG_OUTPUT
  ERROR_VARIABLE  NFCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
IF(NOT NFCONFIG_RESULT EQUAL 0)
  MESSAGE(FATAL_ERROR "nf-config failed: ${NFCONFIG_ERROR}")
ENDIF()
SET (NetCDF_Fortran_PATH "${NFCONFIG_OUTPUT}" CACHE STRING "")

EXECUTE_PROCESS(COMMAND nc-config --prefix
  RESULT_VARIABLE NCCONFIG_RESULT
  OUTPUT_VARIABLE NCCONFIG_OUTPUT
  ERROR_VARIABLE  NCCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
IF(NOT NCCONFIG_RESULT EQUAL 0)
  MESSAGE(FATAL_ERROR "nc-config failed: ${NCCONFIG_ERROR}")
ENDIF()
SET (NetCDF_C_PATH "${NCCONFIG_OUTPUT}" CACHE STRING "")

SET (USE_QUEUING FALSE CACHE BOOL "")
# for standalone HOMME builds:
SET(CPRNC_DIR /global/cfs/cdirs/e3sm/tools/cprnc CACHE FILEPATH "")

SET (HOMMEXX_BFB_TESTING      TRUE    CACHE BOOL "")
SET (BUILD_HOMME_PREQX_KOKKOS TRUE    CACHE BOOL "")
SET (BUILD_HOMME_THETA_KOKKOS TRUE    CACHE BOOL "")
SET (HOMME_TESTING_PROFILE    "short" CACHE STRING "")

SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

IF (MY_COMPILER_TYPE MATCHES "Intel")
    # Works for both 'Intel' (ifort) and 'IntelLLVM' (ifx)
    # Note Intel compilers not yet working for pm-gpu, they may one day be.
    # Enable MKL only when MKLROOT points to an Intel MKL installation.
    # MKLROOT may be set in GNU environments too, so we verify the path.
    IF(DEFINED ENV{MKLROOT})
      SET(_mkl_root "$ENV{MKLROOT}")
      IF(_mkl_root MATCHES "intel" OR _mkl_root MATCHES "oneapi/mkl")
        SET(HOMME_USE_MKL TRUE CACHE BOOL "")
      ELSE()
        MESSAGE(STATUS "MKLROOT set but does not appear to be an Intel MKL path (${_mkl_root}); not enabling HOMME_USE_MKL")
      ENDIF()
      UNSET(_mkl_root)
    ENDIF()

ELSEIF (MY_COMPILER_TYPE STREQUAL "GNU")
    SET(ADD_Fortran_FLAGS "-fbacktrace -ffp-contract=off -g -O0" CACHE STRING "")
    SET(ADD_CXX_FLAGS "-ffp-contract=off -g -O0" CACHE STRING "")
    SET(ADD_C_FLAGS "-ffp-contract=off -g -O0" CACHE STRING "")

ELSEIF (MY_COMPILER_TYPE MATCHES "NVHPC")
    SET(ADD_Fortran_FLAGS "-traceback -g -O0" CACHE STRING "")
    SET(ADD_C_FLAGS       "-g -O0" CACHE STRING "")
    SET(ADD_CXX_FLAGS     "-g -O0" CACHE STRING "")
ELSE()
    MESSAGE(STATUS "Unknown compiler ID: ${MY_COMPILER_TYPE}. No specific flags set.")
ENDIF()

#SET(HOMMEXX_EXEC_SPACE CUDA CACHE STRING "")
SET(HOMMEXX_EXEC_SPACE SERIAL CACHE STRING "")
#SET(HOMMEXX_MPI_ON_DEVICE FALSE CACHE BOOL "")
#SET(HOMMEXX_CUDA_MAX_WARP_PER_TEAM "16" CACHE STRING  "")

SET(BUILD_HOMME_WITHOUT_PIOLIBRARY FALSE CACHE BOOL "")

SET(USE_TRILINOS OFF CACHE BOOL "")

SET(Kokkos_ENABLE_OPENMP OFF CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA OFF CACHE BOOL "")
SET(Kokkos_ENABLE_CUDA_LAMBDA OFF CACHE BOOL "")
SET(Kokkos_ARCH_AMPERE80 ON CACHE BOOL "")
SET(Kokkos_ENABLE_IMPL_CUDA_MALLOC_ASYNC OFF CACHE BOOL "")
#SET(Kokkos_ARCH_ZEN3 ON CACHE BOOL "") # works, and perf same if both AMPERE80 and ZEN3 are on
#SET(Kokkos_ENABLE_CUDA_UVM ON CACHE BOOL "")
SET(Kokkos_ENABLE_EXPLICIT_INSTANTIATION OFF CACHE BOOL "")
#SET(Kokkos_ENABLE_CUDA_ARCH_LINKING OFF CACHE BOOL "")

SET(CXXLIB_SUPPORTED_CACHE FALSE CACHE BOOL "")

SET(ENABLE_OPENMP OFF CACHE BOOL "")
SET(ENABLE_COLUMN_OPENMP OFF CACHE BOOL "")
SET(ENABLE_HORIZ_OPENMP OFF CACHE BOOL "")

SET(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "")

SET(USE_NUM_PROCS 4 CACHE STRING "") # only default

SET(USE_MPIEXEC "srun" CACHE STRING "")
SET(USE_MPI_OPTIONS "-K --cpu-bind=cores" CACHE STRING "")
SET(HOMME_TESTING_TIMELIMIT 1800 CACHE STRING "")
