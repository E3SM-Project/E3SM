# CMake initial cache file
#
# This machine file works with either Intel or gnu
# (selected by which modules are loaded)
#
#
# Perlmutter generic MPI enabled compiler wrappers:
SET (CMAKE_Fortran_COMPILER ftn   CACHE FILEPATH "")
SET (CMAKE_C_COMPILER       cc    CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER     CC    CACHE FILEPATH "")

# For some reason, cmake does not yet know CMAKE_Fortran_COMPILER_ID here, so use a manual approach until can
# fix more elegantly
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

# 1. Detection Logic for the "Big Three" on NERSC
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

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")

EXECUTE_PROCESS(COMMAND nf-config --prefix
  RESULT_VARIABLE NFCONFIG_RESULT
  OUTPUT_VARIABLE NFCONFIG_OUTPUT
  ERROR_VARIABLE  NFCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
SET (NetCDF_Fortran_PATH "${NFCONFIG_OUTPUT}" CACHE STRING "")

EXECUTE_PROCESS(COMMAND nc-config --prefix
  RESULT_VARIABLE NCCONFIG_RESULT
  OUTPUT_VARIABLE NCCONFIG_OUTPUT
  ERROR_VARIABLE  NCCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
SET (NetCDF_C_PATH "${NCCONFIG_OUTPUT}" CACHE STRING "")

SET (USE_QUEUING FALSE CACHE BOOL "")
# for standalone HOMME builds:
SET(CPRNC_DIR /global/cfs/cdirs/e3sm/tools/cprnc CACHE FILEPATH "")

SET (HOMMEXX_BFB_TESTING      TRUE    CACHE BOOL "")
SET (BUILD_HOMME_PREQX_KOKKOS FALSE   CACHE BOOL "")
SET (BUILD_HOMME_THETA_KOKKOS FALSE   CACHE BOOL "")
SET (HOMME_TESTING_PROFILE    "short" CACHE STRING "")

SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

IF (MY_COMPILER_TYPE MATCHES "Intel")
    # Works for both 'Intel' (ifort) and 'IntelLLVM' (ifx)
    SET (ADD_Fortran_FLAGS "-traceback -fp-model strict -O1" CACHE STRING "")
    SET (ADD_C_FLAGS       "-fp-model strict -O1" CACHE STRING "")
    SET (ADD_CXX_FLAGS     "-fp-model strict -O1" CACHE STRING "")

    IF(DEFINED ENV{MKLROOT})
      SET (HOMME_USE_MKL TRUE CACHE BOOL "")
    ENDIF()

ELSEIF (MY_COMPILER_TYPE STREQUAL "GNU")
    # GNU uses -fbacktrace instead of -traceback
    SET(ADD_Fortran_FLAGS "-fbacktrace -g -O0" CACHE STRING "")
    SET(ADD_C_FLAGS       "-g -O0" CACHE STRING "")
    SET(ADD_CXX_FLAGS     "-g -O0" CACHE STRING "")

ELSEIF (MY_COMPILER_TYPE MATCHES "NVHPC")
    # NVIDIA/PGI compilers
    SET(ADD_Fortran_FLAGS "-traceback -O1" CACHE STRING "")
    SET(ADD_C_FLAGS       "-traceback -O1" CACHE STRING "")
    SET(ADD_CXX_FLAGS     "-traceback -O1" CACHE STRING "")

ELSE()
    MESSAGE(STATUS "Unknown compiler ID: ${MY_COMPILER_TYPE}. No specific flags set.")
ENDIF()

SET(USE_MPIEXEC "srun" CACHE STRING "")
SET(USE_MPI_OPTIONS "-K --cpu_bind=cores" CACHE STRING "")
