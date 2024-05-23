# CMake initial cache file for Anvil
#
# 2017/11 MT:  Set NETCDF and PNETCDF paths using same approach as CIME
#              Rely on CIME's FindNetCDF and FindPnetCDF for all remaining config
#
SET (CMAKE_Fortran_COMPILER mpif90   CACHE FILEPATH "")
SET (CMAKE_C_COMPILER       mpicc    CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER     mpicxx   CACHE FILEPATH "")

# Set kokkos arch, to get correct avx flags
SET (Kokkos_ARCH_SKX ON CACHE BOOL "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
#
# anvil module system doesn't set environment variables, but will put
# nc-config in our path.  anvil seperates C and Fortran libraries,
# so we need to set both paths instead of NETCDF_DIR
#
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

#turn on preqxx target and thus strict fpmodel for F vs CXX comparison
SET (BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET (HOMMEXX_BFB_TESTING TRUE CACHE BOOL "")
SET (HOMME_TESTING_PROFILE "short" CACHE STRING "")
SET (BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")

SET (USE_MPIEXEC "mpirun" CACHE STRING "")
SET (USE_MPI_OPTIONS "--bind-to cores" CACHE STRING "")
SET (USE_QUEUING FALSE CACHE BOOL "")
# for standalone HOMME builds:
# SET (CPRNC_DIR /lcrc/group/e3sm/soft/tools/cprnc CACHE FILEPATH "")
OPTION(HOMME_FIND_BLASLAPACK "Whether to use system blas/lapack" ON)
