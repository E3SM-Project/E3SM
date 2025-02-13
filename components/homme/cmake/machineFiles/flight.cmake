# CMake initial cache file copied from homme/cmake/machineFiles/chrysalis.cmake
#
#The netcdf commands work if the following commands are used to setup your environment:
# module purge
# module use /projects/sems/acme-boca-modulefiles/env-module
# module load acme-boca-env
# module load sems-archive-git
# module load sems-archive-cmake
# module load gnu/10.3.1
# module load sems-archive-intel/21.3.0
# module load sems-archive-openmpi/4.1.4
# module load acme-netcdf/4.7.4/acme
#
#
EXECUTE_PROCESS(COMMAND which mpiifort RESULT_VARIABLE MPIIFORT_RESULT)
EXECUTE_PROCESS(COMMAND which ifort    RESULT_VARIABLE IFORT_RESULT)
IF (${MPIIFORT_RESULT} EQUAL 0 AND ${IFORT_RESULT} EQUAL 0)
  MESSAGE(STATUS "Found Intel compiler and Intel MPI: building with mpiifort/mpiicc/mpiicpc")
  SET (CMAKE_Fortran_COMPILER mpiifort CACHE FILEPATH "")
  SET (CMAKE_C_COMPILER       mpiicc   CACHE FILEPATH "")
  SET (CMAKE_CXX_COMPILER     mpiicpc  CACHE FILEPATH "")
ELSE()
  MESSAGE(STATUS "Did not detect ifort or mpiifort: building with mpif90/mpicc/mpicxx")
  SET (CMAKE_Fortran_COMPILER mpif90   CACHE FILEPATH "")
  SET (CMAKE_C_COMPILER       mpicc    CACHE FILEPATH "")
  SET (CMAKE_CXX_COMPILER     mpicxx   CACHE FILEPATH "")
ENDIF()

SET (USE_MPIEXEC "srun" CACHE STRING "")
SET (USE_MPI_OPTIONS "-K --cpu_bind=cores" CACHE STRING "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")

SET (Kokkos_ARCH_NATIVE ON CACHE BOOL "")
SET(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")

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
SET (CPRNC_DIR /projects/ccsm/cprnc/build.toss3 CACHE FILEPATH "")

IF (${IFORT_RESULT} EQUAL 0)
  SET (HOMME_USE_MKL "TRUE" CACHE FILEPATH "")
  # turn on additional intel compiler flags
  SET (ADD_Fortran_FLAGS "-traceback" CACHE STRING "")
  SET (ADD_C_FLAGS       "-traceback" CACHE STRING "")
  SET (ADD_CXX_FLAGS     "-traceback" CACHE STRING "")
ELSE()
  SET (MKLROOT $ENV{MKLROOT} CACHE FILEPATH "")
  SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
ENDIF()