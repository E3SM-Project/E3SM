# CMake initial cache file for Anvil
#
# 2017/11 MT:  Set NETCDF and PNETCDF paths using same approach as CIME
#              Rely on CIME's FindNetCDF and FindPnetCDF for all remaining config
#
# tested 2017/11 with the following in .soft:
#
# +intel-17.0.0
# +mvapich2-2.2-intel-17.0.0-acme
# +netcdf-c-4.4.1-f77-4.4.4-intel-17.0.0-serial
# +pnetcdf-1.7.0-intel-17.0.0-mvapich2-2.2-acme
#
#
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicxx CACHE FILEPATH "")


# anvil uses softenv, which does not set environment variables but sets
# up path for executables. Compute NETCDF_DIR, PNETCDF_DIR
EXECUTE_PROCESS(COMMAND nf-config --prefix
  RESULT_VARIABLE NFCONFIG_RESULT
  OUTPUT_VARIABLE NETCDF_PATH
  ERROR_VARIABLE  NFCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
IF (${NFCONFIG_ERROR})
  MESSAGE(WARNING "nf-config --prefix produced an error. Default linking will be used.")
ELSE()
  MESSAGE("-- NETCDF_PATH = ${NETCDF_PATH}")
  SET (NETCDF_DIR ${NETCDF_PATH} CACHE FILEPATH "")
ENDIF ()

EXECUTE_PROCESS(COMMAND which pnetcdf_version
  COMMAND xargs dirname
  COMMAND xargs dirname
  RESULT_VARIABLE NFCONFIG_RESULT
  OUTPUT_VARIABLE PNETCDF_PATH
  ERROR_VARIABLE  NFCONFIG_ERROR
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
IF (${NFCONFIG_ERROR})
  MESSAGE(WARNING "which pnetcdf_version produced an error. PNETCDF disabled.")
  SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
ELSE()
  MESSAGE("-- PNETCDF_PATH = ${PNETCDF_PATH}")
  SET (PNETCDF_DIR ${PNETCDF_PATH} CACHE FILEPATH "")
ENDIF ()



SET (HOMME_USE_MKL "TRUE" CACHE FILEPATH "") # for Intel
SET (USE_QUEUING FALSE CACHE BOOL "")
# for standalone HOMME builds:
SET (CPRNC_DIR /home/ccsm-data/tools CACHE FILEPATH "")

# suppress "-openmp deprecated" warning number 10411
SET (ADD_Fortran_FLAGS "-diag-disable 10411" CACHE STRING "")
SET (ADD_C_FLAGS "-diag-disable 10411" CACHE STRING "")
SET (ADD_CXX_FLAGS "-diag-disable 10411" CACHE STRING "")



