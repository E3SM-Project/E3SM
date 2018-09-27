# 
# CMake initial cache file for Sandia's redsky
# configurd for redsky's netcdf-intel/4.1 module 
#
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

# Openmpi 1.8 only
SET (USE_MPI_OPTIONS "--map-by node:SPAN" CACHE FILEPATH "")

# Openmpi 1.6
#SET (USE_MPI_OPTIONS "-loadbalance" CACHE FILEPATH "")

# this is ignored if we use FORCE_Fortran_FLAGS
SET (ADD_Fortran_FLAGS "-traceback" CACHE STRING "")

# redsky upgrade 8/2017, need to load sems-netcdf module:
SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (NETCDF_DIR $ENV{SEMS_NETCDF_ROOT} CACHE FILEPATH "")


SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET (CPRNC_DIR /projects/ccsm/cprnc/build.toss3 CACHE FILEPATH "")
