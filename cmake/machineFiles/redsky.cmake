# 
# CMake initial cache file for Sandia's redsky
# configurd for redsky's netcdf-intel/4.1 module 
#
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

SET (USE_MPI_OPTIONS "-loadbalance" CACHE FILEPATH "")

# this is ignored if we use FORCE_Fortran_FLAGS
SET (ADD_Fortran_FLAGS "-traceback" CACHE STRING "")

# override cmake's intel defaults:
# default cmake options for Intel: 
#      -assume byterecl -fp-model precise -ftz -g -O3 
# -O3 causes problems on redsky when openMP is enabled (even for 1 thread)
#
#SET (FORCE_Fortran_FLAGS "-openmp -fp-model fast -ftz -g -O2" CACHE STRING "")
SET (FORCE_Fortran_FLAGS "-openmp -traceback -fp-model precise -ftz -g -O2" CACHE STRING "")


#SET (NETCDF_DIR "/projects/ccsm/netcdf4.2.1.1-intel" CACHE FILEPATH "")
SET (NETCDF_DIR "/opt/netcdf-4.1-intel" CACHE FILEPATH "")
SET (PNETCDF_DIR "/projects/ccsm/pnetcdf1.2.0-intel" CACHE FILEPATH "")
SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

