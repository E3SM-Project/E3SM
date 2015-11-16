# CMake initial cache file for Mac OSX 10.8
# ANL tukey viz cluster with intel 11
# 
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

SET (FORCE_Fortran_FLAGS "-traceback -fp-model fast -ftz -g -O2" CACHE STRING "")
#SET (FORCE_Fortran_FLAGS "-traceback -openmp -fp-model precise -ftz -g -O2" CACHE STRING "")
SET (ENABLE_OPENMP FALSE CACHE BOOL "")


SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (NETCDF_DIR /soft/libraries/unsupported/netcdf-4.1.3 CACHE FILEPATH "")
SET (HDF5_DIR /soft/libraries/hdf5-1.8.9 CACHE FILEPATH "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

