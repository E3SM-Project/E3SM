# CMake initial cache file for Mac OSX 10.8
# ANL mira/cestus machines 
# 
# note: you need to also set this before running CMAKE:
#    setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/soft/libraries/alcf/current/xl/BLAS/lib:/soft/libraries/alcf/current/xl/LAPACK/lib:
#
#

SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

SET (CMAKE_Fortran_FLAGS "-WF,-C! -O2 -Ipreqx_modules" CACHE STRING "")
#SET (FORCE_Fortran_FLAGS "-WF,-C!" CACHE STRING "")
SET (ENABLE_OPENMP FALSE CACHE BOOL "")


#SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (PNETCDF_DIR /soft/libraries/pnetcdf/1.3.1/cnk-xl/current CACHE FILEPATH "")
SET (NETCDF_DIR /soft/libraries/netcdf/4.2.1.1/cnk-xl/V1R2M0-20130417 CACHE FILEPATH "")
SET (HDF5_DIR /soft/libraries/hdf5/1.8.10/cnk-xl/current CACHE FILEPATH "")
SET (ZLIB_DIR /soft/libraries/alcf/current/xl/ZLIB CACHE FILEPATH "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

# add these to LD_LIBRARY_PATH so cmake will find them:
# /soft/libraries/alcf/current/xl/LAPACK/lib 
# /soft/libraries/alcf/current/xl/BLAS/lib 
