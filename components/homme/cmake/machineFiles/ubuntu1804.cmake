# CMake initial cache file for Ubuntu 18.04
# tested with stock gcc/gfortran & MPICH
# assumes self-built netcdf/pnetcdf/cprnc
# Users should modify the paths below to point to these

SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")
SET (PNETCDF_DIR /home/celdred/parallel-netcdf-gnu CACHE FILEPATH "")
SET (NetCDF_C_PATH /home/celdred/netcdf-gnu/netcdf-c CACHE FILEPATH "")
SET (NetCDF_Fortran_PATH /home/celdred/netcdf-gnu/netcdf-f CACHE FILEPATH "")

SET (USE_MPIEXEC "mpiexec.mpich" CACHE STRING "")

SET (NetCDF_C_INCLUDE_DIR /home/celdred/netcdf-gnu/netcdf-c/include CACHE FILEPATH "")
SET (NetCDF_C_LIBRARY /home/celdred/netcdf-gnu/netcdf-c/lib/libnetcdf.so CACHE FILEPATH "")
SET (NetCDF_Fortran_INCLUDE_DIR /home/celdred/netcdf-gnu/netcdf-f/include/ CACHE FILEPATH "")
SET (NetCDF_Fortran_LIBRARY /home/celdred/netcdf-gnu/netcdf-f/lib/libnetcdff.a CACHE FILEPATH "")


SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

SET (CPRNC_DIR /home/celdred/E3SM/components/homme/swe_tests CACHE FILEPATH "")
