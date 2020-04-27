# CMake initial cache file for Linux 64bit RHEL6/CENTOS6
# tested with stock gcc/gfortran & openmpi 
#
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

SET (WITH_PNETCDF TRUE CACHE FILEPATH "")
SET (PNETCDF_DIR /home/celdred/parallel-netcdf-gnu CACHE FILEPATH "")
SET (NetCDF_C_INCLUDE_DIR /usr/include/ CACHE FILEPATH "")
SET (NetCDF_Fortran_INCLUDE_DIR /usr/include/ CACHE FILEPATH "")
SET (NetCDF_C_LIBRARY /usr/lib/x86_64-linux-gnu/libnetcdf.so CACHE FILEPATH "")
SET (NetCDF_Fortran_LIBRARY /usr/lib/x86_64-linux-gnu/libnetcdff.a CACHE FILEPATH "")


# hack until findnetcdf is updated to look for netcdf.mod
# but this is ignored by cprnc
#SET (ADD_Fortran_FLAGS "-I/usr/lib64/gfortran/modules" CACHE STRING "")


SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

