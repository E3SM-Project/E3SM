# CMake initial cache file for Linux 64bit RHEL7/CENTOS7
# tested with stock gcc/gfortran and packages from EPEL:
#    openmpi-devel 
#    blas-devel
#    lapack-devel
#    netcdf-fortran-devel
#
#
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

SET (WITH_PNETCDF FALSE CACHE FILEPATH "")

# hack until findnetcdf is updated to look for netcdf.mod
SET (ADD_Fortran_FLAGS "-I/usr/lib64/gfortran/modules" CACHE STRING "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

