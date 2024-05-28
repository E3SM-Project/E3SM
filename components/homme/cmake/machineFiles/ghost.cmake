
# Compilers
SET (CMAKE_Fortran_COMPILER mpif90 CACHE FILEPATH "")
SET (CMAKE_C_COMPILER mpicc CACHE FILEPATH "")
SET (CMAKE_CXX_COMPILER mpicc CACHE FILEPATH "")

# Flag Settings
set(CMAKE_EXE_LINKER_FLAGS "-ldl" CACHE STRING "")
set(OPT_FLAGS "-O3" CACHE STRING "")
set(ADD_Fortran_FLAGS " -O3" CACHE STRING "")

SET (WITH_PNETCDF TRUE CACHE FILEPATH "")
SET (PnetCDF_PATH /projects/ccsm/src/parallel-netcdf-1.5.0-src/ CACHE FILEPATH "")

# Directory Pointers
#SET (NetCDF_Fortran_DIR /projects/ccsm/src/netcdf-fortran-4.4.0/ CACHE FILEPATH "")
#SET (NetCDF_Fortran_PATH /projects/ccsm/netcdf4.3.2-intel12/ CACHE FILEPATH "")
SET (NetCDF_Fortran_PATH /projects/ccsm/netcdf-legacy/netcdf-fortran-4.4.2/ CACHE FILEPATH "")
#SET (NetCDF_Fortran_INCLUDE_DIR /projects/ccsm/netcdf4.3.2-intel12/include/ CACHE FILEPATH "")
#SET (NetCDF_Fortran_LIBRARY /projects/ccsm/netcdf4.3.2-intel12/lib/libnetcdff.a CACHE FILEPATH "")

#SET (PnetCDF_PATH $ENV{SEMS_NETCDF_ROOT} CACHE FILEPATH "")
#SET (PNETCDF_DIR $ENV{SEMS_NETCDF_ROOT} CACHE FILEPATH "")
SET (HDF5_DIR $ENV{SEMS_HDF5_ROOT} CACHE FILEPATH "")
SET (ZLIB_DIR $ENV{SEMS_ZLIB_ROOT} CACHE FILEPATH "")
#SET (CPRNC_DIR /projects/ccsm/cprnc/build.toss3 CACHE FILEPATH "")

# Turn IO on/off
# SET (BUILD_HOMME_WITHOUT_PIOLIBRARY TRUE CACHE BOOL "")

SET (USE_QUEUING FALSE CACHE BOOL "")
SET (HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")

# Turn on Kokkos build for Theta and Preqx versions 
SET (Kokkos_ARCH_BDW ON CACHE BOOL "")
SET (BUILD_HOMME_THETA_KOKKOS TRUE CACHE BOOL "")
SET (BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
