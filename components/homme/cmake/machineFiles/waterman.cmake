# Note: CMAKE_CXX_COMPILER needs to be set to the path of nvcc_wrapper
# nvcc_wrapper will choose either the Nvidia Cuda compiler or the OpenMP compiler depending on what's being compiled

SET(BUILD_HOMME_PREQX_KOKKOS TRUE CACHE BOOL "")
SET (WITH_PNETCDF FALSE CACHE FILEPATH "")

set (mynetcdf /ascldap/users/onguba/waterman-netcdf4)

SET(NETCDF_DIR ${mynetcdf} CACHE FILEPATH "")
set (ZLIB_DIR $ENV{ZLIB_ROOT} CACHE FILEPATH "")
set (CURL_ROOT $ENV{CURL_ROOT} CACHE FILEPATH "")
set (CURL_LIBRARY -L$ENV{CURL_ROOT}/lib -lcurl CACHE LIST "")

set (NetCDF_Fortran_DIR ${mynetcdf} CACHE FILEPATH "")
set (NetCDF_Fortran_INCLUDE_DIR ${mynetcdf}/include CACHE FILEPATH "")
set (NetCDF_Fortran_LIBRARY ${mynetcdf}/lib -lnetcdff -lhdf5_hl -lhdf5 CACHE LIST "")
set (HDF5_HL_INCLUDE_DIR ${mynetcdf}/include CACHE FILEPATH "")
set (HDF5_HL_LIBRARY ${mynetcdf}/lib -lnetcdff -lhdf5_hl -lhdf5 CACHE LIST "")
set (HDF5_C_INCLUDE_DIR ${mynetcdf}/include CACHE FILEPATH "")
set (HDF5_C_LIBRARY ${mynetcdf}/lib -lnetcdff -lhdf5_hl -lhdf5 CACHE LIST "")

SET (HAVE_EXTRAE TRUE CACHE BOOL "")
SET (Extrae_LIBRARY "-L${mynetcdf}/lib -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -ldl -lz" CACHE STRING "")


SET(HOMME_FIND_BLASLAPACK TRUE CACHE BOOL "")
SET(USE_QUEUING FALSE CACHE BOOL "")
SET(ENABLE_CUDA TRUE CACHE BOOL "")
SET(USE_TRILINOS OFF CACHE BOOL "")

SET(KOKKOS_PATH "/home/onguba/kokkos/build-serial-cuda-nodebug/" CACHE STRING "")

SET(CMAKE_C_COMPILER "mpicc" CACHE STRING "")
SET(CMAKE_CXX_COMPILER "/home/onguba/kokkos/bin/nvcc_wrapper" CACHE STRING "")
SET(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "")

set (ENABLE_COLUMN_OPENMP CACHE BOOL OFF)
set (ENABLE_HORIZ_OPENMP CACHE BOOL OFF)

set (USE_NUM_PROCS 8 CACHE STRING "")




