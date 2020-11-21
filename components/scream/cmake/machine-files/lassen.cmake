# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-v100.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
set(NetCDF_C_PATHS /usr/gdata/climdat/lassen_builds/libs/netcdf/netcdf-c/netcdf-c-install/ CACHE STRING "")
set(NetCDF_Fortran_PATHS /usr/gdata/climdat/lassen_builds/libs/netcdf/netcdf-f/netcdf-f-install/ CACHE STRING "")
set(CMAKE_C_COMPILER /usr/tce/packages/spectrum-mpi/spectrum-mpi-rolling-release-gcc-7.3.1/bin/mpicc CACHE STRING "")
set(BLAS_LIBRARIES /usr/tcetmp/packages/lapack/lapack-3.9.0-gcc-7.3.1/lib/libblas.so CACHE STRING "")
set(LAPACK_LIBRARIES /usr/tcetmp/packages/lapack/lapack-3.9.0-gcc-7.3.1/lib/liblapack.so CACHE STRING "")
