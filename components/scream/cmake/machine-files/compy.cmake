# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/intel-skx.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
#set(NetCDF_DIR /share/apps/netcdf/4.6.3/gcc/8.1.0 CACHE STRING "")
#set(NetCDF_C_PATHS /share/apps/netcdf/4.6.3/gcc/8.1.0/ CACHE STRING "")
#set(NetCDF_Fortran_PATHS /share/apps/netcdf/4.6.3/gcc/8.1.0 CACHE STRING "")
set (NetCDF_PATH /share/apps/netcdf/4.6.3/gcc/8.1.0 CACHE STRING "")

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-p short -n" CACHE STRING "")