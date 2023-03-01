include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# Load all kokkos settings from Ekat's mach file
include (${EKAT_MACH_FILES_PATH}/kokkos/intel-skx.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

set (NetCDF_PATH /share/apps/netcdf/4.6.3/gcc/8.1.0 CACHE STRING "")

#Compy SLURM specific settings
set(EKAT_MPI_NP_FLAG "-p short -n" CACHE STRING "" FORCE)
