include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# Remove this if you are using a resource manager (slurm etc)
set (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES True CACHE BOOL "")

# EKAT MPI settings
set (EKAT_MPIRUN_EXE "mpiexec" CACHE STRING "mpiexec")
set (EKAT_MPI_NP_FLAG "-n" CACHE STRING "-n")

include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake)
