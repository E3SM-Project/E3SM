include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# Copied from summit.cmake, only change is to SCREAM_MACHINE
# and remove support for CPU builds.
#
# Load all kokkos settings from Ekat's mach file
include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-v100.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake) # Unset all EKAT_MPI* params. Must specify them below.

set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

set(EKAT_MPIRUN_EXE "jsrun" CACHE STRING "" FORCE)
set(EKAT_MPI_NP_FLAG "-n" CACHE STRING "" FORCE)
