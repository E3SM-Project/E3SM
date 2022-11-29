# Copied from summit.cmake, only change is to SCREAM_MACHINE
# and remove support for CPU builds.
#
# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-v100.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)

set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

set(SCREAM_MPIRUN_EXE "jsrun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-n" CACHE STRING "")
set(SCREAM_MACHINE "ascent" CACHE STRING "")
