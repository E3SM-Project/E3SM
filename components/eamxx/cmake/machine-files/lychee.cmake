include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(SCREAM_INPUT_ROOT "/home/projects/e3sm/eamxx/data" CACHE STRING "")

set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE)

# Things take forever to build without this
set(SCREAM_SMALL_KERNELS On CACHE BOOL "" FORCE)

# Load H100/HOPPER arch and cuda backend for kokkos
include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-h100.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
