# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/blake.cmake)

set(SCREAM_INPUT_ROOT "/home/projects/e3sm/scream/data" CACHE STRING "")
