# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
set(SCREAM_MACHINE "mappy" CACHE STRING "")

include (${EKAT_MACH_FILES_PATH}/mappy.cmake)
