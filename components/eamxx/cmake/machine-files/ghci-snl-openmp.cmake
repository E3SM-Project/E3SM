# Common settings for our ghci images
include(${CMAKE_CURRENT_LIST_DIR}/ghci.cmake)

# Set SCREAM_MACHINE
set(SCREAM_MACHINE ghci-openmp CACHE STRING "")

# Set OpenMP backend
set(EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
