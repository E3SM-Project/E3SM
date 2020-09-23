# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# Enable Broadwell arch in Kokkos
option(Kokkos_ARCH_SKX "" ON)

set(CMAKE_CXX_FLAGS "-w" CACHE STRING "")
