# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

include (${EKAT_MACH_FILES_PATH}/kokkos/generic.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# Enable AMD Epyc arch in kokkos
set (Kokkos_ARCH_EPYC ON CACHE BOOL "")
