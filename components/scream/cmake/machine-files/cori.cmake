set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

# Load knl arch and openmp backend for kokkos
include (${EKAT_MACH_FILES_PATH}/kokkos/intel-knl.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# Fixes some openmpi link problems we observed on cori. This hack is
# not necessary if CRAYPE_LINK_TYPE=dynamic is in the environment.
set(SCREAM_CORI_HACK True CACHE BOOL "")
