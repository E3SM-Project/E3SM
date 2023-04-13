include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/generic.cmake)

# Additional settings
set(Kokkos_ENABLE_AGGRESSIVE_VECTORIZATION TRUE CACHE BOOL "")
set(Kokkos_ENABLE_OPENMP FALSE CACHE BOOL "")

# Remove this if you are using a resource manager (slurm etc)
set (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES True CACHE BOOL "")
