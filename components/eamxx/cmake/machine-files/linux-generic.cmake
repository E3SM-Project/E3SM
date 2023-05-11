include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# Remove this if you are using a resource manager (slurm etc)
set (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES True CACHE BOOL "")
