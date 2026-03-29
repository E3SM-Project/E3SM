include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

# No resource manager on GitHub-hosted runners; manage resources internally
set (EKAT_TEST_LAUNCHER_MANAGE_RESOURCES True CACHE BOOL "")

set (EKAT_MPIRUN_EXE "mpirun" CACHE STRING "")
set (EKAT_MPI_NP_FLAG "-n" CACHE STRING "")

# Suppress GCC 10+ errors for legacy Fortran argument mismatches in external dependencies
set (CMAKE_Fortran_FLAGS "-fallow-argument-mismatch" CACHE STRING "" FORCE)
