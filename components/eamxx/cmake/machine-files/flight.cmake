include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

# Get AMD arch settings
include(${EKAT_MACH_FILES_PATH}/kokkos/intel-skx.cmake)

# Add OpenMP settings in standalone mode OR e3sm with compile_threaded=ON
if (NOT "${PROJECT_NAME}" STREQUAL "E3SM" OR compile_threaded)
  include(${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
endif()

# Use srun for standalone testing
include(${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)
