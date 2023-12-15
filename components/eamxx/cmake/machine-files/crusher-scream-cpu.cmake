include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

