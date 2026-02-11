include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

include (${EKAT_MACH_FILES_PATH}/kokkos/amd-zen3.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)

include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake)
set(EKAT_MPIRUN_EXE "mpiexec" CACHE STRING "" FORCE)
set(EKAT_MPI_NP_FLAG "-n" CACHE STRING "" FORCE)
set(EKAT_MPI_EXTRA_ARGS "--tag-output --bind-to core" CACHE STRING "")
