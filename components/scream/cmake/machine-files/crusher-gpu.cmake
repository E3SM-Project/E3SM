# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

#serial is needed, but maybe it is always on?
#include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/mi250.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/hip.cmake)

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MACHINE "crusher-gpu" CACHE STRING "")

