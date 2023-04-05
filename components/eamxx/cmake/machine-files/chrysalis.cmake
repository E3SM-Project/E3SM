set (SCREAM_INPUT_ROOT "/lcrc/group/e3sm/ccsm-data/inputdata/" CACHE STRING "" FORCE)

set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MACHINE "chrysalis" CACHE STRING "")



