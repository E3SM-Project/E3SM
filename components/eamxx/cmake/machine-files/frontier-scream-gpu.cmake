set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

include (${EKAT_MACH_FILES_PATH}/kokkos/mi250.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/hip.cmake)

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MACHINE "frontier-scream-gpu" CACHE STRING "")

set(CMAKE_CXX_FLAGS "--offload-arch=gfx90a -fno-gpu-rdc -I$ENV{MPICH_DIR}/include" CACHE STRING "" FORCE)
