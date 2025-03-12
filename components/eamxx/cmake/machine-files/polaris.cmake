include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

if (USE_CUDA)
  include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-a100.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/amd-zen3.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
endif()

include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake)

set(EKAT_MPIRUN_EXE "mpiexec" CACHE STRING "" FORCE)
set(EKAT_MPI_NP_FLAG "-np" CACHE STRING "" FORCE)
set(EKAT_MPI_EXTRA_ARGS "--label --cpu-bind core -d 8" CACHE STRING "")
