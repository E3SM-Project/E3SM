include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

if (DEFINED COMPILER AND NOT COMPILER MATCHES ".*gpu.*")
  include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-v100.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
  set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)
endif()

include (${EKAT_MACH_FILES_PATH}/mpi/other.cmake) # Unset all EKAT_MPI* params. Must specify them below.

set(EKAT_MPIRUN_EXE "jsrun -E LD_PRELOAD=/opt/ibm/spectrum_mpi/lib/pami_490/libpami.so" CACHE STRING "" FORCE)
set(EKAT_MPI_NP_FLAG "-n" CACHE STRING "" FORCE)
