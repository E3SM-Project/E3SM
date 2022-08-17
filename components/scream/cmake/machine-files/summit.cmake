# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

if (DEFINED COMPILER AND NOT COMPILER MATCHES ".*gpu.*")
  include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-v100.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
  set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)
endif()

set(SCREAM_MPIRUN_EXE "jsrun -E LD_PRELOAD=/opt/ibm/spectrum_mpi/lib/pami_490/libpami.so" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-n" CACHE STRING "")
set(SCREAM_MACHINE "summit" CACHE STRING "")
