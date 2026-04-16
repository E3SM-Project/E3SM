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

# Compiler wrappers (Cray PE)
set(SCC "cc" CACHE STRING "" FORCE)
set(SCXX "CC" CACHE STRING "" FORCE)
set(SFC "ftn" CACHE STRING "" FORCE)

# Common flags (from gnugpu.cmake)
set(CMAKE_C_FLAGS " -mcmodel=medium" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS " -mcmodel=medium -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none -fallow-argument-mismatch" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_RELEASE "-O2" CACHE STRING "" FORCE)
set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g" CACHE STRING "" FORCE)

if (USE_CUDA)
  # From polaris_gnugpu.cmake
  set(CMAKE_CUDA_FLAGS "-ccbin CC -O2 -arch sm_80 --use_fast_math" CACHE STRING "" FORCE)
  set(CMAKE_EXE_LINKER_FLAGS " -L$ENV{CUDA_HOME}/lib64/ -lstdc++" CACHE STRING "" FORCE)
endif()
