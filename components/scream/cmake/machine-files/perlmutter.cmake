# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)
IF (USE_CUDA OR KOKKOS_ENABLE_CUDA)
  include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-a100.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
else()  
  include (${EKAT_MACH_FILES_PATH}/kokkos/amd-zen2.cmake)
  #include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
endif()

#option(Kokkos_ARCH_AMPERE80 "" ON)
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    # only works with gnu v10 and above
    set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE)
  endif()
endif()


set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-n" CACHE STRING "")
set(SCREAM_MACHINE "perlmutter" CACHE STRING "")
