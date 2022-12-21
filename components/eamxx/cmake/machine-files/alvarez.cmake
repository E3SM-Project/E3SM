# Load all kokkos settings from Ekat's mach file
set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

#message(STATUS "alvarez PROJECT_NAME=${PROJECT_NAME} USE_CUDA=${USE_CUDA} KOKKOS_ENABLE_CUDA=${KOKKOS_ENABLE_CUDA}")

include (${EKAT_MACH_FILES_PATH}/kokkos/amd-zen3.cmake)
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if (SMP_PRESENT)
    include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
  else()
    include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
  endif()
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
endif()

set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

#message(STATUS "alvarez CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()

set(SCREAM_MPIRUN_EXE "srun" CACHE STRING "")
set(SCREAM_MPI_NP_FLAG "-n" CACHE STRING "")
set(SCREAM_MACHINE "alvarez" CACHE STRING "")
