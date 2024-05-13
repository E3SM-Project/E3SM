include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

message(STATUS "muller-gpu PROJECT_NAME=${PROJECT_NAME} USE_CUDA=${USE_CUDA} KOKKOS_ENABLE_CUDA=${KOKKOS_ENABLE_CUDA}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if (USE_CUDA)
    include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-a100.cmake)
    include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
  else()
    include (${EKAT_MACH_FILES_PATH}/kokkos/amd-zen3.cmake)
    include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
    #include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
  endif()
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-a100.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
endif()

include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

#option(Kokkos_ARCH_AMPERE80 "" ON)
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

message(STATUS "muller-gpu CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()
