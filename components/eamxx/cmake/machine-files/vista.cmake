include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

#message(STATUS "vista PROJECT_NAME=${PROJECT_NAME} USE_CUDA=${USE_CUDA} KOKKOS_ENABLE_CUDA=${KOKKOS_ENABLE_CUDA}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if (USE_CUDA)
    include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-h100.cmake) # H100=Hopper=Ampere90=Hopper90 Kokkos_ARCH_HOPPER90
    include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
  else()
    include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-grace.cmake) # KOKKOS_ARCH_ARMV9_GRACE
    include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
    #include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
  endif()
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/nvidia-h100.cmake)
  include (${EKAT_MACH_FILES_PATH}/kokkos/cuda.cmake)
endif()

include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake) # should be changed to use ibrun

set(EKAT_MPI_EXTRA_ARGS "${EKAT_MPI_EXTRA_ARGS} --gpus-per-task=1" CACHE STRING "" FORCE)

#option(Kokkos_ARCH_AMPERE90 "" ON)
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

#set(CMAKE_CUDA_FLAGS "-allow-unsupported-compiler" CACHE STRING "" FORCE) #ndk try for gcc14

#message(STATUS "vista CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()
