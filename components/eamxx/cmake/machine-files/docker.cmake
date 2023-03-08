include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

#message(STATUS "gcp PROJECT_NAME=${PROJECT_NAME} USE_CUDA=${USE_CUDA} KOKKOS_ENABLE_CUDA=${KOKKOS_ENABLE_CUDA}")
# use default backend?

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

if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()

set(BLAS_LIBRARIES /opt/conda/lib/libblas.so CACHE STRING "")
set(LAPACK_LIBRARIES /opt/conda/lib/liblapack.so CACHE STRING "")
set(SCREAM_INPUT_ROOT "/storage/inputdata/" CACHE STRING "")