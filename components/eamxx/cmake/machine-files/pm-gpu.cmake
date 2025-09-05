include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

set(EKAT_MPI_EXTRA_ARGS "${EKAT_MPI_EXTRA_ARGS} --gpus-per-task=1" CACHE STRING "" FORCE)

#option(Kokkos_ARCH_AMPERE80 "" ON)
set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

option(SCREAM_SMALL_KERNELS "Use small, non-monolothic kokkos kernels for ALL components that support them" ON)

#message(STATUS "pm-gpu CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()
