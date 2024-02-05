include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

# Load knl arch and openmp backend for kokkos
include (${EKAT_MACH_FILES_PATH}/kokkos/intel-knl.cmake)
include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)

include (${EKAT_MACH_FILES_PATH}/mpi/srun.cmake)

if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
       set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()

# Fixes some openmpi link problems we observed on cori. This hack is
# not necessary if CRAYPE_LINK_TYPE=dynamic is in the environment.
set(SCREAM_CORI_HACK True CACHE BOOL "")
