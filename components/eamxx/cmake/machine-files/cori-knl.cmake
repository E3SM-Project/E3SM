set (EKAT_MACH_FILES_PATH ${CMAKE_CURRENT_LIST_DIR}/../../../../externals/ekat/cmake/machine-files)

# Load knl arch and openmp backend for kokkos
include (${EKAT_MACH_FILES_PATH}/kokkos/intel-knl.cmake)

if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if (SMP_PRESENT)
    include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
  else()
    include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
  endif()
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
endif()

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
set(SCREAM_MACHINE "cori-knl" CACHE STRING "")
