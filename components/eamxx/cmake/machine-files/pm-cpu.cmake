include(${CMAKE_CURRENT_LIST_DIR}/common.cmake)
common_setup()

if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if (SMP_PRESENT)
    include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
    #message(STATUS, "pm-cpu openmp SMP_PRESENT=${SMP_PRESENT}")
  else()
    include (${EKAT_MACH_FILES_PATH}/kokkos/serial.cmake)
    #message(STATUS, "pm-cpu serial SMP_PRESENT=${SMP_PRESENT}")
  endif()
else()
  include (${EKAT_MACH_FILES_PATH}/kokkos/openmp.cmake)
endif()

set(CMAKE_CXX_FLAGS "-DTHRUST_IGNORE_CUB_VERSION_CHECK" CACHE STRING "" FORCE)

#message(STATUS "pm-cpu CMAKE_CXX_COMPILER_ID=${CMAKE_CXX_COMPILER_ID} CMAKE_Fortran_COMPILER_VERSION=${CMAKE_Fortran_COMPILER_VERSION}")
if ("${PROJECT_NAME}" STREQUAL "E3SM")
  if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
      set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
    endif()
  endif()
else()
  set(CMAKE_Fortran_FLAGS "-fallow-argument-mismatch"  CACHE STRING "" FORCE) # only works with gnu v10 and above
endif()
