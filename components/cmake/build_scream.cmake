# For now, scream will rely on it's own build of kokkos rather than the
# one in sharedlib.

function(build_scream)

  if (COMP_NAMES MATCHES ".*scream.*")

    message(STATUS "Found scream component")

    include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

    set (CMAKE_C_FLAGS ${CFLAGS})
    set (CMAKE_CXX_FLAGS ${CXXFLAGS})
    set (CMAKE_Fortran_FLAGS ${FFLAGS})
    add_definitions (${CPPDEFS})

    # Include machine file here
    include(${CMAKE_SOURCE_DIR}/scream/cmake/machine-files/${MACH}.cmake)
    add_subdirectory("scream")
  endif()

endfunction(build_scream)
