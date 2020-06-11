# For now, scream will rely on it's own build of kokkos rather than the
# one in sharedlib.

function(build_scream)

  if (COMP_NAMES MATCHES ".*scream.*")
    # Include machine file here
    message("Found scream component")
    include(${CMAKE_SOURCE_DIR}/scream/cmake/machine-files/${MACH}.cmake)

    # Grab CXX compiler from CIME
    include(${CASEROOT}/Macros.cmake)

    # Set build type since scream manages its own flags
    if (DEBUG)
      set(CMAKE_BUILD_TYPE Debug)
    else()
      set(CMAKE_BUILD_TYPE Release)
    endif()

    add_subdirectory("scream")
  endif()

endfunction(build_scream)
