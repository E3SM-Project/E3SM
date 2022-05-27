# For now, scream will rely on it's own build of kokkos rather than the
# one in sharedlib.

function(build_scream)

  if (COMP_NAMES MATCHES ".*scream.*")

    # Include machine file here
    message("Found scream component, USE_CUDA ${USE_CUDA}")
    if (USE_CUDA)
      include(${CMAKE_SOURCE_DIR}/scream/cmake/machine-files/${MACH}.cmake)
    else()
      include(${CMAKE_SOURCE_DIR}/scream/cmake/machine-files/${MACH}-cpu.cmake)
    endif()
    add_subdirectory("scream")
  endif()

endfunction(build_scream)
