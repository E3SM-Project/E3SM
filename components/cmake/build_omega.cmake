function(build_omega)

  if (COMP_NAMES MATCHES ".*omega.*")

     message(STATUS "Found OMEGA component")

     include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

     # Set CIME source path relative to components
     set(CIMESRC_PATH "../cime/src")

     add_subdirectory("omega")

  endif()

endfunction(build_omega)
