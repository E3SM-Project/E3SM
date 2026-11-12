function(build_omega)

  if (COMP_NAMES MATCHES ".*omega.*")

     # Set CIME source path relative to components
     set(CIMESRC_PATH "../cime/src")

     add_subdirectory("omega")

  endif()

endfunction(build_omega)
