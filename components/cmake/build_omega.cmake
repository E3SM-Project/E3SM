function(build_omega)

  if (COMP_NAMES MATCHES ".*omega.*")

     # Load the machine's settings here, as build_model() and build_eamxx() do, so
     # Omega gets the machine's compiler flags. The scope stays within components/omega.
     include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

     # universal.cmake blanks the _RELEASE flags. CIME's CPPDEFS carry NDEBUG, but
     # reach Omega only when EAMxx is in the case, by build_eamxx()'s directory-scope
     # add_compile_definitions, so restore -DNDEBUG for C and C++.
     string(APPEND CMAKE_C_FLAGS_RELEASE   " -DNDEBUG")
     string(APPEND CMAKE_CXX_FLAGS_RELEASE " -DNDEBUG")

     # Set CIME source path relative to components
     set(CIMESRC_PATH "../cime/src")

     add_subdirectory("omega")

  endif()

endfunction(build_omega)
