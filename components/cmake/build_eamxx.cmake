# This is basically just an add_subdirectory call now that eamxx
# uses E3SM/CIME settings.

function(build_eamxx)

  if (COMP_NAMES MATCHES ".*scream.*")

    message(STATUS "Found scream component")

    include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

    add_subdirectory("eamxx")
  endif()

endfunction(build_eamxx)
