function(build_rdycore)

  if (COMP_NAMES MATCHES ".*rdycore.*")

    message(STATUS "Found rdycore component")

    include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

    # Override some RDycore options to simplify the build process.
    set(ENABLE_DRIVER OFF CACHE BOOL "disable RDycore driver" FORCE)
    set(ENABLE_TESTS OFF CACHE BOOL "disable RDycore tests" FORCE)

    # PROJECT_SOURCE_DIR is E3SM/components
    set(EXTERNALS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../externals")
    add_subdirectory(${EXTERNALS_SOURCE_DIR}/rdycore ${CMAKE_BINARY_DIR}/externals/rdycore)
  endif()

endfunction(build_rdycore)
