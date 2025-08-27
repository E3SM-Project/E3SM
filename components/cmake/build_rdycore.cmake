function(build_rdycore)

  if (COMP_NAMES MATCHES ".*rdycore.*")

    message(STATUS "Found rdycore component")

    include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)
    set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${EXTERNALS_SOURCE_DIR}/rdycore/cmake")

    ## FIXME: how much of this is needed??
    ## SCREAM manages its own kokkos settings. We can think about
    ## removing this once the cime cmake_macros aren't using obsolete
    ## Kokkos settings like KOKKOS_OPTIONS.
    #if (KOKKOS_OPTIONS)
    #  unset(KOKKOS_OPTIONS)
    #endif()
    #
    ## Strip leading space, or else add_compile_definitions will add
    ## an empty definition '-D `, which will cause compiler errors
    ## string (REGEX REPLACE "^ " "" CPPDEFS "${CPPDEFS}")
    ## add_compile_definitions expects a list of definitinos, each
    ## one provided without -D. So 1) replace ' ' with ';', and
    ## 2) strip '-D'/
    #string (REPLACE " " ";" CPPDEFS "${CPPDEFS}")
    #string (REPLACE "-D" "" CPPDEFS "${CPPDEFS}")
    #add_compile_definitions ("${CPPDEFS}")
    #
    ## Include machine file here. Prefer a compiler-specific file.
    #set(RDY_MACH_FILE_ROOT ${CMAKE_SOURCE_DIR}/rdycore/cmake/machine-files)
    #if (EXISTS ${RDY_MACH_FILE_ROOT}/${MACH}-${COMPILER}.cmake)
    #  include(${RDY_MACH_FILE_ROOT}/${MACH}-${COMPILER}.cmake)
    #else()
    #  include(${RDY_MACH_FILE_ROOT}/${MACH}.cmake)
    #endif()
    #
    ## The machine files may enable kokkos stuff we don't want
    #if (NOT compile_threaded)
    #  set(Kokkos_ENABLE_OPENMP FALSE)
    #endif()

    # PROJECT_SOURCE_DIR is E3SM/components
    set(EXTERNALS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../externals")
    add_subdirectory(${EXTERNALS_SOURCE_DIR}/rdycore ${CMAKE_BINARY_DIR}/externals/rdycore)
  endif()

endfunction(build_rdycore)
