# For now, scream will rely on it's own build of kokkos rather than the
# one in sharedlib.

function(build_eamxx)

  if (COMP_NAMES MATCHES ".*scream.*")

    message(STATUS "Found scream component")

    include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

    # Transfer CIME-set flags into Cmake-style
    set (CMAKE_C_FLAGS ${CFLAGS})
    set (CMAKE_CXX_FLAGS ${CXXFLAGS})
    set (CMAKE_Fortran_FLAGS ${FFLAGS})

    # SCREAM manages its own kokkos settings. We can think about
    # removing this once the cime cmake_macros aren't using obsolete
    # Kokkos settings like KOKKOS_OPTIONS.
    if (KOKKOS_OPTIONS)
      unset(KOKKOS_OPTIONS)
    endif()

    # Strip leading space, or else add_compile_definitions will add
    # an empty definition '-D `, which will cause compiler errors
    # string (REGEX REPLACE "^ " "" CPPDEFS "${CPPDEFS}")
    # add_compile_definitions expects a list of definitinos, each
    # one provided without -D. So 1) replace ' ' with ';', and
    # 2) strip '-D'/
    string (REPLACE " " ";" CPPDEFS "${CPPDEFS}")
    string (REPLACE "-D" "" CPPDEFS "${CPPDEFS}")
    add_compile_definitions ("${CPPDEFS}")

    # Include machine file here
    include(${CMAKE_SOURCE_DIR}/eamxx/cmake/machine-files/${MACH}.cmake)
    add_subdirectory("eamxx")
  endif()

endfunction(build_eamxx)
