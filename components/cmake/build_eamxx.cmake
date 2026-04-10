# For now, scream will rely on it's own build of kokkos rather than the
# one in sharedlib.

function(build_eamxx)

  if (COMP_NAMES MATCHES ".*scream.*")

    message(STATUS "Found scream component")

    include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

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

    # Include machine file here. Prefer a compiler-specific file.
    # Machine files in this directory are optional for CIME builds; if none
    # is found, fall back to common.cmake which sets SCREAM_MACHINE and
    # includes the EKAT machine file if available.
    set(SCREAM_MACH_FILE_ROOT ${CMAKE_SOURCE_DIR}/eamxx/cmake/machine-files)
    if (EXISTS ${SCREAM_MACH_FILE_ROOT}/${MACH}-${COMPILER}.cmake)
      include(${SCREAM_MACH_FILE_ROOT}/${MACH}-${COMPILER}.cmake)
    elseif (EXISTS ${SCREAM_MACH_FILE_ROOT}/${MACH}.cmake)
      include(${SCREAM_MACH_FILE_ROOT}/${MACH}.cmake)
    else()
      include(${SCREAM_MACH_FILE_ROOT}/common.cmake)
      common_setup()
    endif()

    # The machine files may enable kokkos stuff we don't want
    if (NOT compile_threaded)
      set(Kokkos_ENABLE_OPENMP FALSE)
    endif()

    add_subdirectory("eamxx")
  endif()

endfunction(build_eamxx)
