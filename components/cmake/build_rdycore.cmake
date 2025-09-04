function(build_rdycore)

  if (COMP_NAMES MATCHES ".*rdycore.*")

    # set PETSc environment variables by compiler and machine
    if (NOT COMPILER MATCHES "gnu")
      message(FATAL_ERROR "Unsupported compiler: ${COMPILER} -- RDycore requires GNU compilers.")
    endif()
    set(ENV{PETSC_DIR} "/global/cfs/projectdirs/m4267/petsc/petsc_v3.23.0")
    if (MACH STREQUAL "pm-cpu" OR MACH STREQUAL "pm-gpu")
      set(ENV{PETSC_ARCH} "${MACH}-opt-32bit-gcc-13-2-1-95934b0d393")
    else()
      message(FATAL_ERROR "RDycore is not supported on the machine '${MACH}'")
    endif()

    # use PETSc's bundled version of Kokkos
    set(ENV{Kokkos_ROOT} "$ENV{PETSC_DIR}/$ENV{PETSC_ARCH}")

    message(STATUS "Found rdycore component")

    include(${CMAKE_SOURCE_DIR}/cmake/common_setup.cmake)

    # Override some RDycore options to simplify the build process.
    set(ENABLE_DRIVER OFF CACHE BOOL "disable RDycore driver" FORCE)
    set(ENABLE_TESTS OFF CACHE BOOL "disable RDycore tests" FORCE)

    # PROJECT_SOURCE_DIR is E3SM/components
    set(EXTERNALS_SOURCE_DIR "${PROJECT_SOURCE_DIR}/../externals")
    add_subdirectory(${EXTERNALS_SOURCE_DIR}/rdycore ${CMAKE_BINARY_DIR}/externals/rdycore)

    # add the directory containing rdycore.mod as a global include directory
    get_target_property(RDYCORE_MOD_DIR rdycore_f90 BINARY_DIR)
    include_directories(${RDYCORE_MOD_DIR})
  endif()

endfunction(build_rdycore)
