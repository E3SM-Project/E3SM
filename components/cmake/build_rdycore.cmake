function(build_rdycore)

  if (COMP_NAMES MATCHES ".*rdycore.*")

    # set PETSc include directories
    # FIXME: this isn't the way!
    include_directories("/global/cfs/projectdirs/m4267/petsc/petsc_v3.22.0/include")
    include_directories("/global/cfs/projectdirs/m4267/petsc/petsc_v3.22.0/pm-cpu-hdf5_1_14_3-opt-64bit-gcc-11-2-0-v3.22.0/include")

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
