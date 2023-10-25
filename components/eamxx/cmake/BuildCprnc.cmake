# This cmake utility allows you to build cprnc just by adding
#
#   Include(BuildCprnc)
#   BuildCprnc()
#
# to your CMakeLists.txt.

include (EkatUtils)
macro(BuildCprnc)

  # Make sure this is built only once
  if (NOT TARGET cprnc)
    if (SCREAM_CIME_BUILD)
      string (CONCAT MSG
              "WARNING! By default, scream should not build tests in a CIME build,\n"
              "and cprnc should only be built by scream in case tests are enabled.\n"
              "If you explicitly requested tests to be on in a CIME build,\n"
              "then you can discard this warning. Otherwise, please, contact developers.\n")
      message("${MSG}")
    endif()
    set(BLDROOT ${PROJECT_BINARY_DIR}/externals/cprnc)
    file(WRITE ${BLDROOT}/Macros.cmake
      "
      set(SCC ${CMAKE_C_COMPILER})
      set(SFC ${CMAKE_Fortran_COMPILER})
      set(FFLAGS \"${CMAKE_Fortran_FLAGS}\")
      set(NETCDF_PATH ${NetCDF_Fortran_PATH})
      "
    )
    set(SRC_ROOT ${SCREAM_BASE_DIR}/../..)
    add_subdirectory(${SRC_ROOT}/cime/CIME/non_py/cprnc ${BLDROOT})
    EkatDisableAllWarning(cprnc)

    set(CPRNC_BINARY ${BLDROOT}/cprnc CACHE INTERNAL "")

    configure_file (${SCREAM_BASE_DIR}/cmake/CprncTest.cmake.in
                    ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake @ONLY)
  endif()
endmacro()
