# This cmake utility allows you to build cprnc just by adding
#
#   Include(BuildCprnc)
#   BuildCprnc()
#
# to your CMakeLists.txt.

include (EkatUtils)
macro(BuildCprnc)

  # TODO: handle this more carefully and more gracefully in the future
  # TODO: For now, it is just a hack to get going...
  # find cprnc defined in machine entries
  set(CCSM_CPRNC $ENV{CCSM_CPRNC})
  if(EXISTS "${CCSM_CPRNC}")
    message(STATUS "Path ${CCSM_CPRNC} exists, so we will use it")
    set(CPRNC_BINARY ${CCSM_CPRNC} CACHE INTERNAL "")
    configure_file (${SCREAM_BASE_DIR}/cmake/CprncTest.cmake.in
                    ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake @ONLY)
  else()
    if (NOT "${CCSM_CPRNC}" STREQUAL "")
      message(WARNING "Path ${CCSM_CPRNC} does not exist, so we will try to build it")
    endif()
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
      # Let's make sure cprnc can at least find the FindNetCDF.cmake
      # module that scorpio has.
      list(APPEND CMAKE_MODULE_PATH ${SRC_ROOT}/externals/scorpio/cmake)
      add_subdirectory(${SRC_ROOT}/cime/CIME/non_py/cprnc ${BLDROOT})
      EkatDisableAllWarning(cprnc)

      set(CPRNC_BINARY ${BLDROOT}/cprnc CACHE INTERNAL "")

      configure_file (${SCREAM_BASE_DIR}/cmake/CprncTest.cmake.in
                      ${CMAKE_BINARY_DIR}/bin/CprncTest.cmake @ONLY)
    endif()
  endif()
endmacro()
