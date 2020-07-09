# Use Macros.cmake to get info on netcdf paths.
# Note: the inputs are supposed to be *the name* of the variables storing the result
# Note: Keep this a FUNCTION, not a MACRO, to avoit polluting the calling scope
#       with all the stuff from Macros.cmake
function (GetNetcdfLibs)
  # Sanity check
  if (NOT CIME_BUILD)
    message (FATAL_ERROR "Error! Do not call 'GetNetcdfPaths' in a non-CIME build.\n")
  endif ()

  # Load variables set by CIME
  include(${CASEROOT}/Macros.cmake)

  # Pnetcdf is optional, and only if not running serial
  if (NOT MPILIB STREQUAL mpi-serial)
    if (PNETCDF_PATH)
      find_library(pnetcdf_lib pnetcdf REQUIRED PATHS ${PNETCDF_PATH}/lib)
      set (pnetcdf_lib ${pnetcdf_lib} PARENT_SCOPE)
    endif()
  endif()

  if (NETCDF_C_PATH)
    # Sanity checks
    if (NOT NETCDF_FORTRAN_PATH)
      message(FATAL_ERROR "NETCDF_C_PATH specified without NETCDF_FORTRAN_PATH")
    endif()
    if (NOT EXISTS ${NETCDF_C_PATH}/lib AND NOT EXISTS ${NETCDF_C_PATH}/lib64)
      message(FATAL_ERROR "NETCDF_C_PATH does not contain a lib or lib64 directory")
    endif ()
    if (NOT EXISTS ${NETCDF_FORTRAN_PATH}/lib AND NOT EXISTS ${NETCDF_FORTRAN_PATH}/lib64)
      message(FATAL_ERROR "NETCDF_FORTRAN_PATH does not contain a lib or lib64 directory")
    endif ()

    # Find the libraries
    find_library(netcdf_c_lib netcdf  REQUIRED PATHS ${NETCDF_C_PATH}/lib ${NETCDF_C_PATH}/lib64)
    find_library(netcdf_f_lib netcdff REQUIRED PATHS ${NETCDF_FORTRAN_PATH}/lib ${NETCDF_FORTRAN_PATH}/lib64)

  elseif (NETCDF_FORTRAN_PATH)
    message(FATAL_ERROR "NETCDF_FORTRAN_PATH specified without NETCDF_C_PATH")
  elseif (NETCDF_PATH)

    # Sanity checks
    if (NOT EXISTS ${NETCDF_PATH}/lib AND NOT EXISTS ${NETCDF_PATH}/lib64)
      message(FATAL_ERROR "NETCDF_PATH does not contain a lib or lib64 directory")
    endif ()

    find_library(netcdf_c_lib netcdf  REQUIRED PATHS ${NETCDF_PATH}/lib ${NETCDF_PATH}/lib64)
    find_library(netcdf_f_lib netcdff REQUIRED PATHS ${NETCDF_PATH}/lib ${NETCDF_PATH}/lib64)
  else()
    message(FATAL_ERROR "NETCDF not found: Define NETCDF_PATH or NETCDF_C_PATH and NETCDF_FORTRAN_PATH in config_machines.xml or config_compilers.xml")
  endif()
  set (netcdf_c_lib ${netcdf_c_lib} PARENT_SCOPE)
  set (netcdf_f_lib ${netcdf_f_lib} PARENT_SCOPE)
endfunction ()

macro (CreateScorpioTarget CLIB FLIB)

  # Some sanity checks
  if (NOT CIME_BUILD)
    message (FATAL_ERROR "Error! Scorpio.cmake currently only works in a CIME build.")
  endif ()

  if (NOT DEFINED INSTALL_SHAREDPATH)
    message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
  endif ()

  # If c lib is requested (and we didn't already parsed this script), create interface lib
  if (${CLIB} AND NOT TARGET scorpio_c)
    # Get GPTL as a target
    include (GPTL)
    CreateGPTLTarget()

    # Get Netcdf libs
    GetNetcdfLibs()

    # Look for pioc lib in the lib subdirectory of the one stored in INSTALL_SHAREDPATH (set by CIME)
    find_library(SCORPIO_C_LIB pioc REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # Create the interface library that scream targets can link to
    add_library(scorpio_c INTERFACE IMPORTED GLOBAL)
    # set_target_properties (scorpio_c PROPERTIES IMPORTED_LOCATION ${SCORPIO_C_LIB})
    target_link_libraries(scorpio_c INTERFACE "${SCORPIO_C_LIB};scream_gptl;${netcdf_c_lib}")
    if (pnetcdf_lib)
      target_link_libraries(scorpio_c INTERFACE "${pnetcdf_lib}")
    endif ()

    # Update the list of scream tpls
    list(APPEND SCREAM_TPL_LIBRARIES scorpio_c)
    list(APPEND SCREAM_TPL_INCLUDE_DIRS ${INSTALL_SHAREDPATH}/include)

    # If f lib is requested (and we didn't already parsed this script), create interface lib
    if (${FLIB} AND NOT TARGET scorpio_f)
      # Look for piof lib in the lib subdirectory of the one stored in INSTALL_SHAREDPATH (set by CIME)
      find_library(SCORPIO_F_LIB piof REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

      # Create the interface library that scream targets can link to
      add_library(scorpio_f INTERFACE IMPORTED GLOBAL)
      # set_target_properties (scorpio_f PROPERTIES IMPORTED_LOCATION ${SCORPIO_F_LIB})
      target_link_libraries(scorpio_f INTERFACE ${SCORPIO_F_LIB} ${netecdf_f_lib} ${SCORPIO_C_LIB})
      if (pnetcdf_lib)
        target_link_libraries(scorpio_f INTERFACE "${pnetcdf_lib}")
      endif ()

      # Update the list of scream tpls
      list(APPEND SCREAM_TPL_LIBRARIES scorpio_f)
      list(APPEND SCREAM_TPL_INCLUDE_DIRS ${INSTALL_SHAREDPATH}/include)
    endif ()
  endif ()

endmacro()
