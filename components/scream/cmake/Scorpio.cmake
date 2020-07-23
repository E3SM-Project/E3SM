macro (CreateScorpioTarget CREATE_FLIB)

  # Some sanity checks
  if (NOT CIME_BUILD)
    message (FATAL_ERROR "Error! Scorpio.cmake currently only works in a CIME build.")
  endif ()

  if (NOT DEFINED INSTALL_SHAREDPATH)
    message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
  endif ()

  # If c lib is requested (and we didn't already parsed this script), create interface lib
  if (NOT TARGET scream_pioc)
    # Get GPTL as a target
    include (GPTL)
    CreateGPTLTarget()

    # Get Netcdf libs
    include (GetNetcdfLibs)
    GetNetcdfLibs()

    # Look for pioc lib in the lib subdirectory of the one stored in INSTALL_SHAREDPATH (set by CIME)
    find_library(SCORPIO_C_LIB pioc REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # Create the interface library that scream targets can link to
    add_library(scream_pioc UNKNOWN IMPORTED GLOBAL)
    set_target_properties(scream_pioc PROPERTIES IMPORTED_LOCATION "${SCORPIO_C_LIB}")
    target_link_libraries(scream_pioc INTERFACE "scream_gptl;${netcdf_c_lib}")
    set_target_properties(scream_pioc PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_SHAREDPATH}/include)
    if (pnetcdf_lib)
      target_link_libraries(scream_pioc INTERFACE "${pnetcdf_lib}")
    endif ()
  endif ()

  # If f lib is requested (and we didn't already parsed this script), create interface lib
  if (${CREATE_FLIB} AND NOT TARGET scream_piof)
    # Look for piof lib in the lib subdirectory of the one stored in INSTALL_SHAREDPATH (set by CIME)
    find_library(SCORPIO_F_LIB piof REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # Create the imported library that scream targets can link to
    add_library(scream_piof UNKNOWN IMPORTED GLOBAL)
    set_target_properties(scream_piof PROPERTIES IMPORTED_LOCATION "${SCORPIO_F_LIB}")
    target_link_libraries(scream_piof INTERFACE "${netcdf_f_lib};scream_pioc")
    set_target_properties(scream_piof PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_SHAREDPATH}/include)
  endif ()

endmacro()
