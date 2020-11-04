macro (CreateScorpioTarget CREATE_FLIB)

  # Some sanity checks
  if (NOT SCREAM_CIME_BUILD)
    message (FATAL_ERROR "Error! Scorpio.cmake currently only works in a CIME build.")
  endif ()

  if (NOT DEFINED INSTALL_SHAREDPATH)
    message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
  endif ()

  # If c lib is requested (and we didn't already parsed this script), create interface lib
  if (NOT TARGET pioc)
    # Get GPTL as a target
    include (tpls/GPTL)
    CreateGPTLTarget()

    # Get Netcdf libs
    include (tpls/GetNetcdfLibs)
    GetNetcdfLibs()

    # Look for pioc lib in the lib subdirectory of the one stored in INSTALL_SHAREDPATH (set by CIME)
    find_library(SCORPIO_C_LIB pioc REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # Create the interface library that scream targets can link to
    add_library(pioc UNKNOWN IMPORTED GLOBAL)
    set_target_properties(pioc PROPERTIES IMPORTED_LOCATION "${SCORPIO_C_LIB}")
    target_link_libraries(pioc INTERFACE "${netcdf_c_lib}")
    set_target_properties(pioc PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_SHAREDPATH}/include)
    if (pnetcdf_lib)
      target_link_libraries(pioc INTERFACE "${pnetcdf_lib}")
    endif ()
    target_link_libraries(pioc INTERFACE gptl)
  endif ()

  # If f lib is requested (and we didn't already parsed this script), create interface lib
  if (${CREATE_FLIB} AND NOT TARGET piof)
    # Look for piof lib in the lib subdirectory of the one stored in INSTALL_SHAREDPATH (set by CIME)
    find_library(SCORPIO_F_LIB piof REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # Create the imported library that scream targets can link to
    add_library(piof UNKNOWN IMPORTED GLOBAL)
    set_target_properties(piof PROPERTIES IMPORTED_LOCATION "${SCORPIO_F_LIB}")
    target_link_libraries(piof INTERFACE "${netcdf_f_lib};pioc")
    set_target_properties(piof PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_SHAREDPATH}/include)
  endif ()

endmacro()
