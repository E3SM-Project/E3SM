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
    # Look for pioc lib in the lib subdirectory of the one stored in INSTALL_SHAREDPATH (set by CIME)
    find_library(SCORPIO_C_LIB pioc REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # # Create the interface library that scream targets can link to
    # add_library(scorpio_c INTERFACE)
    # target_link_libraries(scorpio_c INTERFACE ${SCORPIO_C_LIB})
    # target_include_directories(scorpio_c INTERFACE ${INSTALL_SHAREDPATH}/include)

    # # Update the list of scream tpls
    # list(APPEND SCREAM_TPL_LIBRARIES scorpio_c)
    # list(APPEND SCREAM_TPL_INCLUDE_DIRS ${INSTALL_SHAREDPATH}/include)
  endif ()

  # If f lib is requested (and we didn't already parsed this script), create interface lib
  if (${FLIB} AND NOT TARGET scorpio_f)
    # Look for piof lib in the lib subdirectory of the one stored in INSTALL_SHAREDPATH (set by CIME)
    find_library(SCORPIO_F_LIB piof REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # # Create the interface library that scream targets can link to
    # add_library(scorpio_f INTERFACE)
    # target_link_libraries(scorpio_c INTERFACE ${SCORPIO_C_LIB})
    # target_include_directories(scorpio_c INTERFACE ${INSTALL_SHAREDPATH}/include)

    # # Update the list of scream tpls
    # list(APPEND SCREAM_TPL_LIBRARIES scorpio_f)
    # list(APPEND SCREAM_TPL_INCLUDE_DIRS ${INSTALL_SHAREDPATH}/include)
  endif ()

endmacro()
