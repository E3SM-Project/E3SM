macro (CreateGPTLTarget)

  # Some sanity checks
  if (NOT CIME_BUILD)
    message (FATAL_ERROR "Error! GPTL.cmake currently only works in a CIME build.")
  endif ()

  if (NOT DEFINED INSTALL_SHAREDPATH)
    message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
  endif ()

  # If we didn't already parse this script, create interface lib
  if (NOT TARGET scream_gptl)
    find_library(GPTL_LIB gptl REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # Create the interface library that scream targets can link to
    add_library(scream_gptl INTERFACE IMPORTED GLOBAL)
    target_link_libraries (scream_gptl INTERFACE ${GPTL_LIB})

    list(APPEND SCREAM_TPL_LIBRARIES scream_gptl)
    list(APPEND SCREAM_TPL_INCLUDE_DIRS ${INSTALL_SHAREDPATH}/include)
  endif ()
endmacro()
