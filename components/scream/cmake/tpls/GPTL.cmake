macro (CreateGPTLTarget)

  # Some sanity checks
  if (NOT CIME_BUILD)
    message (FATAL_ERROR "Error! GPTL.cmake currently only works in a CIME build.")
  endif ()

  if (NOT DEFINED INSTALL_SHAREDPATH)
    message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
  endif ()

  # If we didn't already parse this script, create imported target
  if (NOT TARGET gptl)
    find_library(GPTL_LIB gptl REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # Create the interface library that scream targets can link to
    add_library(gptl UNKNOWN IMPORTED GLOBAL)
    set_target_properties(gptl  PROPERTIES IMPORTED_LOCATION ${GPTL_LIB})
    set_target_properties(gptl  PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_SHAREDPATH}/include)
  endif ()
endmacro()
