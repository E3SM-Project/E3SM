macro (CreateMctTarget)

  # Some sanity checks
  if (NOT SCREAM_CIME_BUILD)
    message (FATAL_ERROR "Error! You should need the mct target only in CIME builds")
  endif ()
  if (NOT DEFINED INSTALL_SHAREDPATH)
    message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
  endif ()

  if (TARGET mct)
    # We should not call this macro twice
    message (FATAL_ERROR "The mct target was already created!")
  endif()

  # Look for libmct in INSTALL_SHAREDPATH/lib
  find_library(MCT_LIB mct REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

  # Create the interface library, and set target properties
  add_library(mct INTERFACE)
  target_link_libraries(mct INTERFACE ${MCT_LIB})
  target_include_directories(mct INTERFACE ${INSTALL_SHAREDPATH}/include)

  # Link against csm_share
  target_link_libraries(mct INTERFACE csm_share)
endmacro()
