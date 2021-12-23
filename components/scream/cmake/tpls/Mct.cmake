set (SCREAM_TPLS_MODULE_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE INTERNAL "")

macro (CreateMctTarget)

  # Some sanity checks
  if (NOT SCREAM_CIME_BUILD)
    message (FATAL_ERROR "Error! Mct.cmake currently only works in a CIME build.")
  endif ()
  if (NOT DEFINED INSTALL_SHAREDPATH)
    message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
  endif ()

  # If we didn't already parse this script, proceed
  if (NOT TARGET mct)
    # Look for libmct in INSTALL_SHAREDPATH/lib
    find_library(MCT_LIB mct REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

    # Create the imported library that scream targets can link to
    add_library(mct UNKNOWN IMPORTED GLOBAL)
    set_target_properties(mct PROPERTIES IMPORTED_LOCATION ${MCT_LIB})
    set_target_properties(mct PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_SHAREDPATH}/include)

    # Create the csm_share imported target, and link it to mct, so that cmake will correctly
    # attach it to any downstream target linking against mct
    include(${SCREAM_TPLS_MODULE_DIR}/CsmShare.cmake)
    CreateCsmShareTarget()
    target_link_libraries(mct INTERFACE csm_share)
  endif ()
endmacro()
