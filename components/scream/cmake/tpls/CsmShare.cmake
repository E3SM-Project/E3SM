macro (CreateCsmShareTarget)

  # Some sanity checks
  if (NOT CIME_BUILD)
    message (FATAL_ERROR "Error! CsmShare.cmake currently only works in a CIME build.")
  endif ()

  if (NOT DEFINED CSM_SHARE)
    message (FATAL_ERROR "Error! The cmake variable 'CSM_SHARE' is not defined.")
  endif ()

  # If we didn't already parse this script, create interface lib
  if (NOT TARGET scream_csm_share)
    find_library(CSM_SHARE_LIB csm_share REQUIRED PATHS ${CSM_SHARE})

    # Create the imported library that scream targets can link to
    add_library(scream_csm_share UNKNOWN IMPORTED GLOBAL)
    set_target_properties(scream_csm_share PROPERTIES IMPORTED_LOCATION ${CSM_SHARE_LIB})
    set_target_properties(scream_csm_share PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${CSM_SHARE})

    include(Scorpio)
    CreateScorpioTarget(TRUE)
    target_link_libraries(scream_csm_share INTERFACE scream_piof)
  endif ()
endmacro()
