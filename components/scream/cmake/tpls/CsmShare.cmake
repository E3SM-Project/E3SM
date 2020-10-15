set (SCREAM_TPLS_MODULE_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE INTERNAL "")
set (SCREAM_CIME_SRC_DIR ${CMAKE_CURRENT_LIST_DIR}/../../../../cime/src)

macro (CreateCsmShareTarget)

  # Some sanity checks
  if (CIME_BUILD)
    if (NOT DEFINED CSM_SHARE)
      message (FATAL_ERROR "Error! The cmake variable 'CSM_SHARE' is not defined.")
    endif ()

    # If we didn't already parse this script, create imported target
    if (NOT TARGET csm_share)
      find_library(CSM_SHARE_LIB csm_share REQUIRED PATHS ${CSM_SHARE})

      # Create the imported library that scream targets can link to
      add_library(csm_share UNKNOWN IMPORTED GLOBAL)
      set_target_properties(csm_share PROPERTIES IMPORTED_LOCATION ${CSM_SHARE_LIB})
      set_target_properties(csm_share PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${CSM_SHARE})

      include(${SCREAM_TPLS_MODULE_DIR}/Scorpio.cmake)
      CreateScorpioTarget(TRUE)
      target_link_libraries(csm_share INTERFACE piof)
    endif ()
  else ()
    message (FATAL_ERROR "Error! CsmShare not yet supported in standalone mode.")
  endif ()

endmacro()
