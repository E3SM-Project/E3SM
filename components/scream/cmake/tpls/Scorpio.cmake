# If this is a CIME build, create IMPORTED target to wrap scorpio libs.
# Otherwise, simply add scorpio subdirectory.
set (E3SM_EXTERNALS_DIR ${CMAKE_CURRENT_LIST_DIR}/../../../../externals CACHE INTERNAL "")

set (SCREAM_TPLS_MODULE_DIR ${CMAKE_CURRENT_LIST_DIR} CACHE INTERNAL "")

macro (CreateScorpioTarget CREATE_FLIB)

  # If we didn't already parsed this script, proceed
  if (NOT TARGET pioc)

    if (CIME_BUILD)
      # For CIME builds, we simply wrap the already built pioc/piof libs into a cmake target
      if (NOT DEFINED INSTALL_SHAREDPATH)
        message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
      endif ()

      set (SCORPIO_LIB_DIR ${INSTALL_SHAREDPATH}/lib)
      set (SCORPIO_INC_DIR ${INSTALL_SHAREDPATH}/include)

      # Look for pioc in INSTALL_SHAREDPATH/lib
      find_library(SCORPIO_C_LIB pioc REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

      # Create imported target
      add_library(pioc UNKNOWN IMPORTED GLOBAL)
      set_target_properties(pioc PROPERTIES
                IMPORTED_LOCATION "${SCORPIO_C_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_SHAREDPATH}/include)

      # Look for pioc deps, and attach them to the pioc target, so that cmake will
      # propagate them to any downstream target linking against pioc
      include (${SCREAM_TPLS_MODULE_DIR}/GPTL.cmake)
      CreateGPTLTarget()
      include (${SCREAM_TPLS_MODULE_DIR}/GetNetcdfLibs.cmake)
      GetNetcdfLibs()
      target_link_libraries(pioc INTERFACE "gptl;${netcdf_c_lib}")
      if (pnetcdf_lib)
        target_link_libraries(pioc INTERFACE "${pnetcdf_lib}")
      endif ()

      # If f lib is requested (and we didn't already parsed this script), proceed
      if (CREATE_FLIB AND NOT TARGET piof)
        # Look for piof lib in INSTALL_SHAREDPATH/lib
        find_library(SCORPIO_F_LIB piof REQUIRED PATHS ${INSTALL_SHAREDPATH}/lib)

        # Create the imported library that scream targets can link to
        add_library(piof UNKNOWN IMPORTED GLOBAL)
        set_target_properties(piof PROPERTIES
                IMPORTED_LOCATION "${SCORPIO_F_LIB}"
                INTERFACE_INCLUDE_DIRECTORIES ${INSTALL_SHAREDPATH}/include)
        # Link pioc and netcdf-fortran, so cmake will propagate them to any downstream
        # target linking against piof
        target_link_libraries(piof INTERFACE "${netcdf_f_lib};pioc")
      endif ()
    else ()
      # Not a CIME build. Add scorpio as a subdir
      add_subdirectory (${E3SM_EXTERNALS_DIR}/scorpio ${CMAKE_BINARY_DIR}/externals/scorpio)
      EkatDisableAllWarning(pioc)
      EkatDisableAllWarning(gptl)
      if (CREATE_FLIB)
        EkatDisableAllWarning(pioc)
      endif()
    endif()
  endif ()

endmacro()
