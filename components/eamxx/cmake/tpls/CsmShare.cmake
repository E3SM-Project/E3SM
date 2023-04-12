macro (CreateCsmShareTarget)
  if (TARGET csm_share)
    message (FATAL_ERROR "Error! The target csm_share already exists!")
  endif()

  if (SCREAM_CIME_BUILD)
    # Some sanity checks
    if (NOT DEFINED INSTALL_SHAREDPATH)
      message (FATAL_ERROR "Error! The cmake variable 'INSTALL_SHAREDPATH' is not defined.")
    endif ()
    if (NOT DEFINED COMP_INTERFACE)
      message (FATAL_ERROR "Error! The cmake variable 'COMP_INTERFACE' is not defined.")
    endif ()
    if (NOT DEFINED NINST_VALUE)
      message (FATAL_ERROR "Error! The cmake variable 'NINST_VALUE' is not defined.")
    endif ()

    # If we didn't already parse this script, create imported target
    if (NOT TARGET csm_share)

      # Build the name of the path where libcsm_share should be located
      if (USE_ESMF_LIB)
        set(ESMFDIR "esmf")
      else()
        set(ESMFDIR "noesmf")
      endif()
      set(CSM_SHARE "${INSTALL_SHAREDPATH}/${COMP_INTERFACE}/${ESMFDIR}/${NINST_VALUE}/csm_share")

      # Look for libcsm_share in the complex path we built above
      find_library(CSM_SHARE_LIB csm_share REQUIRED PATHS ${CSM_SHARE})

      # Create the interface library, and set target properties
      add_library (csm_share INTERFACE)
      target_link_libraries (csm_share INTERFACE ${CSM_SHARE_LIB})
      target_include_directories(csm_share INTERFACE ${CSM_SHARE})

      # Link against piof
      target_link_libraries(csm_share INTERFACE piof)
    endif ()
  else()
    # Build csm_share library manually

    # Set variables needed for processing genf90 templates
    set(CIMEROOT ${SCREAM_BASE_DIR}/../../cime)
    list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/CIME/non_py/src/CMake)
    set(GENF90 ${CIMEROOT}/CIME/non_py/externals/genf90/genf90.pl)
    set(ENABLE_GENF90 True)
    include(genf90_utils)
    include(Sourcelist_utils)

    # GENF90_SOURCE lists source files we will need to run through the genf90 perl script
    set (GENF90_SOURCE 
        ${SCREAM_BASE_DIR}/../../share/util/shr_infnan_mod.F90.in
        ${SCREAM_BASE_DIR}/../../share/util/shr_assert_mod.F90.in
    )
    # FORTRAN_SOURCE lists the source files we want to build that do NOT need to be run through the genf90
    # perl script. We will append to this list below with our processed genf90 files.
    set (FORTRAN_SOURCE
        ${SCREAM_BASE_DIR}/../../share/util/shr_abort_mod.F90
        ${SCREAM_BASE_DIR}/../../share/util/shr_const_mod.F90
        ${SCREAM_BASE_DIR}/../../share/util/shr_file_mod.F90
        ${SCREAM_BASE_DIR}/../../share/util/shr_kind_mod.F90
        ${SCREAM_BASE_DIR}/../../share/util/shr_log_mod.F90
        ${SCREAM_BASE_DIR}/../../share/util/shr_mpi_mod.F90
        ${SCREAM_BASE_DIR}/../../share/util/shr_orb_mod.F90
        ${SCREAM_BASE_DIR}/../../share/util/shr_strconvert_mod.F90
        ${SCREAM_BASE_DIR}/../../share/util/shr_sys_mod.F90
    )
    # Process genf90 template files. This adds a custom command (and hence target) for each f90 source
    # that needs to be built from the genf90 template files listed in GENF90_SOURCE.
    foreach (SRC_FILE ${GENF90_SOURCE})
        string(REPLACE ".in" "" SRC_FILE_STRIPPED ${SRC_FILE})
        get_filename_component(BASENAME ${SRC_FILE_STRIPPED} NAME)
        set(SRC_FILE_OUT "${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}")
        add_custom_command (
            OUTPUT ${SRC_FILE_OUT}
            COMMAND ${GENF90} ${SRC_FILE} > ${SRC_FILE_OUT}
            DEPENDS ${SRC_FILE} genf90)
        list(APPEND FORTRAN_SOURCE ${SRC_FILE_OUT})
    endforeach ()
    set(share_sources ${FORTRAN_SOURCE})

    add_library(csm_share ${share_sources})

    # These CPP macros are needed in shr_infnan_mod
    target_compile_definitions(csm_share PUBLIC
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CXX_COMPILER_ID:GNU>>:CPRGNU>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CXX_COMPILER_ID:Intel>>:CPRINTEL>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CXX_COMPILER_ID:Clang>>:CPRCRAY>
      $<$<AND:$<COMPILE_LANGUAGE:Fortran>,$<CXX_COMPILER_ID:clang>>:CPRCRAY>)


    if (${CMAKE_SYSTEM} MATCHES "Linux")
      target_compile_definitions(csm_share PUBLIC LINUX)
    endif()
    set_target_properties(csm_share PROPERTIES
      Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
    target_include_directories(csm_share PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/modules)
  endif ()
endmacro()
