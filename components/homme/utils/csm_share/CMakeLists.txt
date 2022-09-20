# Make sure we create the target once
if (NOT TARGET csm_share)
  set (SHARE_DIR ${HOMME_SOURCE_DIR}/../../share/util)

  set (CSM_SHARE_SOURCES)

  # First, handle the files that need to go through genf90

  # Set variables needed for processing genf90 templates
  set(CIMEROOT ${HOMME_SOURCE_DIR}/../../cime)
  list(APPEND CMAKE_MODULE_PATH ${CIMEROOT}/CIME/non_py/src/CMake)
  set(GENF90 ${CIMEROOT}/CIME/non_py/externals/genf90/genf90.pl)
  set(ENABLE_GENF90 True)
  include(genf90_utils)
  include(Sourcelist_utils)

  # GENF90_SOURCE lists source files we will need to run through the genf90 perl script
  set (GENF90_SOURCE
    ${SHARE_DIR}/shr_infnan_mod.F90.in
  )

  # Process genf90 template files. For each src file, this adds a custom command
  # (and hence target) that needs to be built from the genf90 template files
  foreach (SRC_FILE ${GENF90_SOURCE})
    string(REPLACE ".in" "" SRC_FILE_STRIPPED ${SRC_FILE})
    get_filename_component(BASENAME ${SRC_FILE_STRIPPED} NAME)
    set(SRC_FILE_OUT "${CMAKE_CURRENT_BINARY_DIR}/${BASENAME}")
    add_custom_command (
      OUTPUT ${SRC_FILE_OUT}
      COMMAND ${GENF90} ${SRC_FILE} > ${SRC_FILE_OUT}
      DEPENDS ${SRC_FILE} genf90)
    list(APPEND CSM_SHARE_SOURCES ${SRC_FILE_OUT})
  endforeach ()

  # "normal" sources
  list (APPEND CSM_SHARE_SOURCES
    ${SHARE_DIR}/shr_abort_mod.F90
    ${SHARE_DIR}/shr_const_mod.F90
    ${SHARE_DIR}/shr_file_mod.F90
    ${SHARE_DIR}/shr_kind_mod.F90
    ${SHARE_DIR}/shr_log_mod.F90
    ${SHARE_DIR}/shr_mpi_mod.F90
    ${SHARE_DIR}/shr_spfn_mod.F90
    ${SHARE_DIR}/shr_strconvert_mod.F90
    ${SHARE_DIR}/shr_sys_mod.F90
    ${SHARE_DIR}/shr_vmath_mod.F90
  )

  add_library(csm_share ${CSM_SHARE_SOURCES})
  if (${CMAKE_SYSTEM} MATCHES "Linux")
    target_compile_definitions(csm_share PUBLIC LINUX)
  endif()
  set_target_properties(csm_share PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
  target_include_directories(csm_share PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/modules)

endif()
