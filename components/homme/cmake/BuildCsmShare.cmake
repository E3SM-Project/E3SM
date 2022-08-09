# For HOMME standalone builds, we need to build an "economy"
# version of the csm_share library. This is only needed for
# standalone tests!
macro (BuildCsmShare)
  if (NOT HOMME_BUILD_EXECS)
    message (FATAL_ERROR "Error! You should build csm_share only for standalone tests.\n")
  endif()

  # Make sure we create the target once
  if (NOT TARGET csm_share)
    set (SHARE_DIR ${HOMME_SOURCE_DIR}/../../share/util)

    set (CSM_SHARE_SOURCES
      ${SHARE_DIR}/shr_const_mod.F90
      ${SHARE_DIR}/shr_file_mod.F90
      ${SHARE_DIR}/shr_kind_mod.F90
      ${SHARE_DIR}/shr_mpi_mod.F90
      ${SHARE_DIR}/shr_spfn_mod.F90
      ${SHARE_DIR}/shr_sys_mod.F90
      ${SHARE_DIR}/shr_vmath_mod.F90
    )
    add_library(csm_share ${CSM_SHARE_SOURCES})
  endif()
endmacro (BuildCsmShare)
