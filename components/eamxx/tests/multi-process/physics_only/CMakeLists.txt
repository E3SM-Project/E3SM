# Note RRMTMGP only works with double-precision, so only compile tests for DP
if (SCREAM_DOUBLE_PRECISION)
  add_subdirectory(shoc_cld_p3_rrtmgp)
  add_subdirectory(shoc_cld_spa_p3_rrtmgp)
    if (SCREAM_ENABLE_MAM)
      add_subdirectory(mam/optics_rrtmgp)
      add_subdirectory(mam/shoc_aci)
    endif()
endif()

add_subdirectory (atm_proc_subcycling)
add_subdirectory (shoc_p3_nudging)
