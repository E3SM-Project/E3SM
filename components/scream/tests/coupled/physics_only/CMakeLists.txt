# Note RRMTMGP only works with double-precision, so only compile tests for DP
if (SCREAM_DOUBLE_PRECISION)
  add_subdirectory(shoc_cld_p3_rrtmgp)
  add_subdirectory(shoc_cld_spa_p3_rrtmgp)
endif()

add_subdirectory (atm_proc_subcycling)
