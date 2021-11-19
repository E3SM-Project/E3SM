# Note RRMTMGP only works with double-precision, so only compile tests for DP
if (SCREAM_DOUBLE_PRECISION)
  add_subdirectory(shoc_cld_p3_rrtmgp)
  # Only compile the test with SPA if the SPA data directory has been defined.
  if (DEFINED SCREAM_SPA_DATA_DIR)
    add_subdirectory(shoc_cld_spa_p3_rrtmgp)
  endif()
endif()

