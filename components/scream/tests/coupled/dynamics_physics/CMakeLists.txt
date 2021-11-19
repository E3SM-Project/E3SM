## Note that currently the HOMME dycore only works in double precision
if (SCREAM_DOUBLE_PRECISION)
  if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
    add_subdirectory(homme_shoc_cld_p3_rrtmgp)
    # Only compile the test with SPA if the SPA data directory has been defined.
    if (DEFINED SCREAM_SPA_DATA_DIR)
      add_subdirectory(homme_shoc_cld_spa_p3_rrtmgp)
    endif()
  endif()
endif()

