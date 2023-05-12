## Note that currently the HOMME dycore only works in double precision
if (SCREAM_DOUBLE_PRECISION)
  if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
    add_subdirectory(homme_shoc_cld_p3_rrtmgp)
    add_subdirectory(model_restart)
    add_subdirectory(homme_shoc_cld_spa_p3_rrtmgp)
    add_subdirectory(homme_shoc_cld_spa_p3_rrtmgp_128levels)
  endif()
endif()

