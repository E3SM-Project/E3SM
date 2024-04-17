## Note that currently the HOMME dycore only works in double precision
if (SCREAM_DOUBLE_PRECISION)
  if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
    add_subdirectory(homme_shoc_cld_p3_rrtmgp)
    add_subdirectory(homme_shoc_cld_p3_rrtmgp_pg2)
    add_subdirectory(model_restart)
    add_subdirectory(homme_shoc_cld_spa_p3_rrtmgp)
    add_subdirectory(homme_shoc_cld_spa_p3_rrtmgp_128levels)
    add_subdirectory(homme_shoc_cld_spa_p3_rrtmgp_dp)
    if (SCREAM_ENABLE_MAM)
      # Once the mam4xx aerosol microphysics AtmosphereProcess is running, the
      # corresponding test here needs to be reworked with valid aerosol
      # initial conditions.
      #add_subdirectory(homme_mam4xx_pg2)
      add_subdirectory(mam/homme_shoc_cld_p3_mam_optics_rrtmgp)
    endif()
  endif()
endif()

