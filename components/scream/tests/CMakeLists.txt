# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT ${SCREAM_BASELINES_ONLY})
  if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
    add_subdirectory(scream_homme_dyn_ut_nlev72_qsize4)
    add_subdirectory(scream_homme_dyn_p3)
  endif()

  add_subdirectory(scream_p3)
  add_subdirectory(scream_zm)
  add_subdirectory(scream_p3_shoc)
  add_subdirectory(rrtmgp)
endif()
