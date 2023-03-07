# NOTE: if you have baseline-type tests, add the subdirectory OUTSIDE the following if statement
if (NOT SCREAM_BASELINES_ONLY)
  add_subdirectory(physics_only)
  if (NOT "${SCREAM_DYNAMICS_DYCORE}" STREQUAL "NONE")
    add_subdirectory(dynamics_physics)
  endif()
endif()
