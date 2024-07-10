add_subdirectory(physics_only)
if (NOT "${SCREAM_DYNAMICS_DYCORE}" STREQUAL "NONE")
  add_subdirectory(dynamics_physics)
endif()
