if ("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "HOMME")
  add_subdirectory(homme)
  set (SCREAM_HAS_HOMME TRUE CACHE INTERNAL "Signal that Homme is compiled")
elseif("${SCREAM_DYNAMICS_DYCORE}" STREQUAL "NONE")
  message ("SCREAM_DYNAMICS_DYCORE set to 'NONE'. Scream won't enable any test/code that needs dynamics.\n")
else()
  message (FATAL_ERROR "The dynamics dycore '${SCREAM_DYNAMICS_DYCORE}' is currently not supported."
                       "If you wish to disable dynamics, do not specify anything for SCREAM_DYNAMICS_DYCORE, or set it to 'NONE'.")
endif()
