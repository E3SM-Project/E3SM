add_subdirectory(share)

if (NOT $ENV{SCREAM_FAKE_ONLY})
  add_subdirectory(dynamics)
  add_subdirectory(physics)
  add_subdirectory(control)
  if (PROJECT_NAME STREQUAL "E3SM")
    add_subdirectory(mct_coupling)
  endif()
elseif ($ENV{SCREAM_FORCE_RUN_DIFF})
  # Cannot test diff handling without p3
  add_subdirectory(physics)
endif()
