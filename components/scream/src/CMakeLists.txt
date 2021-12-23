add_subdirectory(share)

if (NOT DEFINED ENV{SCREAM_FAKE_ONLY})
  add_subdirectory(dynamics)
  add_subdirectory(physics)
  add_subdirectory(control)
  if (PROJECT_NAME STREQUAL "E3SM")
    add_subdirectory(mct_coupling)
  endif()
elseif (SCREAM_FAKE_NEED_P3)
  # Cannot test diff handling without p3
  add_subdirectory(physics)
endif()
