# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  include(ScreamUtils)

  CreateUnitTest(dummy_diag "dummy_diag_test.cpp" "scream_share;diagnostics")
  # Test potential temperature diagnostic
  CreateUnitTest(potential_temperature "potential_temperature_test.cpp" "scream_share;diagnostics")

endif()
