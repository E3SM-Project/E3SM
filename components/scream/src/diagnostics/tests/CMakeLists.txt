# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  include(ScreamUtils)

  set( NEED_LIBS scream_share diagnostics )

  # Test potential temperature diagnostic
  CreateUnitTest(potential_temperature "potential_temperature_test.cpp" "${NEED_LIBS}" LABELS "diagnostics" )

endif()
