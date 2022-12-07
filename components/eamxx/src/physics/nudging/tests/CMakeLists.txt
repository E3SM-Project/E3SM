if (NOT SCREAM_BASELINES_ONLY)
  include(ScreamUtils)

  set( NEED_LIBS scream_share nudging physics_share scream_io)

  CreateUnitTest(nudging_tests "nudging_tests.cpp" "${NEED_LIBS}" LABELS "physics_nudging" )

endif()
