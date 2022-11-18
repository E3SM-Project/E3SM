if (NOT SCREAM_BASELINES_ONLY)
  include(ScreamUtils)

  CreateUnitTest(nudging "nudging_tests.cpp" "${NEED_LIBS}" LABELS "physics_nudging" )

endif()
