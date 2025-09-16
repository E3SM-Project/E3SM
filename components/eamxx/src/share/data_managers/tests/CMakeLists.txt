# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test field manager
  CreateUnitTest(field_mgr "field_mgr_tests.cpp" LIBS eamxx_data_managers)
endif()
