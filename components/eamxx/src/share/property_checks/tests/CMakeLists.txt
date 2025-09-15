# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test property checks
  CreateUnitTest(property_checks "property_checks.cpp")
endif()
