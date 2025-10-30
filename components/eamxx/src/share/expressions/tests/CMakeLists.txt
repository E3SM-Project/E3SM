# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test expression templates
  CreateUnitTest(expressions "expressions.cpp" LIBS eamxx_field_math)
endif()
