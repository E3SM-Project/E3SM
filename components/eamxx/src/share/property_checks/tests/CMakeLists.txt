# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test FieldNaNCheck
  CreateUnitTest(field_nan_check "field_nan_check_tests.cpp"
    LIBS eamxx_property_checks)

  # Test FieldWithinIntervalCheck
  CreateUnitTest(field_within_interval_check "field_within_interval_check_tests.cpp"
    LIBS eamxx_property_checks)

  # Test FieldLowerBoundCheck
  CreateUnitTest(field_lower_bound_check "field_lower_bound_check_tests.cpp"
    LIBS eamxx_property_checks)

  # Test FieldUpperBoundCheck
  CreateUnitTest(field_upper_bound_check "field_upper_bound_check_tests.cpp"
    LIBS eamxx_property_checks)
endif()
