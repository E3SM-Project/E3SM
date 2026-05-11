# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test FieldNaNCheck
  CreateUnitTest(field_nan_check
    SOURCES field_nan_check_tests.cpp
    LIBS eamxx_property_checks
  )

  # Test FieldWithinIntervalCheck
  CreateUnitTest(field_within_interval_check
    SOURCES field_within_interval_check_tests.cpp
    LIBS eamxx_property_checks
  )

  # Test FieldLowerBoundCheck
  CreateUnitTest(field_lower_bound_check
    SOURCES field_lower_bound_check_tests.cpp
    LIBS eamxx_property_checks
  )

  # Test FieldUpperBoundCheck
  CreateUnitTest(field_upper_bound_check
    SOURCES field_upper_bound_check_tests.cpp
    LIBS eamxx_property_checks
  )
endif()
