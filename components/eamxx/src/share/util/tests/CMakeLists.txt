# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test column ops
  CreateUnitTest(column_ops
    SOURCES column_ops.cpp
    LIBS eamxx_utils
  )

  # Test TimeStamp
  CreateUnitTest(time_stamp
    SOURCES time_stamp_tests.cpp
    LIBS eamxx_utils
  )

  # Test bfb hash
  CreateUnitTest(bfb_hash
    SOURCES bfbhash_tests.cpp
    LIBS eamxx_utils
  )

  # Test array utils
  CreateUnitTest(array_utils
    SOURCES array_utils_tests.cpp
    LIBS eamxx_utils
  )

  # Test combine operations
  CreateUnitTest(combine_ops
    SOURCES combine_ops_tests.cpp
    LIBS eamxx_utils
  )

  # Test miscellanea utils
  CreateUnitTest(misc_utils
    SOURCES misc_utils_tests.cpp
    LIBS eamxx_utils
  )
endif()
