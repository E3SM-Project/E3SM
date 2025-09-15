# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_ONLY_GENERATE_BASELINES)
  include(ScreamUtils)

  # Test TimeStamp
  CreateUnitTest(time_stamp "time_stamp_tests.cpp" LIBS eamxx_utils)

  # Test bfb hash
  CreateUnitTest(bfb_hash "bfbhash_tests.cpp" LIBS eamxx_utils)

  # Test array utils
  CreateUnitTest(array_utils "array_utils_tests.cpp" LIBS eamxx_utils)

  # Test combine operations
  CreateUnitTest(combine_ops "combine_ops_tests.cpp" LIBS eamxx_utils)

  # Test miscellanea utils
  CreateUnitTest(misc_utils "misc_utils_tests.cpp" LIBS eamxx_utils)
endif()
