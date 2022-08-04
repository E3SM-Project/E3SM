include(ScreamUtils)

# These tests are "fake" tests, in the sense that they are placeholder
# that guarantee a certain behavior from the ctest executable.
# We use them to verify that the enclosing testing scripts can correctly
# capture a failure from our tests, as well as make ctest spread the
# tests correctly across available resources.

# NOTE: we don't need any scream libs for these test, but scream's CreateUnitTest
# has libs as a required arg. So use the raw Ekat version (at least to create
# the test executable).

# Test to ensure that a build failure is detected by our testing scripts
EkatCreateUnitTest(build_fail build_fail.cpp
  LABELS "fail")

# Test to ensure that a test failure is detected by our testing scripts
EkatCreateUnitTest(test_fail test_fail.cpp
  LABELS "fail")

# Tests to ensure the testing infrastructure is correctly spreading
# concurrent tests across the available resources
if (SCREAM_TEST_MAX_TOTAL_THREADS GREATER_EQUAL 16)
  EkatCreateUnitTestExec(resource_spread resource_spread.cpp)

  CreateUnitTestFromExec(resource_spread_thread resource_spread
    PRINT_OMP_AFFINITY
    THREADS 1 4 1
    MPI_RANKS 4)

  CreateUnitTestFromExec(resource_spread_rank resource_spread
    PRINT_OMP_AFFINITY
    THREADS 4
    MPI_RANKS 1 4 1)
endif()

add_subdirectory(test_level_check)
