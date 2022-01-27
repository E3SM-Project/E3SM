include(ScreamUtils)

# NOTE: we don't need any scream libs for these test, but scream's CreateUnitTest
# has libs as a required arg. So use the raw Ekat version (at least to create
# the test executable).

# Test to ensure that a build failure is detected by our testing scripts
EkatCreateUnitTest(fake_test fake_test.cpp)

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
