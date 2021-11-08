include(ScreamUtils)

set(NEED_LIBS ekat scream_share)

# Add the fake test
CreateUnitTest(fake_test "fake_test.cpp" "${NEED_LIBS}")

if (SCREAM_TEST_MAX_TOTAL_THREADS GREATER_EQUAL 16)
  CreateUnitTest(resource_spread_thread "resource_spread.cpp" "${NEED_LIBS}"
    PRINT_OMP_AFFINITY
    THREADS 1 4 1
    MPI_RANKS 4)

  CreateUnitTest(resource_spread_rank "" ""
    TEST_EXE resource_spread_thread
    PRINT_OMP_AFFINITY
    THREADS 4
    MPI_RANKS 1 4 1)
endif()
