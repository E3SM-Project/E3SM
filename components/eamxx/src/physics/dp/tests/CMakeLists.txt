INCLUDE (ScreamUtils)

SET (NEED_LIBS dp physics_share scream_share)
set(DP_TESTS_SRCS
  dp_unit_tests.cpp
    ) # DP_TESTS_SRCS

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  CreateUnitTest(dp_tests    "${DP_TESTS_SRCS}" "${NEED_LIBS}"    THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC} DEP dp_tests_ut_np1_omp1)
endif()
