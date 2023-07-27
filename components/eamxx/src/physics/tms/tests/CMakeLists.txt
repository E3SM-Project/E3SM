INCLUDE (ScreamUtils)

SET(NEED_LIBS tms physics_share scream_share)
set(TMS_TESTS_SRCS
    compute_tms_tests.cpp
) # TMS_TESTS_SRCS

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  CreateUnitTest(tms_tests "${TMS_TESTS_SRCS}" "${NEED_LIBS}" THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC} DEP tms_tests_ut_np1_omp1)
endif()
