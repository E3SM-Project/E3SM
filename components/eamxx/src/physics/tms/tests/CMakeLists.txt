INCLUDE (ScreamUtils)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  CreateUnitTest(tms_tests compute_tms_tests.cpp
    LIBS tms
    THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
  )
endif()
