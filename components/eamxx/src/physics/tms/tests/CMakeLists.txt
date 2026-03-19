include(ScreamUtils)

add_subdirectory(infra)

set(TMS_TESTS_SRCS
  compute_tms_tests.cpp
) # TMS_TESTS_SRCS

# All tests should understand the same baseline args
if (SCREAM_ENABLE_BASELINE_TESTS)
  if (SCREAM_ONLY_GENERATE_BASELINES)
    set(BASELINE_FILE_ARG "-g -b ${SCREAM_BASELINES_DIR}/data")
    # We don't want to do thread spreads when generating. That
    # could cause race conditions in the file system.
    set(TMS_THREADS "${SCREAM_TEST_MAX_THREADS}")
  else()
    set(BASELINE_FILE_ARG "-c -b ${SCREAM_BASELINES_DIR}/data")
    set(TMS_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
  endif()
else()
  set(BASELINE_FILE_ARG "-n") # no baselines
  set(TMS_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
endif()

CreateUnitTest(tms_tests "${TMS_TESTS_SRCS}"
  LIBS tms tms_test_infra
  EXE_ARGS "--args ${BASELINE_FILE_ARG}"
  THREADS ${TMS_THREADS}
  LABELS "tms;physics;baseline_gen;baseline_cmp")
