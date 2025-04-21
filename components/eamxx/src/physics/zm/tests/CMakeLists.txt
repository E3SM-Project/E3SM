include(ScreamUtils)

add_subdirectory(infra)

set(ZM_TESTS_SRCS
  zm_test_find_mse_max.cpp
) # ZM_TESTS_SRCS

# All tests should understand the same baseline args
if (SCREAM_ENABLE_BASELINE_TESTS)
  if (SCREAM_ONLY_GENERATE_BASELINES)
    set(BASELINE_FILE_ARG "-g -b ${SCREAM_BASELINES_DIR}/data")
    # We don't want to do thread spreads when generating. That
    # could cause race conditions in the file system.
    set(ZM_THREADS "${SCREAM_TEST_MAX_THREADS}")
  else()
    set(BASELINE_FILE_ARG "-c -b ${SCREAM_BASELINES_DIR}/data")
    set(GZM_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
  endif()
else()
  set(BASELINE_FILE_ARG "-n") # no baselines
  set(GZM_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
endif()

CreateUnitTest(zm_tests "${ZM_TESTS_SRCS}"
  LIBS zm zm_test_infra
  EXE_ARGS "--args ${BASELINE_FILE_ARG}"
  THREADS ${ZM_THREADS}
  LABELS "zm;physics;baseline_gen;baseline_cmp")
