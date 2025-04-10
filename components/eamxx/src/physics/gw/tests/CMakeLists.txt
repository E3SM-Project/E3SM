include(ScreamUtils)

add_subdirectory(infra)

set(GW_TESTS_SRCS
  
) # GW_TESTS_SRCS

# All tests should understand the same baseline args
if (SCREAM_ENABLE_BASELINE_TESTS)
  if (SCREAM_ONLY_GENERATE_BASELINES)
    set(BASELINE_FILE_ARG "-g -b ${SCREAM_BASELINES_DIR}/data")
    # We don't want to do thread spreads when generating. That
    # could cause race conditions in the file system.
    set(GW_THREADS "${SCREAM_TEST_MAX_THREADS}")
  else()
    set(BASELINE_FILE_ARG "-c -b ${SCREAM_BASELINES_DIR}/data")
    set(GW_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
  endif()
else()
  set(BASELINE_FILE_ARG "-n") # no baselines
  set(GW_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
endif()

CreateUnitTest(gw_tests "${GW_TESTS_SRCS}"
  LIBS gw gw_test_infra
  EXE_ARGS "--args ${BASELINE_FILE_ARG}"
  THREADS ${GW_THREADS}
  LABELS "gw;physics;baseline_gen;baseline_cmp")
