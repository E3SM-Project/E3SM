include(ScreamUtils)

add_subdirectory(infra)

CreateUnitTest(tms_tests
  SOURCES compute_tms_tests.cpp
  LIBS tms tms_test_infra
  EXE_ARGS "--args ${SCREAM_BASELINE_FILE_ARG}"
  THREADS ${SCREAM_BASELINE_TEST_THREADS}
  LABELS tms physics baseline_gen baseline_cmp
)
