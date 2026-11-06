include(ScreamUtils)

add_subdirectory(infra)

CreateUnitTest(rayleigh_friction_tests
  SOURCES rayleigh_friction_init_tests.cpp
          rayleigh_friction_tend_tests.cpp
  LIBS rayleigh_friction rayleigh_friction_test_infra
  EXE_ARGS "--args ${SCREAM_BASELINE_FILE_ARG}"
  THREADS ${SCREAM_BASELINE_TEST_THREADS}
  LABELS rayleigh_friction physics baseline_gen baseline_cmp
)
