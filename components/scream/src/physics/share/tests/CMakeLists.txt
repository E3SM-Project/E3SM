include(ScreamUtils)

set(NEED_LIBS physics_share scream_share)
set(PHYSICS_TESTS_SRCS
  physics_saturation_unit_tests.cpp
  physics_test_data_unit_tests.cpp
)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  CreateUnitTest(physics_tests "${PHYSICS_TESTS_SRCS}" "${NEED_LIBS}" THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC} DEP physics_tests_ut_np1_omp1)
endif()

