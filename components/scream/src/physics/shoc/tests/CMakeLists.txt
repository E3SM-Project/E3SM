INCLUDE (ScreamUtils)

SET (NEED_LIBS shoc physics_common scream_control scream_share)
set(SHOC_TESTS_SRCS
    shoc_tests.cpp
    shoc_grid_tests.cpp
    shoc_unit_tests.cpp)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  CreateUnitTest(shoc_tests "${SHOC_TESTS_SRCS}" "${NEED_LIBS}" THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC} DEP shoc_tests_ut_np1_omp1)
  # CreateUnitTest(shoc_stand_alone "shoc_stand_alone.cpp" "${NEED_LIBS}")
endif()

CreateUnitTest(shoc_run_and_cmp "shoc_run_and_cmp.cpp" "${NEED_LIBS}" THREADS ${SCREAM_TEST_MAX_THREADS} EXE_ARGS "${SCREAM_TEST_DATA_DIR}/shoc_run_and_cmp.baseline" DEP shoc_tests_ut_np1_omp1 OPTIONAL EXCLUDE_MAIN_CPP)

# By default, baselines should be created using all fortran
add_custom_target(shoc_baseline
  COMMAND $<TARGET_FILE:shoc_run_and_cmp> -f -g ${SCREAM_TEST_DATA_DIR}/shoc_run_and_cmp.baseline)

add_dependencies(baseline shoc_baseline)

