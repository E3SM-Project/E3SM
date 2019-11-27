include(ScreamUtils)

set(NEED_LIBS p3 scream_share)
set(P3_TESTS_SRCS
  p3_tests.cpp
  p3_unit_tests.cpp
  p3_ice_tables_unit_tests.cpp
  p3_table3_unit_tests.cpp
  p3_find_unit_tests.cpp
  p3_upwind_unit_tests.cpp
  p3_cloud_sed_unit_tests.cpp
  p3_ice_sed_unit_tests.cpp
  p3_dsd2_unit_tests.cpp
)
CreateUnitTest(p3_tests "${P3_TESTS_SRCS}" "${NEED_LIBS}" THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC} DEP p3_tests_ut_np1_omp1)

CreateUnitTest(p3_run_and_cmp "p3_run_and_cmp.cpp" "${NEED_LIBS}" THREADS ${SCREAM_TEST_MAX_THREADS} EXE_ARGS "${SCREAM_TEST_DATA_DIR}/p3_run_and_cmp.baseline" DEP p3_tests_ut_np1_omp1 OPTIONAL EXCLUDE_MAIN_CPP)

# By default, baselines should be created using all fortran
add_custom_target(p3_baseline
  COMMAND $<TARGET_FILE:p3_run_and_cmp> -f -g ${SCREAM_TEST_DATA_DIR}/p3_run_and_cmp.baseline)

add_dependencies(baseline p3_baseline)

configure_file(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4 p3_lookup_table_1.dat-v4 COPYONLY)


