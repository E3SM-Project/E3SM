include(ScreamUtils)

set(NEED_LIBS p3 physics_share scream_share)
set(P3_TESTS_SRCS
  p3_tests.cpp
  p3_unit_tests.cpp
  p3_ice_tables_unit_tests.cpp
  p3_table3_unit_tests.cpp
  p3_back_to_cell_average_unit_tests.cpp
  p3_prevent_ice_overdepletion_unit_tests.cpp
  p3_find_unit_tests.cpp
  p3_upwind_unit_tests.cpp
  p3_calc_rime_density_unit_tests.cpp
  p3_cldliq_imm_freezing_unit_tests.cpp
  p3_rain_imm_freezing_unit_tests.cpp
  p3_droplet_self_coll_unit_tests.cpp
  p3_cloud_sed_unit_tests.cpp
  p3_cloud_rain_acc_unit_tests.cpp
  p3_ice_sed_unit_tests.cpp
  p3_ice_collection_unit_tests.cpp
  p3_rain_sed_unit_tests.cpp
  p3_dsd2_unit_tests.cpp
  p3_rain_self_collection_tests.cpp
  p3_autoconversion_unit_tests.cpp
  p3_ice_relaxation_timescale_unit_tests.cpp
  p3_calc_liq_relaxation_timescale_unit_tests.cpp
  p3_ice_nucleation_unit_tests.cpp
  p3_ice_melting_unit_tests.cpp
  p3_ice_cldliq_wet_growth_unit_tests.cpp
  p3_get_latent_heat_unit_tests.cpp
  p3_subgrid_variance_scaling_unit_tests.cpp
  p3_check_values_unit_tests.cpp
  p3_incloud_mixingratios_unit_tests.cpp
  p3_main_unit_tests.cpp
)

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT ${SCREAM_BASELINES_ONLY})
  CreateUnitTest(p3_tests "${P3_TESTS_SRCS}" "${NEED_LIBS}" THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC} DEP p3_tests_ut_np1_omp1)
endif()

CreateUnitTest(p3_run_and_cmp_cxx "p3_run_and_cmp.cpp" "${NEED_LIBS}" THREADS ${SCREAM_TEST_MAX_THREADS} EXE_ARGS "${SCREAM_TEST_DATA_DIR}/p3_run_and_cmp.baseline" DEP p3_tests_ut_np1_omp1 EXCLUDE_MAIN_CPP)

CreateUnitTest(p3_run_and_cmp_f90 "p3_run_and_cmp.cpp" "${NEED_LIBS}" THREADS ${SCREAM_TEST_MAX_THREADS} EXE_ARGS "-f ${SCREAM_TEST_DATA_DIR}/p3_run_and_cmp.baseline" DEP p3_tests_ut_np1_omp1 EXCLUDE_MAIN_CPP)

# By default, baselines should be created using all fortran
add_custom_target(p3_baseline_f90
  COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=${SCREAM_TEST_MAX_THREADS} $<TARGET_FILE:p3_run_and_cmp_f90> -f -g ${SCREAM_TEST_DATA_DIR}/p3_run_and_cmp.baseline)

add_custom_target(p3_baseline_cxx
  COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=${SCREAM_TEST_MAX_THREADS} $<TARGET_FILE:p3_run_and_cmp_cxx> -g ${SCREAM_TEST_DATA_DIR}/p3_run_and_cmp.baseline)

add_dependencies(baseline p3_baseline_f90)

add_dependencies(baseline_cxx p3_baseline_cxx)

configure_file(${SCREAM_DATA_DIR}/p3_lookup_table_1.dat-v4 data/p3_lookup_table_1.dat-v4 COPYONLY)
