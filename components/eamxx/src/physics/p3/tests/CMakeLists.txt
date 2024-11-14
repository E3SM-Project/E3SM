include(ScreamUtils)

add_subdirectory(infra)

set(P3_TESTS_SRCS
    p3_tests.cpp
    p3_unit_tests.cpp
    p3_ice_tables_unit_tests.cpp
    p3_table3_unit_tests.cpp
    p3_back_to_cell_average_unit_tests.cpp
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
    p3_evaporate_rain_unit_tests.cpp
    p3_ice_cldliq_wet_growth_unit_tests.cpp
    p3_subgrid_variance_scaling_unit_tests.cpp
    p3_check_values_unit_tests.cpp
    p3_incloud_mixingratios_unit_tests.cpp
    p3_main_unit_tests.cpp
    p3_ice_supersat_conservation_tests.cpp
    p3_nc_conservation_tests.cpp
    p3_nr_conservation_tests.cpp
    p3_ni_conservation_tests.cpp
    p3_ice_deposition_sublimation_tests.cpp
    p3_prevent_liq_supersaturation_tests.cpp
    ) # P3_TESTS_SRCS

if (SCREAM_DEBUG AND NOT SCREAM_TEST_SIZE STREQUAL "SHORT")
  set (FORCE_RUN_DIFF_FAILS "WILL_FAIL")
else ()
  set (FORCE_RUN_DIFF_FAILS "")
endif()

# All tests should understand the same baseline args
if (SCREAM_ENABLE_BASELINE_TESTS)
  if (SCREAM_ONLY_GENERATE_BASELINES)
    set(BASELINE_FILE_ARG "-g -b ${SCREAM_BASELINES_DIR}/data")
    # We don't want to do thread spreads when generating. That
    # could cause race conditions in the file system.
    set(P3_THREADS "${SCREAM_TEST_MAX_THREADS}")
  else()
    set(BASELINE_FILE_ARG "-c -b ${SCREAM_BASELINES_DIR}/data")
    set(P3_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
  endif()
else()
  set(BASELINE_FILE_ARG "-n") # no baselines
  set(P3_THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC})
endif()

CreateUnitTest(p3_tests "${P3_TESTS_SRCS}"
  LIBS p3 p3_test_infra
  EXE_ARGS "--args ${BASELINE_FILE_ARG}"
  THREADS ${P3_THREADS}
  LABELS "p3;physics;baseline_gen;baseline_cmp")

# Make sure that a diff in the two implementation triggers a failed test (in debug only)
# No need to run lots of different thread counts.
if (SCREAM_ENABLE_BASELINE_TESTS)
  CreateUnitTest (p3_tests_fail p3_rain_sed_unit_tests.cpp
    LIBS p3 p3_test_infra
    EXE_ARGS "--args ${BASELINE_FILE_ARG}"
    COMPILER_CXX_DEFS SCREAM_FORCE_RUN_DIFF
    LABELS "p3;physics;fail"
    ${FORCE_RUN_DIFF_FAILS})
endif()

# If small kernels are ON, we don't need a separate executable to test them.
# Also, we never want to generate baselines with this separate executable
if (NOT SCREAM_P3_SMALL_KERNELS AND NOT SCREAM_ONLY_GENERATE_BASELINES)
  # Note: Only the p3_main test does something different when
  # small kernels are on. The SK dispatch routines are mostly trivial
  # and it's not worth adding tons of test infrastructure to support
  # BFB unit tests for these.
  CreateUnitTest(p3_sk_tests "p3_main_unit_tests.cpp"
    LIBS p3_sk p3_test_infra
    EXE_ARGS "--args ${BASELINE_FILE_ARG}"
    THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
    LABELS "p3_sk;physics;baseline_cmp")
endif()

# Note: the baseline_gen label label is really only used if SCREAM_ONLY_GENERATE_BASELINES=ON, but no harm adding it
CreateUnitTest(p3_run_and_cmp "p3_run_and_cmp.cpp"
  LIBS p3 p3_test_infra
  EXCLUDE_MAIN_CPP
  THREADS ${SCREAM_TEST_MAX_THREADS}
  EXE_ARGS "${BASELINE_FILE_ARG}"
  LABELS "p3;physics;baseline_gen;baseline_cmp")

# Make sure that a diff from baselines triggers a failed test (in debug only)
if (SCREAM_ENABLE_BASELINE_TESTS)
  CreateUnitTest(p3_run_and_cmp_fail "p3_run_and_cmp.cpp"
    LIBS p3 p3_test_infra
    COMPILER_CXX_DEFS SCREAM_FORCE_RUN_DIFF
    THREADS ${SCREAM_TEST_MAX_THREADS}
    EXE_ARGS "${BASELINE_FILE_ARG}"
    LABELS "p3;physics;fail"
    EXCLUDE_MAIN_CPP
    ${FORCE_RUN_DIFF_FAILS})
endif()
