include(ScreamUtils)

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
    p3_get_latent_heat_unit_tests.cpp
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

# NOTE: tests inside this if statement won't be built in a baselines-only build
if (NOT SCREAM_BASELINES_ONLY)
  CreateUnitTest(p3_tests "${P3_TESTS_SRCS}"
      LIBS p3
      THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
      LABELS "p3;physics")

  # Make sure that a diff in the two implementation triggers a failed test (in debug only)
  CreateUnitTest (p3_tests_fail p3_rain_sed_unit_tests.cpp
      LIBS p3
      COMPILER_CXX_DEFS SCREAM_FORCE_RUN_DIFF
      THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
      LABELS "p3;physics;fail"
      ${FORCE_RUN_DIFF_FAILS})

  if (NOT SCREAM_SMALL_KERNELS)
     CreateUnitTest(p3_sk_tests "${P3_TESTS_SRCS}"
        LIBS p3_sk
        THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
        LABELS "p3_sk;physics")

  # Make sure that a diff in the two implementation triggers a failed test (in debug only)
     CreateUnitTest (p3_sk_tests_fail p3_rain_sed_unit_tests.cpp
        LIBS p3_sk
        COMPILER_CXX_DEFS SCREAM_FORCE_RUN_DIFF
        THREADS 1 ${SCREAM_TEST_MAX_THREADS} ${SCREAM_TEST_THREAD_INC}
        LABELS "p3_sk;physics;fail"
        ${FORCE_RUN_DIFF_FAILS})
  endif()
endif()

if (SCREAM_ENABLE_BASELINE_TESTS)
  set(BASELINE_FILE_ARG "-b ${SCREAM_TEST_DATA_DIR}/p3_run_and_cmp.baseline")

  CreateUnitTestExec(p3_run_and_cmp "p3_run_and_cmp.cpp"
      LIBS p3
      EXCLUDE_MAIN_CPP)

  CreateUnitTestFromExec(p3_run_and_cmp_cxx p3_run_and_cmp
      THREADS ${SCREAM_TEST_MAX_THREADS}
      EXE_ARGS "${BASELINE_FILE_ARG}"
      LABELS "p3;physics")

  CreateUnitTestFromExec(p3_run_and_cmp_f90 p3_run_and_cmp
      THREADS ${SCREAM_TEST_MAX_THREADS}
      EXE_ARGS "-f ${BASELINE_FILE_ARG}"
      LABELS "p3;physics")

  # Make sure that a diff from baselines triggers a failed test (in debug only)
  CreateUnitTest(p3_run_and_cmp_cxx_fail "p3_run_and_cmp.cpp"
      LIBS p3
      COMPILER_CXX_DEFS SCREAM_FORCE_RUN_DIFF
      THREADS ${SCREAM_TEST_MAX_THREADS}
      EXE_ARGS "${BASELINE_FILE_ARG}"
      LABELS "p3;physics;fail"
      EXCLUDE_MAIN_CPP
      ${FORCE_RUN_DIFF_FAILS})

  #
  # Use fake tests to generate shell commands to generate baselines
  #
  CreateUnitTestFromExec(p3_baseline_f90_fake p3_run_and_cmp
    THREADS ${SCREAM_TEST_MAX_THREADS}
    EXE_ARGS "-f -g ${BASELINE_FILE_ARG}"
    PROPERTIES DISABLED True)

  CreateUnitTestFromExec(p3_baseline_cxx_fake p3_run_and_cmp
    THREADS ${SCREAM_TEST_MAX_THREADS}
    EXE_ARGS "-g ${BASELINE_FILE_ARG}"
    PROPERTIES DISABLED True)

  if (SCREAM_TEST_MAX_THREADS GREATER 1)
    get_test_property(p3_baseline_f90_fake_omp${SCREAM_TEST_MAX_THREADS} FULL_TEST_COMMAND P3_F90_GEN)
    get_test_property(p3_baseline_cxx_fake_omp${SCREAM_TEST_MAX_THREADS} FULL_TEST_COMMAND P3_CXX_GEN)
  else()
    get_test_property(p3_baseline_f90_fake FULL_TEST_COMMAND P3_F90_GEN)
    get_test_property(p3_baseline_cxx_fake FULL_TEST_COMMAND P3_CXX_GEN)
  endif()

  if (P3_F90_GEN STREQUAL "NOTFOUND")
    message(FATAL_ERROR "Could not get FULL_TEST_COMMAND for p3_baseline fake test")
  endif()

  separate_arguments(P3_F90_GEN_ARGS UNIX_COMMAND "${P3_F90_GEN}")
  separate_arguments(P3_CXX_GEN_ARGS UNIX_COMMAND "${P3_CXX_GEN}")

  add_custom_target(p3_baseline_f90
    COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=${SCREAM_TEST_MAX_THREADS} ${P3_F90_GEN_ARGS})

  add_custom_target(p3_baseline_cxx
    COMMAND ${CMAKE_COMMAND} -E env OMP_NUM_THREADS=${SCREAM_TEST_MAX_THREADS} ${P3_CXX_GEN_ARGS})

  # By default, baselines should be created using all fortran (make baseline). If the user wants
  # to use CXX to generate their baselines, they should use "make baseline_cxx".
  add_dependencies(baseline     p3_baseline_f90)
  add_dependencies(baseline_cxx p3_baseline_cxx)
endif()
