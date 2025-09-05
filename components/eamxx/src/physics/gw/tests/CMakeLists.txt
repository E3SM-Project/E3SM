include(ScreamUtils)

add_subdirectory(infra)

set(GW_TESTS_SRCS
    gw_gwd_compute_tendencies_from_stress_divergence_tests.cpp
    gw_gw_prof_tests.cpp
    gw_momentum_energy_conservation_tests.cpp
    gw_gwd_compute_stress_profiles_and_diffusivities_tests.cpp
    gw_gwd_project_tau_tests.cpp
    gw_gwd_precalc_rhoi_tests.cpp
    gw_gw_drag_prof_tests.cpp
    gw_gw_front_project_winds_tests.cpp
    gw_gw_front_gw_sources_tests.cpp
    gw_gw_cm_src_tests.cpp
    gw_gw_convect_project_winds_tests.cpp
    gw_gw_heating_depth_tests.cpp
    gw_gw_storm_speed_tests.cpp
    gw_gw_convect_gw_sources_tests.cpp
    gw_gw_beres_src_tests.cpp
    gw_gw_ediff_tests.cpp
    gw_gw_diff_tend_tests.cpp
    gw_gw_oro_src_tests.cpp
    gw_vd_lu_decomp_tests.cpp
    gw_vd_lu_solve_tests.cpp
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
