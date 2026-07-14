include(ScreamUtils)

add_subdirectory(infra)

set(GW_TESTS_SRCS
    gw_compute_tendencies_from_stress_divergence_tests.cpp
    gw_prof_tests.cpp
    gw_momentum_energy_conservation_tests.cpp
    gw_compute_stress_profiles_and_diffusivities_tests.cpp
    gw_project_tau_tests.cpp
    gw_precalc_rhoi_tests.cpp
    gw_drag_prof_tests.cpp
    gw_front_project_winds_tests.cpp
    gw_front_gw_sources_tests.cpp
    gw_cm_src_tests.cpp
    gw_convect_project_winds_tests.cpp
    gw_heating_depth_tests.cpp
    gw_storm_speed_tests.cpp
    gw_convect_gw_sources_tests.cpp
    gw_beres_src_tests.cpp
    gw_ediff_tests.cpp
    gw_diff_tend_tests.cpp
    gw_oro_src_tests.cpp
    gw_vd_lu_decomp_tests.cpp
    gw_vd_lu_solve_tests.cpp
) # GW_TESTS_SRCS

CreateUnitTest(gw_tests
  SOURCES ${GW_TESTS_SRCS}
  LIBS gw gw_test_infra
  EXE_ARGS "--args ${SCREAM_BASELINE_FILE_ARG}"
  THREADS ${SCREAM_BASELINE_TEST_THREADS}
  LABELS gw physics baseline_gen baseline_cmp
)
