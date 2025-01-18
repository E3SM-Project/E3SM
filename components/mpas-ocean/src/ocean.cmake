
# build_options.mk stuff handled here
list(APPEND CPPDEFS "-DCORE_OCEAN")
list(APPEND CPPDEFS "-DEXCLUDE_INIT_MODE")
list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_ocean/shared") # Only need this for '#include "../inc/core_variables.inc"' to work
list(APPEND INCLUDES "core_ocean/gotm/include")

# check if lapack is linked
if (LAPACK_FOUND AND BLAS_FOUND)
  list(APPEND CPPDEFS "-DUSE_LAPACK")
  list(APPEND LIBRARIES BLAS::BLAS LAPACK::LAPACK)
endif()

if (USE_PETSC)
  list(APPEND CPPDEFS "-DUSE_PETSC")
  list(APPEND INCLUDES "${PETSC_INCLUDES}")
  list(APPEND LIBRARIES "${PETSC_LIBRARIES}")
endif()

# driver (files live in E3SM)
list(APPEND RAW_SOURCES
  ../../mpas-ocean/driver/ocn_comp_mct.F
  ../../mpas-ocean/driver/mpaso_cpl_indices.F
  ../../mpas-ocean/driver/mpaso_mct_vars.F
)

# dycore
list(APPEND RAW_SOURCES
  core_ocean/mode_forward/mpas_ocn_forward_mode.F
  core_ocean/mode_forward/mpas_ocn_time_integration.F
  core_ocean/mode_forward/mpas_ocn_time_integration_rk4.F
  core_ocean/mode_forward/mpas_ocn_time_integration_split.F
  core_ocean/mode_forward/mpas_ocn_time_integration_si.F
  core_ocean/mode_forward/mpas_ocn_time_integration_lts.F
  core_ocean/mode_forward/mpas_ocn_time_integration_fblts.F
  core_ocean/mode_forward/mpas_ocn_time_integration_split_ab2.F

  core_ocean/mode_analysis/mpas_ocn_analysis_mode.F

  core_ocean/shared/mpas_ocn_init_routines.F
  core_ocean/shared/mpas_ocn_gm.F
  core_ocean/shared/mpas_ocn_submesoscale_eddies.F
  core_ocean/shared/mpas_ocn_eddy_parameterization_helpers.F
  core_ocean/shared/mpas_ocn_diagnostics.F
  core_ocean/shared/mpas_ocn_diagnostics_variables.F
  core_ocean/shared/mpas_ocn_mesh.F
  core_ocean/shared/mpas_ocn_thick_ale.F
  core_ocean/shared/mpas_ocn_equation_of_state.F
  core_ocean/shared/mpas_ocn_equation_of_state_jm.F
  core_ocean/shared/mpas_ocn_equation_of_state_linear.F
  core_ocean/shared/mpas_ocn_equation_of_state_wright.F
  core_ocean/shared/mpas_ocn_thick_hadv.F
  core_ocean/shared/mpas_ocn_thick_vadv.F
  core_ocean/shared/mpas_ocn_thick_surface_flux.F
  core_ocean/shared/mpas_ocn_vel_hadv_coriolis.F
  core_ocean/shared/mpas_ocn_vel_vadv.F
  core_ocean/shared/mpas_ocn_vel_hmix.F
  core_ocean/shared/mpas_ocn_vel_hmix_del2.F
  core_ocean/shared/mpas_ocn_vel_hmix_leith.F
  core_ocean/shared/mpas_ocn_vel_hmix_del4.F
  core_ocean/shared/mpas_ocn_vel_forcing.F
  core_ocean/shared/mpas_ocn_vel_forcing_surface_stress.F
  core_ocean/shared/mpas_ocn_vel_forcing_explicit_bottom_drag.F
  core_ocean/shared/mpas_ocn_vel_pressure_grad.F
  core_ocean/shared/mpas_ocn_vel_forcing_topographic_wave_drag.F
  core_ocean/shared/mpas_ocn_vertical_advection.F
  core_ocean/shared/mpas_ocn_vertical_regrid.F
  core_ocean/shared/mpas_ocn_vertical_remap.F
  core_ocean/shared/mpas_ocn_vmix.F
  core_ocean/shared/mpas_ocn_vmix_coefs_redi.F
  core_ocean/shared/mpas_ocn_vmix_cvmix.F
  core_ocean/shared/mpas_ocn_vmix_gotm.F
  core_ocean/shared/mpas_ocn_vel_self_attraction_loading.F
  core_ocean/shared/mpas_ocn_tendency.F
  core_ocean/shared/mpas_ocn_tracer_hmix.F
  core_ocean/shared/mpas_ocn_tracer_hmix_del2.F
  core_ocean/shared/mpas_ocn_tracer_hmix_del4.F
  core_ocean/shared/mpas_ocn_tracer_hmix_redi.F
  core_ocean/shared/mpas_ocn_tracer_advection.F
  core_ocean/shared/mpas_ocn_tracer_advection_mono.F
  core_ocean/shared/mpas_ocn_tracer_advection_std.F
  core_ocean/shared/mpas_ocn_tracer_advection_vert.F
  core_ocean/shared/mpas_ocn_tracer_advection_shared.F
  core_ocean/shared/mpas_ocn_tracer_nonlocalflux.F
  core_ocean/shared/mpas_ocn_tracer_short_wave_absorption.F
  core_ocean/shared/mpas_ocn_tracer_short_wave_absorption_jerlov.F
  core_ocean/shared/mpas_ocn_tracer_short_wave_absorption_variable.F
  core_ocean/shared/mpas_ocn_tracer_surface_restoring.F
  core_ocean/shared/mpas_ocn_tracer_interior_restoring.F
  core_ocean/shared/mpas_ocn_tracer_exponential_decay.F
  core_ocean/shared/mpas_ocn_tracer_ideal_age.F
  core_ocean/shared/mpas_ocn_tracer_CFC.F
  core_ocean/shared/mpas_ocn_tracer_TTD.F
  core_ocean/shared/mpas_ocn_tracer_ecosys.F
  core_ocean/shared/mpas_ocn_tracer_DMS.F
  core_ocean/shared/mpas_ocn_tracer_MacroMolecules.F
  core_ocean/shared/mpas_ocn_transport_tests.F
  core_ocean/shared/mpas_ocn_high_freq_thickness_hmix_del2.F
  core_ocean/shared/mpas_ocn_tracer_surface_flux_to_tend.F
  core_ocean/shared/mpas_ocn_test.F
  core_ocean/shared/mpas_ocn_constants.F
  core_ocean/shared/mpas_ocn_config.F
  core_ocean/shared/mpas_ocn_forcing.F
  core_ocean/shared/mpas_ocn_surface_bulk_forcing.F
  core_ocean/shared/mpas_ocn_surface_land_ice_fluxes.F
  core_ocean/shared/mpas_ocn_effective_density_in_land_ice.F
  core_ocean/shared/mpas_ocn_frazil_forcing.F
  core_ocean/shared/mpas_ocn_tidal_forcing.F
  core_ocean/shared/mpas_ocn_time_average_coupled.F
  core_ocean/shared/mpas_ocn_framework_forcing.F
  core_ocean/shared/mpas_ocn_time_varying_forcing.F
  core_ocean/shared/mpas_ocn_wetting_drying.F
  core_ocean/shared/mpas_ocn_vel_tidal_potential.F
  core_ocean/shared/mpas_ocn_stokes_drift.F
  core_ocean/shared/mpas_ocn_manufactured_solution.F
  core_ocean/shared/mpas_ocn_subgrid.F
  core_ocean/shared/mpas_ocn_scaled_dismf.F
)

set(OCEAN_DRIVER
  core_ocean/driver/mpas_ocn_core.F
  core_ocean/driver/mpas_ocn_core_interface.F
)
list(APPEND RAW_SOURCES ${OCEAN_DRIVER})
list(APPEND DISABLE_QSMP ${OCEAN_DRIVER})

# Add CVMix
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/core_ocean/cvmix/.git)
  message(FATAL_ERROR "Missing core_ocean/cvmix/.git, did you forget to 'git submodule update --init --recursive' ?")
endif()
set(CVMIX_FILES
  core_ocean/cvmix/src/shared/cvmix_kinds_and_types.F90
  core_ocean/cvmix/src/shared/cvmix_background.F90
  core_ocean/cvmix/src/shared/cvmix_convection.F90
  core_ocean/cvmix/src/shared/cvmix_ddiff.F90
  core_ocean/cvmix/src/shared/cvmix_kpp.F90
  core_ocean/cvmix/src/shared/cvmix_math.F90
  core_ocean/cvmix/src/shared/cvmix_put_get.F90
  core_ocean/cvmix/src/shared/cvmix_shear.F90
  core_ocean/cvmix/src/shared/cvmix_tidal.F90
  core_ocean/cvmix/src/shared/cvmix_utils.F90
)

# Add BGC
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/core_ocean/BGC/.git)
  message(FATAL_ERROR "Missing core_ocean/BGC/.git, did you forget to 'git submodule update --init --recursive' ?")
endif()
set(BGC_FILES
  core_ocean/BGC/DMS_mod.F90
  core_ocean/BGC/DMS_parms.F90
  core_ocean/BGC/MACROS_mod.F90
  core_ocean/BGC/MACROS_parms.F90
)

# Add MARBL
if (NOT EXISTS core_ocean/MARBL/.git)
  message(FATAL "Missing core_ocean/MARBL/.git, did you forget to 'git submodule update --init --recursive' ?")
endif()
set(MARBL_FILES
  core_ocean/MARBL/src/marbl_ciso_diagnostics_mod.F90
  core_ocean/MARBL/src/marbl_ciso_init_mod.F90
  core_ocean/MARBL/src/marbl_ciso_interior_tendency_mod.F90
  core_ocean/MARBL/src/marbl_ciso_surface_flux_mod.F90
  core_ocean/MARBL/src/marbl_co2calc_mod.F90
  core_ocean/MARBL/src/marbl_constants_mod.F90
  core_ocean/MARBL/src/marbl_debug_mod.F90
  core_ocean/MARBL/src/marbl_diagnostics_mod.F90
  core_ocean/MARBL/src/marbl_diagnostics_share_mod.F90
  core_ocean/MARBL/src/marbl_glo_avg_mod.F90
  core_ocean/MARBL/src/marbl_init_mod.F90
  core_ocean/MARBL/src/marbl_interface.F90
  core_ocean/MARBL/src/marbl_interface_constants.F90
  core_ocean/MARBL/src/marbl_interface_private_types.F90
  core_ocean/MARBL/src/marbl_interface_public_types.F90
  core_ocean/MARBL/src/marbl_interior_tendency_mod.F90
  core_ocean/MARBL/src/marbl_interior_tendency_share_mod.F90
  core_ocean/MARBL/src/marbl_kinds_mod.F90
  core_ocean/MARBL/src/marbl_logging.F90
  core_ocean/MARBL/src/marbl_nhx_surface_emis_mod.F90
  core_ocean/MARBL/src/marbl_oxygen.F90
  core_ocean/MARBL/src/marbl_pft_mod.F90
  core_ocean/MARBL/src/marbl_restore_mod.F90
  core_ocean/MARBL/src/marbl_saved_state_mod.F90
  core_ocean/MARBL/src/marbl_schmidt_number_mod.F90
  core_ocean/MARBL/src/marbl_settings_mod.F90
  core_ocean/MARBL/src/marbl_surface_flux_mod.F90
  core_ocean/MARBL/src/marbl_surface_flux_share_mod.F90
  core_ocean/MARBL/src/marbl_temperature.F90
  core_ocean/MARBL/src/marbl_timing_mod.F90
  core_ocean/MARBL/src/marbl_utils_mod.F90
)

# Add PPR
if (NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/core_ocean/ppr/.git)
  message(FATAL_ERROR "Missing core_ocean/ppr/.git, did you forget to 'git submodule update --init --recursive' ?")
endif()
set(PPR_FILES
  core_ocean/ppr/src/ppr_1d.F90
)

# Add GOTM
if (NOT EXISTS core_ocean/gotm/.git)
  message(FATAL "Missing core_ocean/gotm/.git, did you forget to 'git submodule update --init --recursive' ?")
endif()
set(GOTM_FILES
  core_ocean/gotm/src/util/adv_center.F90
  core_ocean/gotm/src/util/convert_fluxes.F90
  core_ocean/gotm/src/util/diff_center.F90
  core_ocean/gotm/src/util/diff_face.F90
  core_ocean/gotm/src/util/eqstate.F90
  core_ocean/gotm/src/util/gridinterpol.F90
  core_ocean/gotm/src/util/lagrange.F90
  core_ocean/gotm/src/util/ode_solvers.F90
  core_ocean/gotm/src/util/time.F90
  core_ocean/gotm/src/util/tridiagonal.F90
  core_ocean/gotm/src/util/util.F90
  core_ocean/gotm/src/util/field_manager.F90
  core_ocean/gotm/src/turbulence/algebraiclength.F90
  core_ocean/gotm/src/turbulence/alpha_mnb.F90
  core_ocean/gotm/src/turbulence/cmue_a.F90
  core_ocean/gotm/src/turbulence/cmue_b.F90
  core_ocean/gotm/src/turbulence/cmue_c.F90
  core_ocean/gotm/src/turbulence/cmue_d.F90
  core_ocean/gotm/src/turbulence/cmue_ma.F90
  core_ocean/gotm/src/turbulence/cmue_rf.F90
  core_ocean/gotm/src/turbulence/cmue_sg.F90
  core_ocean/gotm/src/turbulence/compute_cpsi3.F90
  core_ocean/gotm/src/turbulence/compute_rist.F90
  core_ocean/gotm/src/turbulence/dissipationeq.F90
  core_ocean/gotm/src/turbulence/epsbalgebraic.F90
  core_ocean/gotm/src/turbulence/fk_craig.F90
  core_ocean/gotm/src/turbulence/genericeq.F90
  core_ocean/gotm/src/turbulence/internal_wave.F90
  core_ocean/gotm/src/turbulence/ispralength.F90
  core_ocean/gotm/src/turbulence/kbalgebraic.F90
  core_ocean/gotm/src/turbulence/kbeq.F90
  core_ocean/gotm/src/turbulence/lengthscaleeq.F90
  core_ocean/gotm/src/turbulence/potentialml.F90
  core_ocean/gotm/src/turbulence/production.F90
  core_ocean/gotm/src/turbulence/q2over2eq.F90
  core_ocean/gotm/src/turbulence/r_ratio.F90
  core_ocean/gotm/src/turbulence/tkealgebraic.F90
  core_ocean/gotm/src/turbulence/tkeeq.F90
  core_ocean/gotm/src/turbulence/turbulence.F90
  core_ocean/gotm/src/turbulence/variances.F90
)

list(APPEND RAW_SOURCES ${CVMIX_FILES} ${BGC_FILES} ${MARBL_FILES} ${GOTM_FILES} ${PPR_FILES})
list(APPEND NO_PREPROCESS ${CVMIX_FILES} ${BGC_FILES} ${MARBL_FILES} ${GOTM_FILES} ${PPR_FILES})

# Add analysis members
list(APPEND RAW_SOURCES
  core_ocean/analysis_members/mpas_ocn_global_stats.F
  core_ocean/analysis_members/mpas_ocn_okubo_weiss.F
  core_ocean/analysis_members/mpas_ocn_okubo_weiss_eigenvalues.c
  core_ocean/analysis_members/mpas_ocn_layer_volume_weighted_averages.F
  core_ocean/analysis_members/mpas_ocn_surface_area_weighted_averages.F
  core_ocean/analysis_members/mpas_ocn_water_mass_census.F
  core_ocean/analysis_members/mpas_ocn_meridional_heat_transport.F
  core_ocean/analysis_members/mpas_ocn_test_compute_interval.F
  core_ocean/analysis_members/mpas_ocn_high_frequency_output.F
  core_ocean/analysis_members/mpas_ocn_zonal_mean.F
  core_ocean/analysis_members/mpas_ocn_lagrangian_particle_tracking_interpolations.F
  core_ocean/analysis_members/mpas_ocn_particle_list.F
  core_ocean/analysis_members/mpas_ocn_lagrangian_particle_tracking_reset.F
  core_ocean/analysis_members/mpas_ocn_lagrangian_particle_tracking.F
  core_ocean/analysis_members/mpas_ocn_eliassen_palm.F
  core_ocean/analysis_members/mpas_ocn_time_filters.F
  core_ocean/analysis_members/mpas_ocn_mixed_layer_depths.F
  core_ocean/analysis_members/mpas_ocn_pointwise_stats.F
  core_ocean/analysis_members/mpas_ocn_debug_diagnostics.F
  core_ocean/analysis_members/mpas_ocn_time_series_stats.F
  core_ocean/analysis_members/mpas_ocn_regional_stats.F
  core_ocean/analysis_members/mpas_ocn_rpn_calculator.F
  core_ocean/analysis_members/mpas_ocn_transect_transport.F
  core_ocean/analysis_members/mpas_ocn_eddy_product_variables.F
  core_ocean/analysis_members/mpas_ocn_moc_streamfunction.F
  core_ocean/analysis_members/mpas_ocn_ocean_heat_content.F
  core_ocean/analysis_members/mpas_ocn_mixed_layer_heat_budget.F
  core_ocean/analysis_members/mpas_ocn_sediment_flux_index.F
  core_ocean/analysis_members/mpas_ocn_sediment_transport.F
  core_ocean/analysis_members/mpas_ocn_harmonic_analysis.F
  core_ocean/analysis_members/mpas_ocn_conservation_check.F
  core_ocean/analysis_members/mpas_ocn_analysis_driver.F
)

# Generate core input
handle_st_nl_gen(
  "namelist.ocean;namelist.ocean.forward mode=forward;namelist.ocean.analysis mode=analysis"
  "streams.ocean stream_list.ocean. mutable;streams.ocean.forward stream_list.ocean.forward. mutable mode=forward;streams.ocean.analysis stream_list.ocean.analysis. mutable mode=analysis"
  ${CORE_INPUT_DIR} ${CORE_BLDDIR}
)
