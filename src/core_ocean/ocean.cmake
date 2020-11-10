
# build_options.mk stuff handled here
list(APPEND CPPDEFS "-DCORE_OCEAN")
list(APPEND INCLUDES "${CMAKE_BINARY_DIR}/core_ocean/BGC" "${CMAKE_BINARY_DIR}/core_ocean/shared" "${CMAKE_BINARY_DIR}/core_ocean/analysis_members" "${CMAKE_BINARY_DIR}/core_ocean/cvmix" "${CMAKE_BINARY_DIR}/core_ocean/mode_forward" "${CMAKE_BINARY_DIR}/core_ocean/mode_analysis" "${CMAKE_BINARY_DIR}/core_ocean/mode_init")

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

  core_ocean/mode_analysis/mpas_ocn_analysis_mode.F

  core_ocean/mode_init/mpas_ocn_init_mode.F
  core_ocean/mode_init/mpas_ocn_init_spherical_utils.F
  core_ocean/mode_init/mpas_ocn_init_vertical_grids.F
  core_ocean/mode_init/mpas_ocn_init_cell_markers.F
  core_ocean/mode_init/mpas_ocn_init_interpolation.F
  core_ocean/mode_init/mpas_ocn_init_ssh_and_landIcePressure.F
  core_ocean/mode_init/mpas_ocn_init_baroclinic_channel.F
  core_ocean/mode_init/mpas_ocn_init_lock_exchange.F
  core_ocean/mode_init/mpas_ocn_init_dam_break.F
  core_ocean/mode_init/mpas_ocn_init_internal_waves.F
  core_ocean/mode_init/mpas_ocn_init_overflow.F
  core_ocean/mode_init/mpas_ocn_init_cvmix_WSwSBF.F
  core_ocean/mode_init/mpas_ocn_init_iso.F
  core_ocean/mode_init/mpas_ocn_init_soma.F
  core_ocean/mode_init/mpas_ocn_init_ziso.F
  core_ocean/mode_init/mpas_ocn_init_sub_ice_shelf_2D.F
  core_ocean/mode_init/mpas_ocn_init_periodic_planar.F
  core_ocean/mode_init/mpas_ocn_init_ecosys_column.F
  core_ocean/mode_init/mpas_ocn_init_sea_mount.F
  core_ocean/mode_init/mpas_ocn_init_global_ocean.F
  core_ocean/mode_init/mpas_ocn_init_isomip.F
  core_ocean/mode_init/mpas_ocn_init_hurricane.F
  core_ocean/mode_init/mpas_ocn_init_isomip_plus.F
  core_ocean/mode_init/mpas_ocn_init_tidal_boundary.F

  core_ocean/shared/mpas_ocn_init_routines.F
  core_ocean/shared/mpas_ocn_gm.F
  core_ocean/shared/mpas_ocn_diagnostics.F
  core_ocean/shared/mpas_ocn_diagnostics_routines.F
  core_ocean/shared/mpas_ocn_thick_ale.F
  core_ocean/shared/mpas_ocn_equation_of_state.F
  core_ocean/shared/mpas_ocn_equation_of_state_jm.F
  core_ocean/shared/mpas_ocn_equation_of_state_linear.F
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
  core_ocean/shared/mpas_ocn_vmix.F
  core_ocean/shared/mpas_ocn_vmix_coefs_const.F
  core_ocean/shared/mpas_ocn_vmix_coefs_rich.F
  core_ocean/shared/mpas_ocn_vmix_coefs_tanh.F
  core_ocean/shared/mpas_ocn_vmix_coefs_redi.F
  core_ocean/shared/mpas_ocn_vmix_cvmix.F
  core_ocean/shared/mpas_ocn_tendency.F
  core_ocean/shared/mpas_ocn_tracer_hmix.F
  core_ocean/shared/mpas_ocn_tracer_hmix_del2.F
  core_ocean/shared/mpas_ocn_tracer_hmix_del4.F
  core_ocean/shared/mpas_ocn_tracer_hmix_redi.F
  core_ocean/shared/mpas_ocn_tracer_advection.F
  core_ocean/shared/mpas_ocn_tracer_advection_mono.F
  core_ocean/shared/mpas_ocn_tracer_advection_std.F
  core_ocean/shared/mpas_ocn_tracer_nonlocalflux.F
  core_ocean/shared/mpas_ocn_tracer_short_wave_absorption.F
  core_ocean/shared/mpas_ocn_tracer_short_wave_absorption_jerlov.F
  core_ocean/shared/mpas_ocn_tracer_short_wave_absorption_variable.F
  core_ocean/shared/mpas_ocn_tracer_surface_restoring.F
  core_ocean/shared/mpas_ocn_tracer_interior_restoring.F
  core_ocean/shared/mpas_ocn_tracer_exponential_decay.F
  core_ocean/shared/mpas_ocn_tracer_ideal_age.F
  core_ocean/shared/mpas_ocn_tracer_TTD.F
  core_ocean/shared/mpas_ocn_tracer_ecosys.F
  core_ocean/shared/mpas_ocn_tracer_DMS.F
  core_ocean/shared/mpas_ocn_tracer_MacroMolecules.F
  core_ocean/shared/mpas_ocn_high_freq_thickness_hmix_del2.F
  core_ocean/shared/mpas_ocn_tracer_surface_flux_to_tend.F
  core_ocean/shared/mpas_ocn_test.F
  core_ocean/shared/mpas_ocn_constants.F
  core_ocean/shared/mpas_ocn_forcing.F
  core_ocean/shared/mpas_ocn_surface_bulk_forcing.F
  core_ocean/shared/mpas_ocn_surface_land_ice_fluxes.F
  core_ocean/shared/mpas_ocn_effective_density_in_land_ice.F
  core_ocean/shared/mpas_ocn_frazil_forcing.F
  core_ocean/shared/mpas_ocn_tidal_forcing.F
  core_ocean/shared/mpas_ocn_time_average_coupled.F
  core_ocean/shared/mpas_ocn_sea_ice.F
  core_ocean/shared/mpas_ocn_framework_forcing.F
  core_ocean/shared/mpas_ocn_time_varying_forcing.F
  core_ocean/shared/mpas_ocn_wetting_drying.F
  core_ocean/shared/mpas_ocn_tidal_potential_forcing.F
)

set(OCEAN_DRIVER
  core_ocean/driver/mpas_ocn_core.F
  core_ocean/driver/mpas_ocn_core_interface.F
)
list(APPEND RAW_SOURCES ${OCEAN_DRIVER})
list(APPEND DISABLE_QSMP ${OCEAN_DRIVER})

# Get CVMix
execute_process(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/core_ocean/get_cvmix.sh
  WORKING_DIRECTORY ${CORE_BLDDIR})

# Get BGC
execute_process(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/core_ocean/get_BGC.sh
  WORKING_DIRECTORY ${CORE_BLDDIR})

# Add CVMix
set(CVMIX_FILES
  ${CORE_BLDDIR}/cvmix/cvmix_kinds_and_types.F90
  ${CORE_BLDDIR}/cvmix/cvmix_background.F90
  ${CORE_BLDDIR}/cvmix/cvmix_convection.F90
  ${CORE_BLDDIR}/cvmix/cvmix_ddiff.F90
  ${CORE_BLDDIR}/cvmix/cvmix_kpp.F90
  ${CORE_BLDDIR}/cvmix/cvmix_math.F90
  ${CORE_BLDDIR}/cvmix/cvmix_put_get.F90
  ${CORE_BLDDIR}/cvmix/cvmix_shear.F90
  ${CORE_BLDDIR}/cvmix/cvmix_tidal.F90
  ${CORE_BLDDIR}/cvmix/cvmix_utils.F90
)

# Add BGC
set(BGC_FILES
  ${CORE_BLDDIR}/BGC/BGC_mod.F90
  ${CORE_BLDDIR}/BGC/BGC_parms.F90
  ${CORE_BLDDIR}/BGC/DMS_mod.F90
  ${CORE_BLDDIR}/BGC/DMS_parms.F90
  ${CORE_BLDDIR}/BGC/MACROS_mod.F90
  ${CORE_BLDDIR}/BGC/MACROS_parms.F90
  ${CORE_BLDDIR}/BGC/co2calc.F90
)

list(APPEND RAW_SOURCES ${CVMIX_FILES} ${BGC_FILES})
list(APPEND NO_PREPROCESS ${CVMIX_FILES} ${BGC_FILES})

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
  core_ocean/analysis_members/mpas_ocn_analysis_driver.F
)

# add accelerator/gpu flags
list(APPEND ADD_ACC_FLAGS
  core_ocean/shared/mpas_ocn_equation_of_state_jm.f90
  core_ocean/shared/mpas_ocn_mesh.f90
  core_ocean/shared/mpas_ocn_surface_bulk_forcing.f90
  core_ocean/shared/mpas_ocn_surface_land_ice_fluxes.f90
  core_ocean/shared/mpas_ocn_tendency.f90
  core_ocean/shared/mpas_ocn_vel_forcing_explicit_bottom_drag.f90
  core_ocean/shared/mpas_ocn_vel_forcing_surface_stress.f90
  core_ocean/shared/mpas_ocn_vel_hadv_coriolis.f90
  core_ocean/shared/mpas_ocn_vel_hmix_del2.f90
  core_ocean/shared/mpas_ocn_vel_hmix_del4.f90
  core_ocean/shared/mpas_ocn_vel_hmix_leith.f90
  core_ocean/shared/mpas_ocn_vel_pressure_grad.f90
  core_ocean/shared/mpas_ocn_vel_vadv.f90
)

# Generate core input
handle_st_nl_gen(
  "namelist.ocean;namelist.ocean.forward mode=forward;namelist.ocean.analysis mode=analysis;namelist.ocean.init mode=init"
  "streams.ocean stream_list.ocean. mutable;streams.ocean.forward stream_list.ocean.forward. mutable mode=forward;streams.ocean.analysis stream_list.ocean.analysis. mutable mode=analysis;streams.ocean.init stream_list.ocean.init. mutable mode=init"
  ${CORE_INPUT_DIR} ${CORE_BLDDIR}
)
