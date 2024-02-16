      character (len=StrKIND) :: config_ocean_run_mode

      logical :: config_do_restart
      character (len=StrKIND) :: config_restart_timestamp_name
      character (len=StrKIND) :: config_start_time
      character (len=StrKIND) :: config_stop_time
      character (len=StrKIND) :: config_run_duration
      character (len=StrKIND) :: config_calendar_type
      character (len=StrKIND) :: config_output_reference_time

      logical :: config_write_output_on_startup
      integer :: config_pio_num_iotasks
      integer :: config_pio_stride

      integer :: config_num_halos
      character (len=StrKIND) :: config_block_decomp_file_prefix
      integer :: config_number_of_blocks
      logical :: config_explicit_proc_decomp
      character (len=StrKIND) :: config_proc_decomp_file_prefix

      character (len=StrKIND) :: config_dt
      character (len=StrKIND) :: config_time_integrator
      integer :: config_number_of_time_levels

      logical :: config_hmix_scaleWithMesh
      real (kind=RKIND) :: config_maxMeshDensity
      logical :: config_hmix_use_ref_cell_width
      real (kind=RKIND) :: config_hmix_ref_cell_width
      real (kind=RKIND) :: config_apvm_scale_factor

      logical :: config_use_mom_del2
      real (kind=RKIND) :: config_mom_del2
      logical :: config_use_tracer_del2
      real (kind=RKIND) :: config_tracer_del2

      logical :: config_use_mom_del4
      real (kind=RKIND) :: config_mom_del4
      real (kind=RKIND) :: config_mom_del4_div_factor
      logical :: config_use_tracer_del4
      real (kind=RKIND) :: config_tracer_del4

      logical :: config_use_Leith_del2
      real (kind=RKIND) :: config_Leith_parameter
      real (kind=RKIND) :: config_Leith_dx
      real (kind=RKIND) :: config_Leith_visc2_max

      logical :: config_use_Redi
      character (len=StrKIND) :: config_Redi_closure
      real (kind=RKIND) :: config_Redi_constant_kappa
      real (kind=RKIND) :: config_Redi_maximum_slope
      logical :: config_Redi_use_slope_taper
      logical :: config_Redi_use_surface_taper
      logical :: config_Redi_limit_term1
      logical :: config_Redi_use_quasi_monotone_limiter
      real (kind=RKIND) :: config_Redi_quasi_monotone_safety_factor
      integer :: config_Redi_min_layers_diag_terms
      character (len=StrKIND) :: config_Redi_horizontal_taper
      real (kind=RKIND) :: config_Redi_horizontal_ramp_min
      real (kind=RKIND) :: config_Redi_horizontal_ramp_max

      logical :: config_submesoscale_enable
      real (kind=RKIND) :: config_submesoscale_tau
      real (kind=RKIND) :: config_submesoscale_Ce
      real (kind=RKIND) :: config_submesoscale_Lfmin
      real (kind=RKIND) :: config_submesoscale_ds_max

      logical :: config_use_GM
      character (len=StrKIND) :: config_GM_closure
      real (kind=RKIND) :: config_GM_constant_kappa
      real (kind=RKIND) :: config_GM_constant_bclModeSpeed
      character (len=StrKIND) :: config_GM_minBclModeSpeed_method
      real (kind=RKIND) :: config_GM_spatially_variable_min_kappa
      real (kind=RKIND) :: config_GM_spatially_variable_max_kappa
      real (kind=RKIND) :: config_GM_spatially_variable_baroclinic_mode
      real (kind=RKIND) :: config_GM_Visbeck_alpha
      real (kind=RKIND) :: config_GM_Visbeck_max_depth
      real (kind=RKIND) :: config_GM_EG_riMin
      real (kind=RKIND) :: config_GM_EG_kappa_factor
      real (kind=RKIND) :: config_GM_EG_Rossby_factor
      real (kind=RKIND) :: config_GM_EG_Rhines_factor
      character (len=StrKIND) :: config_GM_horizontal_taper
      real (kind=RKIND) :: config_GM_horizontal_ramp_min
      real (kind=RKIND) :: config_GM_horizontal_ramp_max
      real (kind=RKIND) :: config_GMRedi_Rossby_ramp_min
      real (kind=RKIND) :: config_GMRedi_Rossby_ramp_max

      real (kind=RKIND) :: config_eddyMLD_dens_threshold
      real (kind=RKIND) :: config_eddyMLD_reference_depth
      real (kind=RKIND) :: config_eddyMLD_reference_pressure
      logical :: config_eddyMLD_use_old

      logical :: config_use_cvmix
      real (kind=RKIND) :: config_cvmix_prandtl_number
      character (len=StrKIND) :: config_cvmix_background_scheme
      real (kind=RKIND) :: config_cvmix_background_diffusion
      real (kind=RKIND) :: config_cvmix_background_viscosity
      real (kind=RKIND) :: config_cvmix_BryanLewis_bl1
      real (kind=RKIND) :: config_cvmix_BryanLewis_bl2
      real (kind=RKIND) :: config_cvmix_BryanLewis_transitionDepth
      real (kind=RKIND) :: config_cvmix_BryanLewis_transitionWidth
      logical :: config_use_cvmix_convection
      real (kind=RKIND) :: config_cvmix_convective_diffusion
      real (kind=RKIND) :: config_cvmix_convective_viscosity
      logical :: config_cvmix_convective_basedOnBVF
      real (kind=RKIND) :: config_cvmix_convective_triggerBVF
      logical :: config_use_cvmix_shear
      integer :: config_cvmix_num_ri_smooth_loops
      logical :: config_cvmix_use_BLD_smoothing
      character (len=StrKIND) :: config_cvmix_shear_mixing_scheme
      real (kind=RKIND) :: config_cvmix_shear_PP_nu_zero
      real (kind=RKIND) :: config_cvmix_shear_PP_alpha
      real (kind=RKIND) :: config_cvmix_shear_PP_exp
      real (kind=RKIND) :: config_cvmix_shear_KPP_nu_zero
      real (kind=RKIND) :: config_cvmix_shear_KPP_Ri_zero
      real (kind=RKIND) :: config_cvmix_shear_KPP_exp
      logical :: config_use_cvmix_tidal_mixing
      logical :: config_use_cvmix_double_diffusion
      logical :: config_use_cvmix_kpp
      logical :: config_use_cvmix_fixed_boundary_layer
      real (kind=RKIND) :: config_cvmix_kpp_boundary_layer_depth
      real (kind=RKIND) :: config_cvmix_kpp_criticalBulkRichardsonNumber
      character (len=StrKIND) :: config_cvmix_kpp_matching
      logical :: config_cvmix_kpp_EkmanOBL
      logical :: config_cvmix_kpp_MonObOBL
      character (len=StrKIND) :: config_cvmix_kpp_interpolationOMLType
      real (kind=RKIND) :: config_cvmix_kpp_surface_layer_extent
      real (kind=RKIND) :: config_cvmix_kpp_surface_layer_averaging
      real (kind=RKIND) :: configure_cvmix_kpp_minimum_OBL_under_sea_ice
      real (kind=RKIND) :: config_cvmix_kpp_stop_OBL_search
      logical :: config_cvmix_kpp_use_enhanced_diff
      logical :: config_cvmix_kpp_nonlocal_with_implicit_mix
      logical :: config_cvmix_kpp_use_theory_wave
      character (len=StrKIND) :: config_cvmix_kpp_langmuir_mixing_opt
      character (len=StrKIND) :: config_cvmix_kpp_langmuir_entrainment_opt
      logical :: config_cvmix_kpp_use_active_wave

      logical :: config_use_active_wave
      integer :: config_n_stokes_drift_wavenumber_partitions

      logical :: config_use_gotm
      character (len=StrKIND) :: config_gotm_namelist_file
      real (kind=RKIND) :: config_gotm_constant_surface_roughness_length
      real (kind=RKIND) :: config_gotm_constant_bottom_roughness_length
      real (kind=RKIND) :: config_gotm_constant_bottom_drag_coeff

      logical :: config_use_variable_drag
      logical :: config_use_bulk_wind_stress
      logical :: config_use_bulk_thickness_flux
      real (kind=RKIND) :: config_flux_attenuation_coefficient
      real (kind=RKIND) :: config_flux_attenuation_coefficient_runoff

      logical :: config_use_time_varying_atmospheric_forcing
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_type
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_start_time
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_reference_time
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_cycle_start
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_cycle_duration
      character (len=StrKIND) :: config_time_varying_atmospheric_forcing_interval
      real (kind=RKIND) :: config_time_varying_atmospheric_forcing_ramp
      real (kind=RKIND) :: config_time_varying_atmospheric_forcing_ramp_delay
      logical :: config_use_time_varying_land_ice_forcing
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_start_time
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_reference_time
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_cycle_start
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_cycle_duration
      character (len=StrKIND) :: config_time_varying_land_ice_forcing_interval

      real (kind=RKIND) :: config_ssh_grad_relax_timescale
      logical :: config_remove_AIS_coupler_runoff

      character (len=StrKIND) :: config_sw_absorption_type
      integer :: config_jerlov_water_type
      real (kind=RKIND) :: config_surface_buoyancy_depth
      logical :: config_enable_shortwave_energy_fixer

      logical :: config_use_tidal_forcing
      real (kind=RKIND) :: config_use_tidal_forcing_tau
      character (len=StrKIND) :: config_tidal_forcing_type
      character (len=StrKIND) :: config_tidal_forcing_model
      real (kind=RKIND) :: config_tidal_forcing_monochromatic_amp
      real (kind=RKIND) :: config_tidal_forcing_monochromatic_period
      real (kind=RKIND) :: config_tidal_forcing_monochromatic_phaseLag
      real (kind=RKIND) :: config_tidal_forcing_monochromatic_baseline

      logical :: config_use_self_attraction_loading
      real (kind=RKIND) :: config_self_attraction_loading_smoothing_width
      character (len=StrKIND) :: config_mpas_to_grid_weights_file
      character (len=StrKIND) :: config_grid_to_mpas_weights_file
      character (len=StrKIND) :: config_self_attraction_loading_compute_interval
      integer :: config_nLatitude
      integer :: config_nLongitude
      logical :: config_use_parallel_self_attraction_loading
      integer :: config_parallel_self_attraction_loading_order
      integer :: config_parallel_self_attraction_loading_n_cells_per_block
      logical :: config_parallel_self_attraction_loading_bfb

      logical :: config_use_tidal_potential_forcing
      character (len=StrKIND) :: config_tidal_potential_reference_time
      logical :: config_use_tidal_potential_forcing_M2
      logical :: config_use_tidal_potential_forcing_S2
      logical :: config_use_tidal_potential_forcing_N2
      logical :: config_use_tidal_potential_forcing_K2
      logical :: config_use_tidal_potential_forcing_K1
      logical :: config_use_tidal_potential_forcing_O1
      logical :: config_use_tidal_potential_forcing_Q1
      logical :: config_use_tidal_potential_forcing_P1
      real (kind=RKIND) :: config_tidal_potential_ramp
      real (kind=RKIND) :: config_self_attraction_and_loading_beta

      logical :: config_use_frazil_ice_formation
      logical :: config_frazil_in_open_ocean
      logical :: config_frazil_under_land_ice
      real (kind=RKIND) :: config_frazil_heat_of_fusion
      real (kind=RKIND) :: config_frazil_ice_density
      real (kind=RKIND) :: config_frazil_fractional_thickness_limit
      real (kind=RKIND) :: config_specific_heat_sea_water
      real (kind=RKIND) :: config_frazil_maximum_depth
      real (kind=RKIND) :: config_frazil_sea_ice_reference_salinity
      real (kind=RKIND) :: config_frazil_land_ice_reference_salinity
      real (kind=RKIND) :: config_frazil_maximum_freezing_temperature
      logical :: config_frazil_use_surface_pressure

      character (len=StrKIND) :: config_land_ice_flux_mode
      character (len=StrKIND) :: config_land_ice_flux_formulation
      logical :: config_land_ice_flux_useHollandJenkinsAdvDiff
      real (kind=RKIND) :: config_land_ice_flux_attenuation_coefficient
      real (kind=RKIND) :: config_land_ice_flux_boundaryLayerThickness
      real (kind=RKIND) :: config_land_ice_flux_boundaryLayerNeighborWeight
      real (kind=RKIND) :: config_land_ice_flux_cp_ice
      real (kind=RKIND) :: config_land_ice_flux_rho_ice
      real (kind=RKIND) :: config_land_ice_flux_explicit_topDragCoeff
      real (kind=RKIND) :: config_land_ice_flux_ISOMIP_gammaT
      real (kind=RKIND) :: config_land_ice_flux_rms_tidal_velocity
      real (kind=RKIND) :: config_land_ice_flux_jenkins_heat_transfer_coefficient
      real (kind=RKIND) :: config_land_ice_flux_jenkins_salt_transfer_coefficient

      character (len=StrKIND) :: config_vert_advection_method
      integer :: config_vert_remap_order
      integer :: config_vert_remap_interval
      integer :: config_vert_tracer_adv_flux_order
      integer :: config_horiz_tracer_adv_order
      real (kind=RKIND) :: config_coef_3rd_order
      character (len=StrKIND) :: config_flux_limiter
      character (len=StrKIND) :: config_remap_limiter
      character (len=StrKIND) :: config_thickness_flux_type

      character (len=StrKIND) :: config_bottom_drag_mode
      character (len=StrKIND) :: config_implicit_bottom_drag_type
      real (kind=RKIND) :: config_implicit_constant_bottom_drag_coeff
      logical :: config_use_implicit_top_drag
      real (kind=RKIND) :: config_implicit_top_drag_coeff
      real (kind=RKIND) :: config_loglaw_bottom_roughness
      real (kind=RKIND) :: config_loglaw_layer_depth_max
      real (kind=RKIND) :: config_loglaw_bottom_drag_min
      real (kind=RKIND) :: config_loglaw_bottom_drag_max
      real (kind=RKIND) :: config_explicit_bottom_drag_coeff
      logical :: config_use_topographic_wave_drag
      real (kind=RKIND) :: config_topographic_wave_drag_coeff
      character (len=StrKIND) :: config_thickness_drag_type

      real (kind=RKIND) :: config_Rayleigh_damping_coeff
      logical :: config_Rayleigh_damping_depth_variable
      logical :: config_Rayleigh_bottom_friction
      real (kind=RKIND) :: config_Rayleigh_bottom_damping_coeff

      logical :: config_use_vegetation_drag
      real (kind=RKIND) :: config_vegetation_drag_coefficient

      logical :: config_use_wetting_drying
      logical :: config_prevent_drying
      real (kind=RKIND) :: config_drying_min_cell_height
      logical :: config_zero_drying_velocity
      logical :: config_zero_drying_velocity_ramp
      real (kind=RKIND) :: config_zero_drying_velocity_ramp_hmin
      real (kind=RKIND) :: config_zero_drying_velocity_ramp_hmax
      logical :: config_verify_not_dry
      real (kind=RKIND) :: config_drying_safety_height

      real (kind=RKIND) :: config_density0

      integer :: config_dt_scaling_LTS

      character (len=StrKIND) :: config_pressure_gradient_type
      real (kind=RKIND) :: config_common_level_weight
      real (kind=RKIND) :: config_zonal_ssh_grad
      real (kind=RKIND) :: config_meridional_ssh_grad

      character (len=StrKIND) :: config_eos_type
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_0
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_S
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_p
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_pS
      real (kind=RKIND) :: config_open_ocean_freezing_temperature_coeff_mushy_az1_liq
      real (kind=RKIND) :: config_land_ice_cavity_freezing_temperature_coeff_0
      real (kind=RKIND) :: config_land_ice_cavity_freezing_temperature_coeff_S
      real (kind=RKIND) :: config_land_ice_cavity_freezing_temperature_coeff_p
      real (kind=RKIND) :: config_land_ice_cavity_freezing_temperature_coeff_pS

      real (kind=RKIND) :: config_eos_linear_alpha
      real (kind=RKIND) :: config_eos_linear_beta
      real (kind=RKIND) :: config_eos_linear_Tref
      real (kind=RKIND) :: config_eos_linear_Sref
      real (kind=RKIND) :: config_eos_linear_densityref

      real (kind=RKIND) :: config_eos_wright_ref_pressure

      integer :: config_n_ts_iter
      integer :: config_n_bcl_iter_beg
      integer :: config_n_bcl_iter_mid
      integer :: config_n_bcl_iter_end

      character (len=StrKIND) :: config_btr_dt
      integer :: config_n_btr_cor_iter
      logical :: config_vel_correction
      integer :: config_btr_subcycle_loop_factor
      real (kind=RKIND) :: config_btr_gam1_velWt1
      real (kind=RKIND) :: config_btr_gam2_SSHWt1
      real (kind=RKIND) :: config_btr_gam3_velWt2
      logical :: config_btr_solve_SSH2

      character (len=StrKIND) :: config_btr_si_preconditioner
      real (kind=RKIND) :: config_btr_si_tolerance
      integer :: config_n_btr_si_large_iter
      logical :: config_btr_si_partition_match_mode

      character (len=StrKIND) :: config_vert_coord_movement
      character (len=StrKIND) :: config_ALE_thickness_proportionality
      real (kind=RKIND) :: config_vert_taper_weight_depth_1
      real (kind=RKIND) :: config_vert_taper_weight_depth_2
      logical :: config_use_min_max_thickness
      real (kind=RKIND) :: config_min_thickness
      real (kind=RKIND) :: config_max_thickness_factor
      logical :: config_dzdk_positive

      logical :: config_use_freq_filtered_thickness
      real (kind=RKIND) :: config_thickness_filter_timescale
      logical :: config_use_highFreqThick_restore
      real (kind=RKIND) :: config_highFreqThick_restore_time
      logical :: config_use_highFreqThick_del2
      real (kind=RKIND) :: config_highFreqThick_del2

      logical :: config_check_zlevel_consistency
      logical :: config_check_ssh_consistency
      logical :: config_filter_btr_mode
      logical :: config_prescribe_velocity
      logical :: config_prescribe_thickness
      logical :: config_include_KE_vertex
      logical :: config_check_tracer_monotonicity
      logical :: config_compute_active_tracer_budgets
      logical :: config_disable_thick_all_tend
      logical :: config_disable_thick_hadv
      logical :: config_disable_thick_vadv
      logical :: config_disable_thick_sflux
      logical :: config_disable_vel_all_tend
      logical :: config_disable_vel_hadv
      logical :: config_disable_vel_coriolis
      logical :: config_disable_vel_pgrad
      logical :: config_disable_vel_hmix
      logical :: config_disable_vel_surface_stress
      logical :: config_disable_vel_topographic_wave_drag
      logical :: config_disable_vel_explicit_bottom_drag
      logical :: config_disable_vel_vmix
      logical :: config_disable_vel_vadv
      logical :: config_disable_tr_all_tend
      logical :: config_disable_tr_adv
      logical :: config_disable_tr_hmix
      logical :: config_disable_tr_vmix
      logical :: config_disable_tr_sflux
      logical :: config_disable_tr_nonlocalflux
      logical :: config_disable_redi_k33
      logical :: config_read_nearest_restart

      logical :: config_conduct_tests
      logical :: config_test_tensors
      character (len=StrKIND) :: config_tensor_test_function

      integer :: config_transport_tests_vert_levels
      real (kind=RKIND) :: config_transport_tests_temperature
      real (kind=RKIND) :: config_transport_tests_salinity
      integer :: config_transport_tests_flow_id

      integer :: config_vert_levels

      logical :: config_use_manufactured_solution
      real (kind=RKIND) :: config_manufactured_solution_wavelength_x
      real (kind=RKIND) :: config_manufactured_solution_wavelength_y
      real (kind=RKIND) :: config_manufactured_solution_amplitude

      logical :: config_use_activeTracers
      logical :: config_use_activeTracers_surface_bulk_forcing
      logical :: config_use_activeTracers_surface_restoring
      logical :: config_use_activeTracers_interior_restoring
      logical :: config_use_activeTracers_exponential_decay
      logical :: config_use_activeTracers_idealAge_forcing
      logical :: config_use_activeTracers_ttd_forcing
      logical :: config_use_surface_salinity_monthly_restoring
      character (len=StrKIND) :: config_surface_salinity_monthly_restoring_compute_interval
      real (kind=RKIND) :: config_salinity_restoring_constant_piston_velocity
      real (kind=RKIND) :: config_salinity_restoring_max_difference
      logical :: config_salinity_restoring_under_sea_ice

      logical :: config_use_debugTracers
      logical :: config_reset_debugTracers_near_surface
      integer :: config_reset_debugTracers_top_nLayers
      logical :: config_use_debugTracers_surface_bulk_forcing
      logical :: config_use_debugTracers_surface_restoring
      logical :: config_use_debugTracers_interior_restoring
      logical :: config_use_debugTracers_exponential_decay
      logical :: config_use_debugTracers_idealAge_forcing
      logical :: config_use_debugTracers_ttd_forcing

      logical :: config_use_ecosysTracers
      character (len=StrKIND) :: config_ecosys_atm_co2_option
      character (len=StrKIND) :: config_ecosys_atm_alt_co2_option
      logical :: config_ecosys_atm_alt_co2_use_eco
      real (kind=RKIND) :: config_ecosys_atm_co2_constant_value
      logical :: config_use_ecosysTracers_surface_bulk_forcing
      logical :: config_use_ecosysTracers_surface_restoring
      logical :: config_use_ecosysTracers_interior_restoring
      logical :: config_use_ecosysTracers_exponential_decay
      logical :: config_use_ecosysTracers_idealAge_forcing
      logical :: config_use_ecosysTracers_ttd_forcing
      logical :: config_use_ecosysTracers_surface_value
      logical :: config_use_ecosysTracers_river_inputs_from_coupler
      logical :: config_use_ecosysTracers_sea_ice_coupling
      logical :: config_ecosysTracers_diagnostic_fields_level1
      logical :: config_ecosysTracers_diagnostic_fields_level2
      logical :: config_ecosysTracers_diagnostic_fields_level3
      logical :: config_ecosysTracers_diagnostic_fields_level4
      logical :: config_ecosysTracers_diagnostic_fields_level5

      logical :: config_use_DMSTracers
      logical :: config_use_DMSTracers_surface_bulk_forcing
      logical :: config_use_DMSTracers_surface_restoring
      logical :: config_use_DMSTracers_interior_restoring
      logical :: config_use_DMSTracers_exponential_decay
      logical :: config_use_DMSTracers_idealAge_forcing
      logical :: config_use_DMSTracers_ttd_forcing
      logical :: config_use_DMSTracers_surface_value
      logical :: config_use_DMSTracers_sea_ice_coupling

      logical :: config_use_MacroMoleculesTracers
      logical :: config_use_MacroMoleculesTracers_surface_bulk_forcing
      logical :: config_use_MacroMoleculesTracers_surface_restoring
      logical :: config_use_MacroMoleculesTracers_interior_restoring
      logical :: config_use_MacroMoleculesTracers_exponential_decay
      logical :: config_use_MacroMoleculesTracers_idealAge_forcing
      logical :: config_use_MacroMoleculesTracers_ttd_forcing
      logical :: config_use_MacroMoleculesTracers_surface_value
      logical :: config_use_MacroMoleculesTracers_sea_ice_coupling

      logical :: config_use_idealAgeTracers
      logical :: config_use_idealAgeTracers_surface_bulk_forcing
      logical :: config_use_idealAgeTracers_surface_restoring
      logical :: config_use_idealAgeTracers_interior_restoring
      logical :: config_use_idealAgeTracers_exponential_decay
      logical :: config_use_idealAgeTracers_idealAge_forcing
      logical :: config_use_idealAgeTracers_ttd_forcing

      logical :: config_use_CFCTracers
      logical :: config_use_CFCTracers_surface_bulk_forcing
      logical :: config_use_CFCTracers_surface_restoring
      logical :: config_use_CFCTracers_interior_restoring
      logical :: config_use_CFCTracers_exponential_decay
      logical :: config_use_CFCTracers_idealAge_forcing
      logical :: config_use_CFCTracers_ttd_forcing
      logical :: config_use_CFC11
      logical :: config_use_CFC12

      logical :: config_AM_globalStats_enable
      character (len=StrKIND) :: config_AM_globalStats_compute_interval
      logical :: config_AM_globalStats_compute_on_startup
      logical :: config_AM_globalStats_write_on_startup
      logical :: config_AM_globalStats_text_file
      character (len=StrKIND) :: config_AM_globalStats_directory
      character (len=StrKIND) :: config_AM_globalStats_output_stream

      logical :: config_AM_surfaceAreaWeightedAverages_enable
      logical :: config_AM_surfaceAreaWeightedAverages_compute_on_startup
      logical :: config_AM_surfaceAreaWeightedAverages_write_on_startup
      character (len=StrKIND) :: config_AM_surfaceAreaWeightedAverages_compute_interval
      character (len=StrKIND) :: config_AM_surfaceAreaWeightedAverages_output_stream

      logical :: config_AM_waterMassCensus_enable
      character (len=StrKIND) :: config_AM_waterMassCensus_compute_interval
      character (len=StrKIND) :: config_AM_waterMassCensus_output_stream
      logical :: config_AM_waterMassCensus_compute_on_startup
      logical :: config_AM_waterMassCensus_write_on_startup
      real (kind=RKIND) :: config_AM_waterMassCensus_minTemperature
      real (kind=RKIND) :: config_AM_waterMassCensus_maxTemperature
      real (kind=RKIND) :: config_AM_waterMassCensus_minSalinity
      real (kind=RKIND) :: config_AM_waterMassCensus_maxSalinity
      logical :: config_AM_waterMassCensus_compute_predefined_regions
      character (len=StrKIND) :: config_AM_waterMassCensus_region_group

      logical :: config_AM_layerVolumeWeightedAverage_enable
      character (len=StrKIND) :: config_AM_layerVolumeWeightedAverage_compute_interval
      logical :: config_AM_layerVolumeWeightedAverage_compute_on_startup
      logical :: config_AM_layerVolumeWeightedAverage_write_on_startup
      character (len=StrKIND) :: config_AM_layerVolumeWeightedAverage_output_stream

      logical :: config_AM_zonalMean_enable
      logical :: config_AM_zonalMean_compute_on_startup
      logical :: config_AM_zonalMean_write_on_startup
      character (len=StrKIND) :: config_AM_zonalMean_compute_interval
      character (len=StrKIND) :: config_AM_zonalMean_output_stream
      integer :: config_AM_zonalMean_num_bins
      real (kind=RKIND) :: config_AM_zonalMean_min_bin
      real (kind=RKIND) :: config_AM_zonalMean_max_bin

      logical :: config_AM_okuboWeiss_enable
      logical :: config_AM_okuboWeiss_compute_on_startup
      logical :: config_AM_okuboWeiss_write_on_startup
      character (len=StrKIND) :: config_AM_okuboWeiss_compute_interval
      character (len=StrKIND) :: config_AM_okuboWeiss_output_stream
      character (len=StrKIND) :: config_AM_okuboWeiss_directory
      real (kind=RKIND) :: config_AM_okuboWeiss_threshold_value
      real (kind=RKIND) :: config_AM_okuboWeiss_normalization
      real (kind=RKIND) :: config_AM_okuboWeiss_lambda2_normalization
      logical :: config_AM_okuboWeiss_use_lat_lon_coords
      logical :: config_AM_okuboWeiss_compute_eddy_census
      integer :: config_AM_okuboWeiss_eddy_min_cells

      logical :: config_AM_meridionalHeatTransport_enable
      character (len=StrKIND) :: config_AM_meridionalHeatTransport_compute_interval
      logical :: config_AM_meridionalHeatTransport_compute_on_startup
      logical :: config_AM_meridionalHeatTransport_write_on_startup
      character (len=StrKIND) :: config_AM_meridionalHeatTransport_output_stream
      integer :: config_AM_meridionalHeatTransport_num_bins
      real (kind=RKIND) :: config_AM_meridionalHeatTransport_min_bin
      real (kind=RKIND) :: config_AM_meridionalHeatTransport_max_bin
      character (len=StrKIND) :: config_AM_meridionalHeatTransport_region_group

      logical :: config_AM_testComputeInterval_enable
      character (len=StrKIND) :: config_AM_testComputeInterval_compute_interval
      logical :: config_AM_testComputeInterval_compute_on_startup
      logical :: config_AM_testComputeInterval_write_on_startup
      character (len=StrKIND) :: config_AM_testComputeInterval_output_stream

      logical :: config_AM_highFrequencyOutput_enable
      character (len=StrKIND) :: config_AM_highFrequencyOutput_compute_interval
      character (len=StrKIND) :: config_AM_highFrequencyOutput_output_stream
      logical :: config_AM_highFrequencyOutput_compute_on_startup
      logical :: config_AM_highFrequencyOutput_write_on_startup

      logical :: config_AM_timeFilters_enable
      character (len=StrKIND) :: config_AM_timeFilters_compute_interval
      character (len=StrKIND) :: config_AM_timeFilters_output_stream
      character (len=StrKIND) :: config_AM_timeFilters_restart_stream
      logical :: config_AM_timeFilters_compute_on_startup
      logical :: config_AM_timeFilters_write_on_startup
      logical :: config_AM_timeFilters_initialize_filters
      character (len=StrKIND) :: config_AM_timeFilters_tau
      logical :: config_AM_timeFilters_compute_cell_centered_values

      logical :: config_AM_lagrPartTrack_enable
      character (len=StrKIND) :: config_AM_lagrPartTrack_compute_interval
      logical :: config_AM_lagrPartTrack_compute_on_startup
      character (len=StrKIND) :: config_AM_lagrPartTrack_output_stream
      character (len=StrKIND) :: config_AM_lagrPartTrack_restart_stream
      character (len=StrKIND) :: config_AM_lagrPartTrack_input_stream
      logical :: config_AM_lagrPartTrack_write_on_startup
      integer :: config_AM_lagrPartTrack_filter_number
      integer :: config_AM_lagrPartTrack_timeIntegration
      character (len=StrKIND) :: config_AM_lagrPartTrack_reset_criteria
      character (len=StrKIND) :: config_AM_lagrPartTrack_reset_global_timestamp
      character (len=StrKIND) :: config_AM_lagrPartTrack_region_stream
      logical :: config_AM_lagrPartTrack_reset_if_outside_region
      logical :: config_AM_lagrPartTrack_reset_if_inside_region
      logical :: config_AM_lagrPartTrack_sample_horizontal_interp
      logical :: config_AM_lagrPartTrack_sample_temperature
      logical :: config_AM_lagrPartTrack_sample_salinity
      logical :: config_AM_lagrPartTrack_sample_DIC
      logical :: config_AM_lagrPartTrack_sample_ALK
      logical :: config_AM_lagrPartTrack_sample_PO4
      logical :: config_AM_lagrPartTrack_sample_NO3
      logical :: config_AM_lagrPartTrack_sample_SiO3
      logical :: config_AM_lagrPartTrack_sample_NH4
      logical :: config_AM_lagrPartTrack_sample_Fe
      logical :: config_AM_lagrPartTrack_sample_O2

      logical :: config_AM_eliassenPalm_enable
      character (len=StrKIND) :: config_AM_eliassenPalm_compute_interval
      character (len=StrKIND) :: config_AM_eliassenPalm_output_stream
      character (len=StrKIND) :: config_AM_eliassenPalm_restart_stream
      logical :: config_AM_eliassenPalm_compute_on_startup
      logical :: config_AM_eliassenPalm_write_on_startup
      logical :: config_AM_eliassenPalm_debug
      integer :: config_AM_eliassenPalm_nBuoyancyLayers
      real (kind=RKIND) :: config_AM_eliassenPalm_rhomin_buoycoor
      real (kind=RKIND) :: config_AM_eliassenPalm_rhomax_buoycoor

      logical :: config_AM_mixedLayerDepths_enable
      character (len=StrKIND) :: config_AM_mixedLayerDepths_compute_interval
      character (len=StrKIND) :: config_AM_mixedLayerDepths_output_stream
      logical :: config_AM_mixedLayerDepths_write_on_startup
      logical :: config_AM_mixedLayerDepths_compute_on_startup
      logical :: config_AM_mixedLayerDepths_Tthreshold
      real (kind=RKIND) :: config_AM_mixedLayerDepths_crit_temp_threshold
      real (kind=RKIND) :: config_AM_mixedLayerDepths_reference_pressure
      logical :: config_AM_mixedLayerDepths_Tgradient
      logical :: config_AM_mixedLayerDepths_Dgradient
      real (kind=RKIND) :: config_AM_mixedLayerDepths_temp_gradient_threshold
      real (kind=RKIND) :: config_AM_mixedLayerDepths_den_gradient_threshold
      integer :: config_AM_mixedLayerDepths_interp_method

      logical :: config_AM_regionalStatsDaily_enable
      logical :: config_AM_regionalStatsDaily_compute_on_startup
      logical :: config_AM_regionalStatsDaily_write_on_startup
      character (len=StrKIND) :: config_AM_regionalStatsDaily_compute_interval
      character (len=StrKIND) :: config_AM_regionalStatsDaily_output_stream
      character (len=StrKIND) :: config_AM_regionalStatsDaily_restart_stream
      character (len=StrKIND) :: config_AM_regionalStatsDaily_input_stream
      character (len=StrKIND) :: config_AM_regionalStatsDaily_operation
      character (len=StrKIND) :: config_AM_regionalStatsDaily_region_type
      character (len=StrKIND) :: config_AM_regionalStatsDaily_region_group
      character (len=StrKIND) :: config_AM_regionalStatsDaily_1d_weighting_function
      character (len=StrKIND) :: config_AM_regionalStatsDaily_2d_weighting_function
      character (len=StrKIND) :: config_AM_regionalStatsDaily_1d_weighting_field
      character (len=StrKIND) :: config_AM_regionalStatsDaily_2d_weighting_field
      character (len=StrKIND) :: config_AM_regionalStatsDaily_vertical_mask
      character (len=StrKIND) :: config_AM_regionalStatsDaily_vertical_dimension

      logical :: config_AM_regionalStatsWeekly_enable
      logical :: config_AM_regionalStatsWeekly_compute_on_startup
      logical :: config_AM_regionalStatsWeekly_write_on_startup
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_compute_interval
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_output_stream
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_restart_stream
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_input_stream
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_operation
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_region_type
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_region_group
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_1d_weighting_function
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_2d_weighting_function
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_1d_weighting_field
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_2d_weighting_field
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_vertical_mask
      character (len=StrKIND) :: config_AM_regionalStatsWeekly_vertical_dimension

      logical :: config_AM_regionalStatsMonthly_enable
      logical :: config_AM_regionalStatsMonthly_compute_on_startup
      logical :: config_AM_regionalStatsMonthly_write_on_startup
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_compute_interval
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_output_stream
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_restart_stream
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_input_stream
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_operation
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_region_type
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_region_group
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_1d_weighting_function
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_2d_weighting_function
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_1d_weighting_field
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_2d_weighting_field
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_vertical_mask
      character (len=StrKIND) :: config_AM_regionalStatsMonthly_vertical_dimension

      logical :: config_AM_regionalStatsCustom_enable
      logical :: config_AM_regionalStatsCustom_compute_on_startup
      logical :: config_AM_regionalStatsCustom_write_on_startup
      character (len=StrKIND) :: config_AM_regionalStatsCustom_compute_interval
      character (len=StrKIND) :: config_AM_regionalStatsCustom_output_stream
      character (len=StrKIND) :: config_AM_regionalStatsCustom_restart_stream
      character (len=StrKIND) :: config_AM_regionalStatsCustom_input_stream
      character (len=StrKIND) :: config_AM_regionalStatsCustom_operation
      character (len=StrKIND) :: config_AM_regionalStatsCustom_region_type
      character (len=StrKIND) :: config_AM_regionalStatsCustom_region_group
      character (len=StrKIND) :: config_AM_regionalStatsCustom_1d_weighting_function
      character (len=StrKIND) :: config_AM_regionalStatsCustom_2d_weighting_function
      character (len=StrKIND) :: config_AM_regionalStatsCustom_1d_weighting_field
      character (len=StrKIND) :: config_AM_regionalStatsCustom_2d_weighting_field
      character (len=StrKIND) :: config_AM_regionalStatsCustom_vertical_mask
      character (len=StrKIND) :: config_AM_regionalStatsCustom_vertical_dimension

      logical :: config_AM_timeSeriesStatsDaily_enable
      logical :: config_AM_timeSeriesStatsDaily_compute_on_startup
      logical :: config_AM_timeSeriesStatsDaily_write_on_startup
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_compute_interval
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_output_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_restart_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_operation
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_reference_times
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_duration_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_repeat_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_reset_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsDaily_backward_output_offset

      logical :: config_AM_timeSeriesStatsMonthly_enable
      logical :: config_AM_timeSeriesStatsMonthly_compute_on_startup
      logical :: config_AM_timeSeriesStatsMonthly_write_on_startup
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_compute_interval
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_output_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_restart_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_operation
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_reference_times
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_duration_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_repeat_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_reset_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthly_backward_output_offset

      logical :: config_AM_timeSeriesStatsClimatology_enable
      logical :: config_AM_timeSeriesStatsClimatology_compute_on_startup
      logical :: config_AM_timeSeriesStatsClimatology_write_on_startup
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_compute_interval
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_output_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_restart_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_operation
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_reference_times
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_duration_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_repeat_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_reset_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsClimatology_backward_output_offset

      logical :: config_AM_timeSeriesStatsMonthlyMax_enable
      logical :: config_AM_timeSeriesStatsMonthlyMax_compute_on_startup
      logical :: config_AM_timeSeriesStatsMonthlyMax_write_on_startup
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_compute_interval
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_output_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_restart_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_operation
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_reference_times
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_duration_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_repeat_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_reset_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMax_backward_output_offset

      logical :: config_AM_timeSeriesStatsMonthlyMin_enable
      logical :: config_AM_timeSeriesStatsMonthlyMin_compute_on_startup
      logical :: config_AM_timeSeriesStatsMonthlyMin_write_on_startup
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_compute_interval
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_output_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_restart_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_operation
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_reference_times
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_duration_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_repeat_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_reset_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsMonthlyMin_backward_output_offset

      logical :: config_AM_timeSeriesStatsCustom_enable
      logical :: config_AM_timeSeriesStatsCustom_compute_on_startup
      logical :: config_AM_timeSeriesStatsCustom_write_on_startup
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_compute_interval
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_output_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_restart_stream
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_operation
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_reference_times
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_duration_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_repeat_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_reset_intervals
      character (len=StrKIND) :: config_AM_timeSeriesStatsCustom_backward_output_offset

      logical :: config_AM_pointwiseStats_enable
      character (len=StrKIND) :: config_AM_pointwiseStats_compute_interval
      character (len=StrKIND) :: config_AM_pointwiseStats_output_stream
      logical :: config_AM_pointwiseStats_compute_on_startup
      logical :: config_AM_pointwiseStats_write_on_startup

      logical :: config_AM_debugDiagnostics_enable
      character (len=StrKIND) :: config_AM_debugDiagnostics_compute_interval
      character (len=StrKIND) :: config_AM_debugDiagnostics_output_stream
      logical :: config_AM_debugDiagnostics_compute_on_startup
      logical :: config_AM_debugDiagnostics_write_on_startup
      logical :: config_AM_debugDiagnostics_check_state

      logical :: config_AM_rpnCalculator_enable
      logical :: config_AM_rpnCalculator_compute_on_startup
      logical :: config_AM_rpnCalculator_write_on_startup
      character (len=StrKIND) :: config_AM_rpnCalculator_compute_interval
      character (len=StrKIND) :: config_AM_rpnCalculator_output_stream
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_a
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_b
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_c
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_d
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_e
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_f
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_g
      character (len=StrKIND) :: config_AM_rpnCalculator_variable_h
      character (len=StrKIND) :: config_AM_rpnCalculator_expression_1
      character (len=StrKIND) :: config_AM_rpnCalculator_expression_2
      character (len=StrKIND) :: config_AM_rpnCalculator_expression_3
      character (len=StrKIND) :: config_AM_rpnCalculator_expression_4
      character (len=StrKIND) :: config_AM_rpnCalculator_output_name_1
      character (len=StrKIND) :: config_AM_rpnCalculator_output_name_2
      character (len=StrKIND) :: config_AM_rpnCalculator_output_name_3
      character (len=StrKIND) :: config_AM_rpnCalculator_output_name_4

      logical :: config_AM_transectTransport_enable
      character (len=StrKIND) :: config_AM_transectTransport_compute_interval
      character (len=StrKIND) :: config_AM_transectTransport_output_stream
      logical :: config_AM_transectTransport_compute_on_startup
      logical :: config_AM_transectTransport_write_on_startup
      character (len=StrKIND) :: config_AM_transectTransport_transect_group

      logical :: config_AM_eddyProductVariables_enable
      character (len=StrKIND) :: config_AM_eddyProductVariables_compute_interval
      character (len=StrKIND) :: config_AM_eddyProductVariables_output_stream
      logical :: config_AM_eddyProductVariables_compute_on_startup
      logical :: config_AM_eddyProductVariables_write_on_startup

      logical :: config_AM_mocStreamfunction_enable
      character (len=StrKIND) :: config_AM_mocStreamfunction_compute_interval
      character (len=StrKIND) :: config_AM_mocStreamfunction_output_stream
      logical :: config_AM_mocStreamfunction_compute_on_startup
      logical :: config_AM_mocStreamfunction_write_on_startup
      real (kind=RKIND) :: config_AM_mocStreamfunction_min_bin
      real (kind=RKIND) :: config_AM_mocStreamfunction_max_bin
      integer :: config_AM_mocStreamfunction_num_bins
      character (len=StrKIND) :: config_AM_mocStreamfunction_region_group
      character (len=StrKIND) :: config_AM_mocStreamfunction_transect_group

      logical :: config_AM_oceanHeatContent_enable
      character (len=StrKIND) :: config_AM_oceanHeatContent_compute_interval
      character (len=StrKIND) :: config_AM_oceanHeatContent_output_stream
      logical :: config_AM_oceanHeatContent_compute_on_startup
      logical :: config_AM_oceanHeatContent_write_on_startup

      logical :: config_AM_mixedLayerHeatBudget_enable
      character (len=StrKIND) :: config_AM_mixedLayerHeatBudget_compute_interval
      character (len=StrKIND) :: config_AM_mixedLayerHeatBudget_output_stream
      logical :: config_AM_mixedLayerHeatBudget_compute_on_startup
      logical :: config_AM_mixedLayerHeatBudget_write_on_startup

      logical :: config_AM_sedimentFluxIndex_enable
      logical :: config_AM_sedimentFluxIndex_compute_on_startup
      logical :: config_AM_sedimentFluxIndex_write_on_startup
      character (len=StrKIND) :: config_AM_sedimentFluxIndex_compute_interval
      character (len=StrKIND) :: config_AM_sedimentFluxIndex_output_stream
      character (len=StrKIND) :: config_AM_sedimentFluxIndex_directory
      logical :: config_AM_sedimentFluxIndex_use_lat_lon_coords

      logical :: config_AM_sedimentTransport_enable
      logical :: config_AM_sedimentTransport_compute_on_startup
      logical :: config_AM_sedimentTransport_write_on_startup
      character (len=StrKIND) :: config_AM_sedimentTransport_compute_interval
      character (len=StrKIND) :: config_AM_sedimentTransport_output_stream
      character (len=StrKIND) :: config_AM_sedimentTransport_directory
      real (kind=RKIND) :: config_AM_sedimentTransport_grain_size
      character (len=StrKIND) :: config_AM_sedimentTransport_ws_formula
      character (len=StrKIND) :: config_AM_sedimentTransport_bedld_formula
      character (len=StrKIND) :: config_AM_sedimentTransport_SSC_ref_formula
      real (kind=RKIND) :: config_AM_sedimentTransport_drag_coefficient
      real (kind=RKIND) :: config_AM_sedimentTransport_erate
      real (kind=RKIND) :: config_AM_sedimentTransport_tau_ce
      real (kind=RKIND) :: config_AM_sedimentTransport_tau_cd
      real (kind=RKIND) :: config_AM_sedimentTransport_Manning_coef
      real (kind=RKIND) :: config_AM_sedimentTransport_grain_porosity
      real (kind=RKIND) :: config_AM_sedimentTransport_water_density
      real (kind=RKIND) :: config_AM_sedimentTransport_grain_density
      real (kind=RKIND) :: config_AM_sedimentTransport_alpha
      real (kind=RKIND) :: config_AM_sedimentTransport_kinematic_viscosity
      real (kind=RKIND) :: config_AM_sedimentTransport_vertical_diffusion_coefficient
      logical :: config_AM_sedimentTransport_bedload
      logical :: config_AM_sedimentTransport_suspended
      logical :: config_AM_sedimentTransport_use_lat_lon_coords

      logical :: config_AM_harmonicAnalysis_enable
      character (len=StrKIND) :: config_AM_harmonicAnalysis_compute_interval
      real (kind=RKIND) :: config_AM_harmonicAnalysis_start_delay
      real (kind=RKIND) :: config_AM_harmonicAnalysis_duration
      character (len=StrKIND) :: config_AM_harmonicAnalysis_output_stream
      character (len=StrKIND) :: config_AM_harmonicAnalysis_restart_stream
      logical :: config_AM_harmonicAnalysis_compute_on_startup
      logical :: config_AM_harmonicAnalysis_write_on_startup
      logical :: config_AM_harmonicAnalysis_use_M2
      logical :: config_AM_harmonicAnalysis_use_S2
      logical :: config_AM_harmonicAnalysis_use_N2
      logical :: config_AM_harmonicAnalysis_use_K2
      logical :: config_AM_harmonicAnalysis_use_K1
      logical :: config_AM_harmonicAnalysis_use_O1
      logical :: config_AM_harmonicAnalysis_use_Q1
      logical :: config_AM_harmonicAnalysis_use_P1

      logical :: config_AM_conservationCheck_enable
      character (len=StrKIND) :: config_AM_conservationCheck_compute_interval
      character (len=StrKIND) :: config_AM_conservationCheck_output_stream
      logical :: config_AM_conservationCheck_compute_on_startup
      logical :: config_AM_conservationCheck_write_on_startup
      logical :: config_AM_conservationCheck_write_to_logfile
      character (len=StrKIND) :: config_AM_conservationCheck_restart_stream

      integer :: config_baroclinic_channel_vert_levels
      logical :: config_baroclinic_channel_use_distances
      real (kind=RKIND) :: config_baroclinic_channel_surface_temperature
      real (kind=RKIND) :: config_baroclinic_channel_bottom_temperature
      real (kind=RKIND) :: config_baroclinic_channel_temperature_difference
      real (kind=RKIND) :: config_baroclinic_channel_gradient_width_frac
      real (kind=RKIND) :: config_baroclinic_channel_gradient_width_dist
      real (kind=RKIND) :: config_baroclinic_channel_bottom_depth
      real (kind=RKIND) :: config_baroclinic_channel_salinity
      real (kind=RKIND) :: config_baroclinic_channel_coriolis_parameter

      integer :: config_lock_exchange_vert_levels
      real (kind=RKIND) :: config_lock_exchange_bottom_depth
      real (kind=RKIND) :: config_lock_exchange_cold_temperature
      real (kind=RKIND) :: config_lock_exchange_warm_temperature
      character (len=StrKIND) :: config_lock_exchange_direction
      real (kind=RKIND) :: config_lock_exchange_salinity
      character (len=StrKIND) :: config_lock_exchange_layer_type
      real (kind=RKIND) :: config_lock_exchange_isopycnal_min_thickness

      integer :: config_internal_waves_vert_levels
      logical :: config_internal_waves_use_distances
      real (kind=RKIND) :: config_internal_waves_surface_temperature
      real (kind=RKIND) :: config_internal_waves_bottom_temperature
      real (kind=RKIND) :: config_internal_waves_temperature_difference
      real (kind=RKIND) :: config_internal_waves_amplitude_width_frac
      real (kind=RKIND) :: config_internal_waves_amplitude_width_dist
      real (kind=RKIND) :: config_internal_waves_bottom_depth
      real (kind=RKIND) :: config_internal_waves_salinity
      character (len=StrKIND) :: config_internal_waves_layer_type
      real (kind=RKIND) :: config_internal_waves_isopycnal_displacement

      integer :: config_overflow_vert_levels
      logical :: config_overflow_use_distances
      real (kind=RKIND) :: config_overflow_bottom_depth
      real (kind=RKIND) :: config_overflow_ridge_depth
      real (kind=RKIND) :: config_overflow_plug_temperature
      real (kind=RKIND) :: config_overflow_domain_temperature
      real (kind=RKIND) :: config_overflow_salinity
      real (kind=RKIND) :: config_overflow_plug_width_frac
      real (kind=RKIND) :: config_overflow_slope_center_frac
      real (kind=RKIND) :: config_overflow_slope_width_frac
      real (kind=RKIND) :: config_overflow_plug_width_dist
      real (kind=RKIND) :: config_overflow_slope_center_dist
      real (kind=RKIND) :: config_overflow_slope_width_dist
      character (len=StrKIND) :: config_overflow_layer_type
      real (kind=RKIND) :: config_overflow_isopycnal_min_thickness

      integer :: config_dam_break_vert_levels
      real (kind=RKIND) :: config_dam_break_eta0
      real (kind=RKIND) :: config_dam_break_dc
      real (kind=RKIND) :: config_dam_break_R0
      real (kind=RKIND) :: config_dam_break_Xl
      real (kind=RKIND) :: config_dam_break_Yl
      real (kind=RKIND) :: config_dam_break_Inlet

      real (kind=RKIND) :: config_global_ocean_minimum_depth
      character (len=StrKIND) :: config_global_ocean_depth_file
      character (len=StrKIND) :: config_global_ocean_depth_dimname
      character (len=StrKIND) :: config_global_ocean_depth_varname
      real (kind=RKIND) :: config_global_ocean_depth_conversion_factor
      character (len=StrKIND) :: config_global_ocean_temperature_file
      character (len=StrKIND) :: config_global_ocean_salinity_file
      character (len=StrKIND) :: config_global_ocean_tracer_nlat_dimname
      character (len=StrKIND) :: config_global_ocean_tracer_nlon_dimname
      character (len=StrKIND) :: config_global_ocean_tracer_ndepth_dimname
      real (kind=RKIND) :: config_global_ocean_tracer_depth_conversion_factor
      integer :: config_global_ocean_tracer_vert_levels
      character (len=StrKIND) :: config_global_ocean_temperature_varname
      character (len=StrKIND) :: config_global_ocean_salinity_varname
      logical :: config_global_ocean_tracer_latlon_degrees
      character (len=StrKIND) :: config_global_ocean_tracer_lat_varname
      character (len=StrKIND) :: config_global_ocean_tracer_lon_varname
      character (len=StrKIND) :: config_global_ocean_tracer_depth_varname
      character (len=StrKIND) :: config_global_ocean_tracer_method
      integer :: config_global_ocean_smooth_TS_iterations
      character (len=StrKIND) :: config_global_ocean_swData_file
      character (len=StrKIND) :: config_global_ocean_swData_nlat_dimname
      character (len=StrKIND) :: config_global_ocean_swData_nlon_dimname
      character (len=StrKIND) :: config_global_ocean_swData_lat_varname
      character (len=StrKIND) :: config_global_ocean_swData_lon_varname
      logical :: config_global_ocean_swData_latlon_degrees
      character (len=StrKIND) :: config_global_ocean_swData_method
      character (len=StrKIND) :: config_global_ocean_chlorophyll_varname
      character (len=StrKIND) :: config_global_ocean_zenithAngle_varname
      character (len=StrKIND) :: config_global_ocean_clearSky_varname
      real (kind=RKIND) :: config_global_ocean_piston_velocity
      real (kind=RKIND) :: config_global_ocean_interior_restore_rate
      character (len=StrKIND) :: config_global_ocean_topography_source
      character (len=StrKIND) :: config_global_ocean_topography_file
      character (len=StrKIND) :: config_global_ocean_topography_nlat_dimname
      character (len=StrKIND) :: config_global_ocean_topography_nlon_dimname
      logical :: config_global_ocean_topography_latlon_degrees
      character (len=StrKIND) :: config_global_ocean_topography_lat_varname
      character (len=StrKIND) :: config_global_ocean_topography_lon_varname
      character (len=StrKIND) :: config_global_ocean_topography_varname
      logical :: config_global_ocean_topography_has_ocean_frac
      character (len=StrKIND) :: config_global_ocean_topography_ocean_frac_varname
      character (len=StrKIND) :: config_global_ocean_topography_method
      logical :: config_global_ocean_fill_bathymetry_holes
      integer :: config_global_ocean_topography_smooth_iterations
      real (kind=RKIND) :: config_global_ocean_topography_smooth_weight
      logical :: config_global_ocean_deepen_critical_passages
      logical :: config_global_ocean_depress_by_land_ice
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_file
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_nlat_dimname
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_nlon_dimname
      logical :: config_global_ocean_land_ice_topo_latlon_degrees
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_lat_varname
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_lon_varname
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_thickness_varname
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_draft_varname
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_ice_frac_varname
      character (len=StrKIND) :: config_global_ocean_land_ice_topo_grounded_frac_varname
      logical :: config_global_ocean_use_constant_land_ice_cavity_temperature
      real (kind=RKIND) :: config_global_ocean_constant_land_ice_cavity_temperature
      logical :: config_global_ocean_cull_inland_seas
      character (len=StrKIND) :: config_global_ocean_windstress_file
      character (len=StrKIND) :: config_global_ocean_windstress_nlat_dimname
      character (len=StrKIND) :: config_global_ocean_windstress_nlon_dimname
      logical :: config_global_ocean_windstress_latlon_degrees
      character (len=StrKIND) :: config_global_ocean_windstress_lat_varname
      character (len=StrKIND) :: config_global_ocean_windstress_lon_varname
      character (len=StrKIND) :: config_global_ocean_windstress_zonal_varname
      character (len=StrKIND) :: config_global_ocean_windstress_meridional_varname
      character (len=StrKIND) :: config_global_ocean_windstress_method
      real (kind=RKIND) :: config_global_ocean_windstress_conversion_factor
      character (len=StrKIND) :: config_global_ocean_ecosys_file
      character (len=StrKIND) :: config_global_ocean_ecosys_forcing_file
      character (len=StrKIND) :: config_global_ocean_ecosys_nlat_dimname
      character (len=StrKIND) :: config_global_ocean_ecosys_nlon_dimname
      character (len=StrKIND) :: config_global_ocean_ecosys_ndepth_dimname
      real (kind=RKIND) :: config_global_ocean_ecosys_depth_conversion_factor
      integer :: config_global_ocean_ecosys_vert_levels
      character (len=StrKIND) :: config_global_ocean_ecosys_lat_varname
      character (len=StrKIND) :: config_global_ocean_ecosys_lon_varname
      character (len=StrKIND) :: config_global_ocean_ecosys_depth_varname
      logical :: config_global_ocean_ecosys_latlon_degrees
      character (len=StrKIND) :: config_global_ocean_ecosys_method
      character (len=StrKIND) :: config_global_ocean_ecosys_forcing_time_dimname
      integer :: config_global_ocean_smooth_ecosys_iterations

      integer :: config_cvmix_WSwSBF_vert_levels
      real (kind=RKIND) :: config_cvmix_WSwSBF_surface_temperature
      real (kind=RKIND) :: config_cvmix_WSwSBF_surface_salinity
      real (kind=RKIND) :: config_cvmix_WSwSBF_surface_restoring_temperature
      real (kind=RKIND) :: config_cvmix_WSwSBF_surface_restoring_salinity
      real (kind=RKIND) :: config_cvmix_WSwSBF_temperature_piston_velocity
      real (kind=RKIND) :: config_cvmix_WSwSBF_salinity_piston_velocity
      real (kind=RKIND) :: config_cvmix_WSwSBF_sensible_heat_flux
      real (kind=RKIND) :: config_cvmix_WSwSBF_latent_heat_flux
      real (kind=RKIND) :: config_cvmix_WSwSBF_shortwave_heat_flux
      real (kind=RKIND) :: config_cvmix_WSwSBF_rain_flux
      real (kind=RKIND) :: config_cvmix_WSwSBF_evaporation_flux
      real (kind=RKIND) :: config_cvmix_WSwSBF_interior_temperature_restoring_rate
      real (kind=RKIND) :: config_cvmix_WSwSBF_interior_salinity_restoring_rate
      real (kind=RKIND) :: config_cvmix_WSwSBF_temperature_gradient
      real (kind=RKIND) :: config_cvmix_WSwSBF_salinity_gradient
      real (kind=RKIND) :: config_cvmix_WSwSBF_temperature_gradient_mixed_layer
      real (kind=RKIND) :: config_cvmix_WSwSBF_salinity_gradient_mixed_layer
      real (kind=RKIND) :: config_cvmix_WSwSBF_mixed_layer_depth_temperature
      real (kind=RKIND) :: config_cvmix_WSwSBF_mixed_layer_depth_salinity
      real (kind=RKIND) :: config_cvmix_WSwSBF_mixed_layer_temperature_change
      real (kind=RKIND) :: config_cvmix_WSwSBF_mixed_layer_salinity_change
      character (len=StrKIND) :: config_cvmix_WSwSBF_vertical_grid
      real (kind=RKIND) :: config_cvmix_WSwSBF_bottom_depth
      real (kind=RKIND) :: config_cvmix_WSwSBF_max_windstress
      real (kind=RKIND) :: config_cvmix_WSwSBF_coriolis_parameter

      integer :: config_iso_vert_levels
      real (kind=RKIND) :: config_iso_main_channel_depth
      real (kind=RKIND) :: config_iso_north_wall_lat
      real (kind=RKIND) :: config_iso_south_wall_lat
      logical :: config_iso_ridge_flag
      real (kind=RKIND) :: config_iso_ridge_center_lon
      real (kind=RKIND) :: config_iso_ridge_height
      real (kind=RKIND) :: config_iso_ridge_width
      logical :: config_iso_plateau_flag
      real (kind=RKIND) :: config_iso_plateau_center_lon
      real (kind=RKIND) :: config_iso_plateau_center_lat
      real (kind=RKIND) :: config_iso_plateau_height
      real (kind=RKIND) :: config_iso_plateau_radius
      real (kind=RKIND) :: config_iso_plateau_slope_width
      logical :: config_iso_shelf_flag
      real (kind=RKIND) :: config_iso_shelf_depth
      real (kind=RKIND) :: config_iso_shelf_width
      logical :: config_iso_cont_slope_flag
      real (kind=RKIND) :: config_iso_max_cont_slope
      logical :: config_iso_embayment_flag
      real (kind=RKIND) :: config_iso_embayment_center_lon
      real (kind=RKIND) :: config_iso_embayment_center_lat
      real (kind=RKIND) :: config_iso_embayment_radius
      real (kind=RKIND) :: config_iso_embayment_depth
      logical :: config_iso_depression_flag
      real (kind=RKIND) :: config_iso_depression_center_lon
      real (kind=RKIND) :: config_iso_depression_south_lat
      real (kind=RKIND) :: config_iso_depression_north_lat
      real (kind=RKIND) :: config_iso_depression_width
      real (kind=RKIND) :: config_iso_depression_depth
      real (kind=RKIND) :: config_iso_salinity
      real (kind=RKIND) :: config_iso_wind_stress_max
      real (kind=RKIND) :: config_iso_acc_wind
      real (kind=RKIND) :: config_iso_asf_wind
      real (kind=RKIND) :: config_iso_wind_trans
      real (kind=RKIND) :: config_iso_heat_flux_south
      real (kind=RKIND) :: config_iso_heat_flux_middle
      real (kind=RKIND) :: config_iso_heat_flux_north
      real (kind=RKIND) :: config_iso_heat_flux_lat_ss
      real (kind=RKIND) :: config_iso_heat_flux_lat_sm
      real (kind=RKIND) :: config_iso_heat_flux_lat_mn
      real (kind=RKIND) :: config_iso_region1_center_lon
      real (kind=RKIND) :: config_iso_region1_center_lat
      real (kind=RKIND) :: config_iso_region2_center_lon
      real (kind=RKIND) :: config_iso_region2_center_lat
      real (kind=RKIND) :: config_iso_region3_center_lon
      real (kind=RKIND) :: config_iso_region3_center_lat
      real (kind=RKIND) :: config_iso_region4_center_lon
      real (kind=RKIND) :: config_iso_region4_center_lat
      logical :: config_iso_heat_flux_region1_flag
      real (kind=RKIND) :: config_iso_heat_flux_region1
      real (kind=RKIND) :: config_iso_heat_flux_region1_radius
      logical :: config_iso_heat_flux_region2_flag
      real (kind=RKIND) :: config_iso_heat_flux_region2
      real (kind=RKIND) :: config_iso_heat_flux_region2_radius
      real (kind=RKIND) :: config_iso_surface_temperature_piston_velocity
      real (kind=RKIND) :: config_iso_initial_temp_t1
      real (kind=RKIND) :: config_iso_initial_temp_t2
      real (kind=RKIND) :: config_iso_initial_temp_h0
      real (kind=RKIND) :: config_iso_initial_temp_h1
      real (kind=RKIND) :: config_iso_initial_temp_mt
      real (kind=RKIND) :: config_iso_initial_temp_latS
      real (kind=RKIND) :: config_iso_initial_temp_latN
      real (kind=RKIND) :: config_iso_temperature_sponge_t1
      real (kind=RKIND) :: config_iso_temperature_sponge_h1
      real (kind=RKIND) :: config_iso_temperature_sponge_l1
      real (kind=RKIND) :: config_iso_temperature_sponge_tau1
      logical :: config_iso_temperature_restore_region1_flag
      real (kind=RKIND) :: config_iso_temperature_restore_t1
      real (kind=RKIND) :: config_iso_temperature_restore_lcx1
      real (kind=RKIND) :: config_iso_temperature_restore_lcy1
      logical :: config_iso_temperature_restore_region2_flag
      real (kind=RKIND) :: config_iso_temperature_restore_t2
      real (kind=RKIND) :: config_iso_temperature_restore_lcx2
      real (kind=RKIND) :: config_iso_temperature_restore_lcy2
      logical :: config_iso_temperature_restore_region3_flag
      real (kind=RKIND) :: config_iso_temperature_restore_t3
      real (kind=RKIND) :: config_iso_temperature_restore_lcx3
      real (kind=RKIND) :: config_iso_temperature_restore_lcy3
      logical :: config_iso_temperature_restore_region4_flag
      real (kind=RKIND) :: config_iso_temperature_restore_t4
      real (kind=RKIND) :: config_iso_temperature_restore_lcx4
      real (kind=RKIND) :: config_iso_temperature_restore_lcy4

      integer :: config_soma_vert_levels
      real (kind=RKIND) :: config_soma_domain_width
      real (kind=RKIND) :: config_soma_center_latitude
      real (kind=RKIND) :: config_soma_center_longitude
      real (kind=RKIND) :: config_soma_phi
      real (kind=RKIND) :: config_soma_bottom_depth
      real (kind=RKIND) :: config_soma_shelf_width
      real (kind=RKIND) :: config_soma_shelf_depth
      real (kind=RKIND) :: config_soma_ref_density
      real (kind=RKIND) :: config_soma_density_difference
      real (kind=RKIND) :: config_soma_thermocline_depth
      real (kind=RKIND) :: config_soma_density_difference_linear
      real (kind=RKIND) :: config_soma_surface_temperature
      real (kind=RKIND) :: config_soma_surface_salinity
      logical :: config_soma_use_surface_temp_restoring
      real (kind=RKIND) :: config_soma_surface_temp_restoring_at_center_latitude
      real (kind=RKIND) :: config_soma_surface_temp_restoring_latitude_gradient
      real (kind=RKIND) :: config_soma_restoring_temp_piston_vel

      integer :: config_ziso_vert_levels
      logical :: config_ziso_add_easterly_wind_stress_ASF
      real (kind=RKIND) :: config_ziso_wind_transition_position
      real (kind=RKIND) :: config_ziso_antarctic_shelf_front_width
      real (kind=RKIND) :: config_ziso_wind_stress_shelf_front_max
      logical :: config_ziso_use_slopping_bathymetry
      real (kind=RKIND) :: config_ziso_meridional_extent
      real (kind=RKIND) :: config_ziso_zonal_extent
      real (kind=RKIND) :: config_ziso_bottom_depth
      real (kind=RKIND) :: config_ziso_shelf_depth
      real (kind=RKIND) :: config_ziso_slope_half_width
      real (kind=RKIND) :: config_ziso_slope_center_position
      real (kind=RKIND) :: config_ziso_reference_coriolis
      real (kind=RKIND) :: config_ziso_coriolis_gradient
      real (kind=RKIND) :: config_ziso_wind_stress_max
      real (kind=RKIND) :: config_ziso_mean_restoring_temp
      real (kind=RKIND) :: config_ziso_restoring_temp_dev_ta
      real (kind=RKIND) :: config_ziso_restoring_temp_dev_tb
      real (kind=RKIND) :: config_ziso_restoring_temp_tau
      real (kind=RKIND) :: config_ziso_restoring_temp_piston_vel
      real (kind=RKIND) :: config_ziso_restoring_temp_ze
      real (kind=RKIND) :: config_ziso_restoring_sponge_l
      real (kind=RKIND) :: config_ziso_initial_temp_t1
      real (kind=RKIND) :: config_ziso_initial_temp_t2
      real (kind=RKIND) :: config_ziso_initial_temp_h1
      real (kind=RKIND) :: config_ziso_initial_temp_mt
      logical :: config_ziso_frazil_enable
      real (kind=RKIND) :: config_ziso_frazil_temperature_anomaly

      integer :: config_sub_ice_shelf_2D_vert_levels
      real (kind=RKIND) :: config_sub_ice_shelf_2D_bottom_depth
      real (kind=RKIND) :: config_sub_ice_shelf_2D_cavity_thickness
      real (kind=RKIND) :: config_sub_ice_shelf_2D_slope_height
      real (kind=RKIND) :: config_sub_ice_shelf_2D_edge_width
      real (kind=RKIND) :: config_sub_ice_shelf_2D_y1
      real (kind=RKIND) :: config_sub_ice_shelf_2D_y2
      real (kind=RKIND) :: config_sub_ice_shelf_2D_temperature
      real (kind=RKIND) :: config_sub_ice_shelf_2D_surface_salinity
      real (kind=RKIND) :: config_sub_ice_shelf_2D_bottom_salinity

      integer :: config_periodic_planar_vert_levels
      real (kind=RKIND) :: config_periodic_planar_bottom_depth
      real (kind=RKIND) :: config_periodic_planar_velocity_strength

      integer :: config_ecosys_column_vert_levels
      character (len=StrKIND) :: config_ecosys_column_vertical_grid
      character (len=StrKIND) :: config_ecosys_column_TS_filename
      character (len=StrKIND) :: config_ecosys_column_ecosys_filename
      real (kind=RKIND) :: config_ecosys_column_bottom_depth

      integer :: config_sea_mount_vert_levels
      character (len=StrKIND) :: config_sea_mount_layer_type
      character (len=StrKIND) :: config_sea_mount_stratification_type
      real (kind=RKIND) :: config_sea_mount_density_coef_linear
      real (kind=RKIND) :: config_sea_mount_density_coef_exp
      real (kind=RKIND) :: config_sea_mount_density_gradient_linear
      real (kind=RKIND) :: config_sea_mount_density_gradient_exp
      real (kind=RKIND) :: config_sea_mount_density_depth_linear
      real (kind=RKIND) :: config_sea_mount_density_depth_exp
      real (kind=RKIND) :: config_sea_mount_density_ref
      real (kind=RKIND) :: config_sea_mount_density_Tref
      real (kind=RKIND) :: config_sea_mount_density_alpha
      real (kind=RKIND) :: config_sea_mount_bottom_depth
      real (kind=RKIND) :: config_sea_mount_height
      real (kind=RKIND) :: config_sea_mount_radius
      real (kind=RKIND) :: config_sea_mount_width
      real (kind=RKIND) :: config_sea_mount_salinity
      real (kind=RKIND) :: config_sea_mount_coriolis_parameter

      integer :: config_isomip_vert_levels
      character (len=StrKIND) :: config_isomip_vertical_level_distribution
      real (kind=RKIND) :: config_isomip_bottom_depth
      real (kind=RKIND) :: config_isomip_temperature
      real (kind=RKIND) :: config_isomip_salinity
      real (kind=RKIND) :: config_isomip_restoring_temperature
      real (kind=RKIND) :: config_isomip_temperature_piston_velocity
      real (kind=RKIND) :: config_isomip_restoring_salinity
      real (kind=RKIND) :: config_isomip_salinity_piston_velocity
      real (kind=RKIND) :: config_isomip_coriolis_parameter
      real (kind=RKIND) :: config_isomip_southern_boundary
      real (kind=RKIND) :: config_isomip_northern_boundary
      real (kind=RKIND) :: config_isomip_western_boundary
      real (kind=RKIND) :: config_isomip_eastern_boundary
      real (kind=RKIND) :: config_isomip_y1
      real (kind=RKIND) :: config_isomip_z1
      real (kind=RKIND) :: config_isomip_ice_fraction1
      real (kind=RKIND) :: config_isomip_y2
      real (kind=RKIND) :: config_isomip_z2
      real (kind=RKIND) :: config_isomip_ice_fraction2
      real (kind=RKIND) :: config_isomip_y3
      real (kind=RKIND) :: config_isomip_z3
      real (kind=RKIND) :: config_isomip_ice_fraction3

      integer :: config_isomip_plus_vert_levels
      character (len=StrKIND) :: config_isomip_plus_vertical_level_distribution
      real (kind=RKIND) :: config_isomip_plus_max_bottom_depth
      integer :: config_isomip_plus_minimum_levels
      real (kind=RKIND) :: config_isomip_plus_min_column_thickness
      real (kind=RKIND) :: config_isomip_plus_min_ocean_fraction
      character (len=StrKIND) :: config_isomip_plus_topography_file
      real (kind=RKIND) :: config_isomip_plus_init_top_temp
      real (kind=RKIND) :: config_isomip_plus_init_bot_temp
      real (kind=RKIND) :: config_isomip_plus_init_top_sal
      real (kind=RKIND) :: config_isomip_plus_init_bot_sal
      real (kind=RKIND) :: config_isomip_plus_restore_top_temp
      real (kind=RKIND) :: config_isomip_plus_restore_bot_temp
      real (kind=RKIND) :: config_isomip_plus_restore_top_sal
      real (kind=RKIND) :: config_isomip_plus_restore_bot_sal
      real (kind=RKIND) :: config_isomip_plus_restore_rate
      real (kind=RKIND) :: config_isomip_plus_restore_evap_rate
      real (kind=RKIND) :: config_isomip_plus_restore_xMin
      real (kind=RKIND) :: config_isomip_plus_restore_xMax
      real (kind=RKIND) :: config_isomip_plus_coriolis_parameter
      real (kind=RKIND) :: config_isomip_plus_effective_density

      integer :: config_hurricane_vert_levels
      real (kind=RKIND) :: config_hurricane_min_depth
      real (kind=RKIND) :: config_hurricane_max_depth
      real (kind=RKIND) :: config_hurricane_gaussian_hump_amplitude
      logical :: config_hurricane_use_gaussian_hump
      real (kind=RKIND) :: config_hurricane_gaussian_lon_center
      real (kind=RKIND) :: config_hurricane_gaussian_lat_center
      real (kind=RKIND) :: config_hurricane_gaussian_width
      real (kind=RKIND) :: config_hurricane_gaussian_amplitude
      real (kind=RKIND) :: config_hurricane_gaussian_slr_amp
      real (kind=RKIND) :: config_hurricane_land_z_limit
      real (kind=RKIND) :: config_hurricane_marsh_z_limit
      real (kind=RKIND) :: config_hurricane_land_drag
      real (kind=RKIND) :: config_hurricane_marsh_drag
      real (kind=RKIND) :: config_hurricane_channel_drag
      real (kind=RKIND) :: config_hurricane_sea_level_rise_adjustment

      integer :: config_tidal_boundary_vert_levels
      integer :: config_tidal_boundary_min_vert_levels
      character (len=StrKIND) :: config_tidal_boundary_layer_type
      real (kind=RKIND) :: config_tidal_boundary_right_bottom_depth
      logical :: config_tidal_start_dry
      logical :: config_tidal_boundary_use_distances
      real (kind=RKIND) :: config_tidal_boundary_left_value
      real (kind=RKIND) :: config_tidal_boundary_right_value
      real (kind=RKIND) :: config_tidal_boundary_left_bottom_depth
      real (kind=RKIND) :: config_tidal_boundary_salinity
      real (kind=RKIND) :: config_tidal_boundary_domain_temperature
      real (kind=RKIND) :: config_tidal_boundary_plug_temperature
      real (kind=RKIND) :: config_tidal_boundary_plug_width_frac
      real (kind=RKIND) :: config_tidal_forcing_left_Cd_or_n
      real (kind=RKIND) :: config_tidal_forcing_right_Cd_or_n
      logical :: config_use_idealized_transect
      real (kind=RKIND) :: config_idealized_transect_Lshore
      real (kind=RKIND) :: config_idealized_transect_Sshore
      real (kind=RKIND) :: config_idealized_transect_Lcoast
      real (kind=RKIND) :: config_idealized_transect_Scoast
      real (kind=RKIND) :: config_idealized_transect_Lmarsh
      real (kind=RKIND) :: config_idealized_transect_Smarsh
      real (kind=RKIND) :: config_idealized_transect_roughness
      real (kind=RKIND) :: config_idealized_transect_roughness_marsh
      real (kind=RKIND) :: config_idealized_vegetation_diameter
      real (kind=RKIND) :: config_idealized_vegetation_height
      real (kind=RKIND) :: config_idealized_vegetation_density

      real (kind=RKIND) :: config_cosine_bell_temperature
      real (kind=RKIND) :: config_cosine_bell_salinity
      real (kind=RKIND) :: config_cosine_bell_lat_center
      real (kind=RKIND) :: config_cosine_bell_lon_center
      real (kind=RKIND) :: config_cosine_bell_psi0
      real (kind=RKIND) :: config_cosine_bell_radius
      real (kind=RKIND) :: config_cosine_bell_vel_pd

      integer :: config_mixed_layer_eddy_vert_levels
      real (kind=RKIND) :: config_mixed_layer_eddy_bottom_depth
      real (kind=RKIND) :: config_mixed_layer_eddy_mixed_layer_depth
      real (kind=RKIND) :: config_mixed_layer_eddy_base_temperature
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_stratification_mixed_layer
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_stratification_interior
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_horizontal_gradient
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_front_width
      real (kind=RKIND) :: config_mixed_layer_eddy_temperature_perturbation_magnitude
      real (kind=RKIND) :: config_mixed_layer_eddy_salinity
      logical :: config_mixed_layer_eddy_two_fronts
      real (kind=RKIND) :: config_mixed_layer_eddy_restoring_width
      real (kind=RKIND) :: config_mixed_layer_eddy_restoring_tau
      real (kind=RKIND) :: config_mixed_layer_eddy_heat_flux
      real (kind=RKIND) :: config_mixed_layer_eddy_evaporation_flux
      real (kind=RKIND) :: config_mixed_layer_eddy_wind_stress_zonal
      real (kind=RKIND) :: config_mixed_layer_eddy_wind_stress_meridional
      real (kind=RKIND) :: config_mixed_layer_eddy_coriolis_parameter

      integer :: config_test_sht_function_option
      real (kind=RKIND) :: config_test_sht_cosine_bell_lat_center
      real (kind=RKIND) :: config_test_sht_cosine_bell_lon_center
      real (kind=RKIND) :: config_test_sht_cosine_bell_psi0
      real (kind=RKIND) :: config_test_sht_cosine_bell_radius
      real (kind=RKIND) :: config_test_sht_function3_cell_width_equator
      real (kind=RKIND) :: config_test_sht_function3_cell_width_pole
      real (kind=RKIND) :: config_test_sht_function3_lat_transition_start
      real (kind=RKIND) :: config_test_sht_function3_lat_transition_width
      integer :: config_test_sht_n_iterations

      integer :: config_parabolic_bowl_vert_levels
      real (kind=RKIND) :: config_parabolic_bowl_Coriolis_parameter
      real (kind=RKIND) :: config_parabolic_bowl_eta0
      real (kind=RKIND) :: config_parabolic_bowl_b0
      real (kind=RKIND) :: config_parabolic_bowl_omega
      real (kind=RKIND) :: config_parabolic_bowl_gravity
      logical :: config_parabolic_bowl_adjust_domain_center

      logical :: config_alter_ICs_for_pcs
      character (len=StrKIND) :: config_pc_alteration_type
      real (kind=RKIND) :: config_min_pc_fraction

      character (len=StrKIND) :: config_init_configuration
      logical :: config_expand_sphere
      logical :: config_realistic_coriolis_parameter
      logical :: config_write_cull_cell_mask
      character (len=StrKIND) :: config_vertical_grid

      real (kind=RKIND) :: config_1dCVTgenerator_stretch1
      real (kind=RKIND) :: config_1dCVTgenerator_stretch2
      real (kind=RKIND) :: config_1dCVTgenerator_dzSeed

      character (len=StrKIND) :: config_init_vertical_grid_type

      integer :: config_rx1_outer_iter_count
      integer :: config_rx1_inner_iter_count
      real (kind=RKIND) :: config_rx1_init_inner_weight
      real (kind=RKIND) :: config_rx1_max
      real (kind=RKIND) :: config_rx1_horiz_smooth_weight
      real (kind=RKIND) :: config_rx1_vert_smooth_weight
      real (kind=RKIND) :: config_rx1_slope_weight
      real (kind=RKIND) :: config_rx1_zstar_weight
      integer :: config_rx1_horiz_smooth_open_ocean_cells
      integer :: config_rx1_min_levels
      real (kind=RKIND) :: config_rx1_min_layer_thickness

