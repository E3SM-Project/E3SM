
EAMxx runtime configurable parameters
=====================================

# Atmosphere Processes Parameters

## sc_import

* sc_import::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* sc_import::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* sc_import::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* sc_import::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* sc_import::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* sc_import::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## sc_export

* sc_export::prescribed_constants::fields:
    - description: **MISSING**
    - type: array(string)
* sc_export::prescribed_constants::values:
    - description: **MISSING**
    - type: array(real)

* sc_export::prescribed_from_file::fields:
    - description: **MISSING**
    - type: array(string)
* sc_export::prescribed_from_file::files:
    - description: **MISSING**
    - type: array(string)

* sc_export::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* sc_export::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* sc_export::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* sc_export::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* sc_export::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* sc_export::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## homme

* homme::Moisture:
    - description: **MISSING**
    - type: **MISSING**
* homme::BfbHash:
    - description: **MISSING**
    - type: integer
* homme::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* homme::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* homme::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* homme::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* homme::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* homme::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## p3

* p3::do_prescribed_ccn:
    - description: **MISSING**
    - type: **MISSING**
* p3::do_predict_nc:
    - description: **MISSING**
    - type: **MISSING**
* p3::enable_column_conservation_checks:
    - description: **MISSING**
    - type: **MISSING**
* p3::tables:
    - description: **MISSING**
    - type: array(file)
* p3::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* p3::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* p3::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* p3::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* p3::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* p3::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## shoc

* shoc::enable_column_conservation_checks:
    - description: **MISSING**
    - type: **MISSING**
* shoc::check_flux_state_consistency:
    - description: **MISSING**
    - type: **MISSING**
* shoc::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* shoc::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* shoc::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* shoc::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* shoc::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* shoc::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## cldFraction

* cldFraction::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* cldFraction::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* cldFraction::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* cldFraction::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* cldFraction::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* cldFraction::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## testOnly

* testOnly::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* testOnly::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* testOnly::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* testOnly::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* testOnly::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* testOnly::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## spa

* spa::spa_remap_file:
    - description: **MISSING**
    - type: file
* spa::spa_data_file:
    - description: **MISSING**
    - type: file
* spa::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* spa::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* spa::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* spa::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* spa::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* spa::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## rrtmgp

* rrtmgp::rrtmgp_coefficients_file_sw:
    - description: **MISSING**
    - type: file
* rrtmgp::rrtmgp_coefficients_file_lw:
    - description: **MISSING**
    - type: file
* rrtmgp::rrtmgp_cloud_optics_file_sw:
    - description: **MISSING**
    - type: file
* rrtmgp::rrtmgp_cloud_optics_file_lw:
    - description: **MISSING**
    - type: file
* rrtmgp::column_chunk_size:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::active_gases:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::ch4vmr:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::co2vmr:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::n2ovmr:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::f11vmr:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::f12vmr:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::n2vmr:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::covmr:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::orbital_year:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::orbital_eccentricity:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::orbital_obliquity:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::orbital_mvelp:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::rad_frequency:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::do_aerosol_rad:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::enable_column_conservation_checks:
    - description: **MISSING**
    - type: **MISSING**
* rrtmgp::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* rrtmgp::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* rrtmgp::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* rrtmgp::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* rrtmgp::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* rrtmgp::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## mac_aero_mic

* mac_aero_mic::atm_procs_list:
    - description: **MISSING**
    - type: **MISSING**
* mac_aero_mic::number_of_subcycles:
    - description: **MISSING**
    - type: **MISSING**
* mac_aero_mic::Type:
    - description: **MISSING**
    - type: **MISSING**
* mac_aero_mic::schedule_type:
    - description: **MISSING**
    - type: **MISSING**
    - valid values: Sequential
* mac_aero_mic::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* mac_aero_mic::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* mac_aero_mic::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* mac_aero_mic::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* mac_aero_mic::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

## physics

* physics::atm_procs_list:
    - description: **MISSING**
    - type: **MISSING**
* physics::Type:
    - description: **MISSING**
    - type: **MISSING**
* physics::schedule_type:
    - description: **MISSING**
    - type: **MISSING**
    - valid values: Sequential
* physics::number_of_subcycles:
    - description: how many times to subcycle this atm process
    - type: **MISSING**
    - constraints: gt 0
* physics::enable_precondition_checks:
    - description: **MISSING**
    - type: logical
* physics::enable_postcondition_checks:
    - description: **MISSING**
    - type: logical
* physics::repair_log_level:
    - description: **MISSING**
    - type: string
    - valid values: trace,debug,info,warn
* physics::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* physics::compute_tendencies:
    - description: list of computed fields for which this process will back out tendencies
    - type: array(string)

# Initial Conditions Parameters

* initial_conditions::Filename:
    - description: **MISSING**
    - type: file
* initial_conditions::topography_filename:
    - description: **MISSING**
    - type: file
* initial_conditions::phis:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::restart_casename:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::surf_evap:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::precip_liq_surf_mass:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::precip_ice_surf_mass:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::cldfrac_liq:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::sgs_buoy_flux:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::eddy_diff_mom:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::T_prev_micro_step:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::qv_prev_micro_step:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::qr:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::nr:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::qm:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::bm:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::ni_activated:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::nc_nuceat_tend:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::tke:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::sfc_alb_dir_vis:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::sfc_alb_dir_nir:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::sfc_alb_dif_vis:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::sfc_alb_dif_nir:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::surf_sens_flux:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::surf_lw_flux_up:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::surf_mom_flux:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::qc:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::qi:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::nc:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::ni:
    - description: **MISSING**
    - type: **MISSING**
* initial_conditions::o3_volume_mix_ratio:
    - description: **MISSING**
    - type: **MISSING**

# Atmosphere Driver Parameters

* driver_options::atmosphere_dag_verbosity_level:
    - description: **MISSING**
    - type: **MISSING**
* driver_options::atm_log_level:
    - description: **MISSING**
    - type: **MISSING**
    - valid values: trace,debug,info,warn,error
* driver_options::output_to_screen:
    - description: **MISSING**
    - type: logical
* driver_options::mass_column_conservation_error_tolerance:
    - description: **MISSING**
    - type: **MISSING**
* driver_options::energy_column_conservation_error_tolerance:
    - description: **MISSING**
    - type: **MISSING**
* driver_options::column_conservation_checks_fail_handling_type:
    - description: **MISSING**
    - type: **MISSING**
* driver_options::check_all_computed_fields_for_nans:
    - description: **MISSING**
    - type: logical

# Scorpio Parameters

* Scorpio::output_yaml_files:
    - description: **MISSING**
    - type: array(string)
* Scorpio::model_restart::filename_prefix:
    - description: **MISSING**
    - type: **MISSING**


# Homme namelist

* ctl_nl::cubed_sphere_map:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::disable_diagnostics:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::dt_remap_factor:
    - description: **MISSING**
    - type: **MISSING**
    - constraints: ge 1
* ctl_nl::dt_tracer_factor:
    - description: **MISSING**
    - type: **MISSING**
    - constraints: ge 1
* ctl_nl::hv_ref_profiles:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::hypervis_order:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::hypervis_scaling:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::hypervis_subcycle:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::hypervis_subcycle_tom:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::hypervis_subcycle_q:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::nu:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::nu_top:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::pgrad_correction:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::se_ftype:
    - description: **MISSING**
    - type: **MISSING**
    - valid values: 0,2
* ctl_nl::se_geometry:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::se_limiter_option:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::se_ne:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::se_ne_x:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::se_ne_y:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::se_nsplit:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::se_partmethod:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::se_topology:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::se_tstep:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::statefreq:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::theta_advect_form:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::theta_hydrostatic_mode:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::tstep_type:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::vert_remap_q_alg:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::transport_alg:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::vtheta_thresh:
    - description: **MISSING**
    - type: **MISSING**
* ctl_nl::internal_diagnostics_level:
    - description: **MISSING**
    - type: integer
* ctl_nl::mesh_file:
    - description: **MISSING**
    - type: file

