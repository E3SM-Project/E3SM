!--------------------------------------------------------------------------------------------------
! $Id$
!==================================================================================================
!
!       ########  ###       ###    ### #########  #########           ###     ######### ###########
!     ###    ### ###       ###    ### ###    ### ###    ###        ### ###   ###    ###    ###
!    ###        ###       ###    ### ###    ### ###    ###       ###   ###  ###    ###    ###
!   ###        ###       ###    ### #########  #########       ########### #########     ###
!  ###        ###       ###    ### ###    ### ###    ###      ###     ### ###           ###
! ###    ### ###       ###    ### ###    ### ###    ###      ###     ### ###           ###
! ########  ########## ########  #########  #########       ###     ### ###       ###########
!
! The CLUBB API serves as the doorway through which external models can interact with CLUBB.
!
!               PLEASE REMEMBER, IF ANY CODE IS CHANGED IN THIS DOCUMENT,
!                   THE CHANGES MUST BE PROPOGATED TO ALL HOST MODELS.
!
module clubb_api_module


  use mt95, only : &
    assignment( = ), &
    genrand_state, & ! Internal representation of the RNG state.
    genrand_srepr, & ! Public representation of the RNG state. Should be used to save the RNG state
    genrand_intg, &
    genrand_init_api => genrand_init

  use array_index, only : &
    hydromet_list, &
    hydromet_tol, & ! Tolerance values for all hydrometeors    [units vary]
    iiNg, & ! Hydrometeor array index for graupel concentration, Ng
    iiNi, & ! Hydrometeor array index for ice concentration, Ni
    iiNr, & ! Hydrometeor array index for rain drop concentration, Nr
    iiNs, & ! Hydrometeor array index for snow concentration, Ns
    iirg, & ! Hydrometeor array index for graupel mixing ratio, rg
    iiri, & ! Hydrometeor array index for ice mixing ratio, ri
    iirr, & ! Hydrometeor array index for rain water mixing ratio, rr
    iirs, & ! Hydrometeor array index for snow mixing ratio, rs
    iiPDF_chi, &
    iiPDF_rr,  &
    iiPDF_w,   &
    iiPDF_Nr,  &
    iiPDF_ri,  &
    iiPDF_Ni,  &
    iiPDF_Ncn, &
    iiPDF_rs,  &
    iiPDF_Ns,  &
    iiPDF_rg,  &
    iiPDF_Ng,  &
    iisclr_rt, &
    iisclr_thl, &
    iisclr_CO2, &
    iiedsclr_rt, &
    iiedsclr_thl, &
    iiedsclr_CO2, &
    l_frozen_hm, & ! if true, then the hydrometeor is frozen; otherwise liquid
    l_mix_rat_hm ! if true, then the quantity is a hydrometeor mixing ratio

  use clubb_precision, only : &
    time_precision, &
    core_rknd, &
    stat_nknd, &
    stat_rknd, &
    dp  ! Double Precision

  use constants_clubb, only : &
    cloud_frac_min, & ! Threshold for cloud fractions
    cm3_per_m3, & ! Cubic centimeters per cubic meter
    Cp, & ! Dry air specific heat at constant p [J/kg/K]
    em_min, & ! Minimum value for em (turbulence kinetic energy)
    ep, & ! ep  = 0.622  [-]
    fstderr, & ! Fortran file unit I/O constant
    fstdout, & ! Fortran file unit I/O constant
    grav, & ! Gravitational acceleration     [m/s^2]
    Ls, & ! Latent heat of sublimation          [J/kg]
    Lv, & ! Latent heat of vaporization         [J/kg]
    Lf, & ! Latent heat of fusion               [J/kg]
    pi, & ! The ratio of radii to their circumference
    pi_dp, & ! pi in double precision
    radians_per_deg_dp, &
    Rd, & ! Dry air gas constant                [J/kg/K]
    Rv, & ! Water vapor gas constant            [J/kg/K]
    sec_per_day, & ! Seconds in a day.
    sec_per_hr, &  ! Seconds in an hour.
    sec_per_min, & ! Seconds in a minute.
    T_freeze_K, & ! Freezing point of water [K]
    var_length, & ! Maximum variable name length in CLUBB GrADS or netCDF output
    zero, & ! 0.0_core_rknd
    zero_threshold, & ! Defining a threshold on a physical quantity to be 0.
    ! Tolerances
    Nc_tol, & ! Tolerance value for N_c  [#/kg]
    Ng_tol, & ! Tolerance value for N_s [#/kg]
    Ni_tol, & ! Tolerance value for N_i [#/kg]
    Nr_tol, & ! Tolerance value for N_r [#/kg]
    Ns_tol, & ! Tolerance value for N_s [#/kg]
    rg_tol, & ! Tolerance value for r_g [kg/kg]
    rho_lw, &
    ri_tol, & ! Tolerance value for r_i [kg/kg]
    rr_tol, & ! Tolerance value for r_r [kg/kg]
    rs_tol, & ! Tolerance value for r_s [kg/kg]
    rt_tol, & ! [kg/kg]
    thl_tol, & ! [K]
    w_tol_sqd ! [m^2/s^2]

  use corr_varnce_module, only : &
      corr_array_n_cloud, & ! Variable(s)
      corr_array_n_below, &
      pdf_dim,        &
      hmp2_ip_on_hmm2_ip, &
      Ncnp2_on_Ncnm2,     &
      hmp2_ip_on_hmm2_ip_slope_type,      & ! Types
      hmp2_ip_on_hmm2_ip_intrcpt_type

  use error_code, only: &
      clubb_at_least_debug_level,  & ! Procedure
      err_code,                    & ! Error Indicator
      clubb_no_error,              & ! Constants
      clubb_fatal_error

  use hydromet_pdf_parameter_module, only : &
    hydromet_pdf_parameter, &
    precipitation_fractions

  use model_flags, only : &
      clubb_config_flags_type, & ! Type
      iiPDF_ADG1, &
      iiPDF_new_hybrid, &
      ipdf_pre_advance_fields, &
      ipdf_post_advance_fields, &
      l_use_boussinesq    ! Use Boussinesq form of predictive equations (default is Anelastic).

  use parameters_model, only : &
    hydromet_dim    ! Number of hydrometeor species

  use parameters_tunable, only : &
    params_list,         & ! Variable(s)
    nu_vertical_res_dep    ! Type(s)

  use parameter_indices, only:  &
    nparams, & ! Variable(s)
    iC1, iC1b, iC1c, &
    iC2rt, iC2thl, iC2rtthl, iC4, iC_uu_shr, iC_uu_buoy, &
    iC6rt, iC6rtb, iC6rtc, iC6thl, iC6thlb, iC6thlc, &
    iC7, iC7b, iC7c, iC8, iC8b, iC10, iC11, iC11b, iC11c, &
    iC12, iC13, iC14, iC_wp3_pr_turb, iC_wp3_pr_dfsn, iC_wp2_splat, &
    iC6rt_Lscale0, iC6thl_Lscale0, &
    iC7_Lscale0, iwpxp_L_thresh, ic_K, ic_K1, inu1, &
    ic_K2, inu2, ic_K6, inu6, ic_K8, inu8, ic_K9, inu9, &
    inu10, ic_K_hm, ic_K_hmb, iK_hm_min_coef, inu_hm, &
    islope_coef_spread_DG_means_w, ipdf_component_stdev_factor_w, &
    icoef_spread_DG_means_rt, icoef_spread_DG_means_thl, &
    ibeta, igamma_coef, igamma_coefb, igamma_coefc, ilmin_coef, &
    iomicron, izeta_vrnce_rat, iupsilon_precip_frac_rat, &
    ilambda0_stability_coef, imult_coef, itaumin, itaumax, imu, &
    iLscale_mu_coef, iLscale_pert_coef, ialpha_corr, iSkw_denom_coef, &
    ic_K10, ic_K10h, ithlp2_rad_coef, ithlp2_rad_cloud_frac_thresh, &
    iup2_sfc_coef, iSkw_max_mag, iC_invrs_tau_bkgnd, &
    iC_invrs_tau_sfc, iC_invrs_tau_shear, iC_invrs_tau_N2, &
    iC_invrs_tau_N2_wp2, iC_invrs_tau_N2_xp2, iC_invrs_tau_N2_wpxp, &
    iC_invrs_tau_N2_clear_wp3, ialtitude_threshold, irtp2_clip_coef, &
    iRichardson_num_min, iRichardson_num_max, ia3_coef_min, iCx_min, iCx_max


  use pdf_parameter_module, only : &
! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
    num_pdf_params, &
!#endif
    pdf_parameter, &
    implicit_coefs_terms

  use sponge_layer_damping, only : &
    thlm_sponge_damp_settings,    & ! Variable(s)
    rtm_sponge_damp_settings,     &
    uv_sponge_damp_settings,      &
    wp2_sponge_damp_settings,     &
    wp3_sponge_damp_settings,     &
    up2_vp2_sponge_damp_settings, &
    thlm_sponge_damp_profile,     &
    rtm_sponge_damp_profile,      &
    uv_sponge_damp_profile,       &
    wp2_sponge_damp_profile,      &
    wp3_sponge_damp_profile,      &
    up2_vp2_sponge_damp_profile

  use stat_file_module, only : &
    clubb_i, &    ! Used to output multiple columns
    clubb_j       ! The indices must not exceed nlon (for i) or nlat (for j).

  use stats_rad_zm_module, only : &
    nvarmax_rad_zm ! Maximum variables allowed

  use stats_rad_zt_module, only : &
    nvarmax_rad_zt  ! Maximum variables allowed

  use stats_variables, only : &
    l_stats_last, & ! Last time step of output period
    stats_tsamp, & ! Sampling interval   [s]
    stats_tout, & ! Output interval     [s]
    l_output_rad_files, & ! Flag to turn off radiation statistics output
    l_stats, & ! Main flag to turn statistics on/off
    l_stats_samp, & ! Sample flag for current time step
    l_grads, & ! Output to GrADS format
    fname_rad_zt, & ! Name of the stats file for the stats_zt radiation grid fields
    fname_rad_zm, & ! Name of the stats file for the stats_zm radiation grid fields
    fname_sfc, & ! Name of the stats file for surface only fields
    l_netcdf, & ! Output to NetCDF format
    ! These are used in CAM only
    ztscr01, ztscr02, ztscr03, &
    ztscr04, ztscr05, ztscr06, &
    ztscr07, ztscr08, ztscr09, &
    ztscr10, ztscr11, ztscr12, &
    ztscr13, ztscr14, ztscr15, &
    ztscr16, ztscr17, ztscr18, &
    ztscr19, ztscr20, ztscr21, &
    zmscr01, zmscr02, zmscr03, &
    zmscr04, zmscr05, zmscr06, &
    zmscr07, zmscr08, zmscr09, &
    zmscr10, zmscr11, zmscr12, &
    zmscr13, zmscr14, zmscr15, &
    zmscr16, zmscr17

  use stats_zm_module, only : &
    nvarmax_zm ! Maximum variables allowed

  use stats_zt_module, only : &
    nvarmax_zt ! Maximum variables allowed

  use stats_sfc_module, only : &
    nvarmax_sfc

  use grid_class, only: grid ! Type

  use stats_type, only: stats ! Type

  implicit none

  private

  public &
    ! To Implement CLUBB:
    set_default_parameters_api, & ! Procedure(s)
    read_parameters_api, &
    setup_clubb_core_api, &
        ! CLUBB can be set more specifically using these flags:
        iiPDF_ADG1, &
        iiPDF_new_hybrid, &
        ipdf_pre_advance_fields, &
        ipdf_post_advance_fields, &
        l_use_boussinesq, &
        ! The parameters of CLUBB can be retrieved and tuned using these indices:
        iC1, iC1b, iC1c, &
        iC2rt, iC2thl, iC2rtthl, iC4, iC_uu_shr, iC_uu_buoy, &
        iC6rt, iC6rtb, iC6rtc, iC6thl, iC6thlb, iC6thlc, &
        iC7, iC7b, iC7c, iC8, iC8b, iC10, iC11, iC11b, iC11c, &
        iC12, iC13, iC14, iC_wp3_pr_turb, iC_wp3_pr_dfsn, iC_wp2_splat, & 
        iC6rt_Lscale0, iC6thl_Lscale0, &
        iC7_Lscale0, iwpxp_L_thresh, ic_K, ic_K1, inu1, &
        ic_K2, inu2, ic_K6, inu6, ic_K8, inu8, ic_K9, inu9, &
        inu10, ic_K_hm, ic_K_hmb, iK_hm_min_coef, inu_hm, &
        islope_coef_spread_DG_means_w, ipdf_component_stdev_factor_w, &
        icoef_spread_DG_means_rt, icoef_spread_DG_means_thl, &
        ibeta, igamma_coef, igamma_coefb, igamma_coefc, ilmin_coef, &
        iomicron, izeta_vrnce_rat, iupsilon_precip_frac_rat, &
        ilambda0_stability_coef, imult_coef, itaumin, itaumax, imu, &
        iLscale_mu_coef, iLscale_pert_coef, ialpha_corr, iSkw_denom_coef, &
        ic_K10, ic_K10h, ithlp2_rad_coef, ithlp2_rad_cloud_frac_thresh, &
        iup2_sfc_coef, iSkw_max_mag, iC_invrs_tau_bkgnd, &
        iC_invrs_tau_sfc, iC_invrs_tau_shear, iC_invrs_tau_N2, &
        iC_invrs_tau_N2_wp2, iC_invrs_tau_N2_xp2, iC_invrs_tau_N2_wpxp, &
        iC_invrs_tau_N2_clear_wp3, ialtitude_threshold, irtp2_clip_coef, &
        iRichardson_num_min, iRichardson_num_max, ia3_coef_min, iCx_min, iCx_max



  public &
        advance_clubb_core_api, &
        advance_clubb_core_api_single_col, &
        advance_clubb_core_api_multi_col, &
        pdf_parameter, &
        implicit_coefs_terms, &
        ! A hydromet array is required, and these variables are required for a hydromet array:
        hydromet_list, &
        hydromet_tol, &
        hydromet_dim, &
        iiNg, &
        iiNi, &
        iiNr, &
        iiNs, &
        iirg, &
        iiri, &
        iirr, &
        iirs, &
        iisclr_rt, &
        iisclr_thl, &
        iisclr_CO2, &
        iiedsclr_rt, &
        iiedsclr_thl, &
        iiedsclr_CO2, &
        l_frozen_hm, &
        l_mix_rat_hm, &
    cleanup_clubb_core_api

  public &
    ! To Implement SILHS:
    setup_corr_varnce_array_api, &
    setup_pdf_parameters_api, &
    hydromet_pdf_parameter, &
    init_pdf_hydromet_arrays_api, &
    ! generate_silhs_sample - SILHS API
    genrand_init_api, & ! if you are doing restarts)
    genrand_state, &
    genrand_srepr, &
    genrand_intg, &
    ! To use the results, you will need these variables:
    corr_array_n_cloud, &
    corr_array_n_below, &
    pdf_dim,        &
    iiPDF_chi,          &
    iiPDF_rr,           &
    iiPDF_w,            &
    iiPDF_Nr,           &
    iiPDF_ri,           &
    iiPDF_Ni,           &
    iiPDF_Ncn,          &
    iiPDF_rs,           &
    iiPDF_Ns,           &
    iiPDF_rg,           &
    iiPDF_Ng,           &
    hmp2_ip_on_hmm2_ip, &
    Ncnp2_on_Ncnm2,     &
    hmp2_ip_on_hmm2_ip_slope_type,      & ! Types
    hmp2_ip_on_hmm2_ip_intrcpt_type, &
    grid, &
    stats

  public &
    ! To Interact With CLUBB's Grid:
    ! For Varying Grids
    setup_grid_heights_api, &    ! if heights vary with time
    setup_grid_api

  public &
    ! To Obtain More Output from CLUBB for Diagnostics:
    stats_begin_timestep_api, &
    stats_end_timestep_api, &
    stats_finalize_api, &
    stats_init_api, &
    l_stats, &
    l_stats_last, &
    l_stats_samp, &
    stats_tsamp, &
    stats_tout

  public :: &
    calculate_thlp2_rad_api, params_list, &
    update_xp2_mc_api, sat_mixrat_liq_api

  public :: &
    ! To Convert Between Common CLUBB-related quantities:
    lin_interpolate_two_points_api, & ! OR
    lin_interpolate_on_grid_api, &
    T_in_K2thlm_api, &
    thlm2T_in_K_api, &
    zm2zt_api, &
    zt2zm_api

  public &
    ! To Check For and Handle CLUBB's Errors:
    calculate_spurious_source_api, &
    clubb_at_least_debug_level_api, &
    clubb_fatal_error, &
    clubb_no_error, &
    fill_holes_driver_api, & ! OR
    fill_holes_vertical_api, &
    fill_holes_hydromet_api, &
    set_clubb_debug_level_api, &
    vertical_integral_api

  public &
    ! Constants That May be Helpful:
    cloud_frac_min, &
    cm3_per_m3, &
    core_rknd, &
    Cp, &
    dp, &
    em_min, &
    ep, &
    fstderr, &
    fstdout, &
    grav, &
    Lf, &
    Ls, &
    Lv, &
    pi_dp, &
    pi, &
    radians_per_deg_dp, &
    Rd, &
    Rv, &
    sec_per_day, &
    sec_per_hr, &
    sec_per_min, &
    T_freeze_K, &
    time_precision, &
    var_length, &
    zero_threshold, &
    zero, &
    ! Tolerances
    Nc_tol, &
    Ng_tol, &
    Ni_tol, &
    Nr_tol, &
    Ns_tol, &
    rg_tol, &
    rho_lw, &
    ri_tol, &
    rr_tol, &
    rs_tol, &
    rt_tol, &
    thl_tol, &
    w_tol_sqd

  public &
    ! Attempt to Not Use the Following:
! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
    pack_pdf_params_api, &
    unpack_pdf_params_api, &
    num_pdf_params, &
!#endif
    init_pdf_params_api, &
    init_precip_fracs_api, &
    precipitation_fractions, &
    init_pdf_implicit_coefs_terms_api, &
    adj_low_res_nu_api, &
    assignment( = ), &
    clubb_i, &
    clubb_j, &
    compute_current_date_api, &
    fname_rad_zm, &
    fname_rad_zt, &
    fname_sfc, &
    gregorian2julian_day_api, &
    l_grads, &
    l_netcdf, &
    l_output_rad_files, &
    leap_year_api, &
    nvarmax_rad_zm, &
    nvarmax_rad_zt, &
    nvarmax_sfc, &
    nvarmax_zm, &
    nvarmax_zt
    public &
    nparams, &
    nu_vertical_res_dep, &
    setup_parameters_api, &
    stat_nknd, &
    stat_rknd, &
    stats_accumulate_hydromet_api, &
    stats_init_rad_zm_api, &
    stats_init_rad_zt_api, &
    stats_init_sfc_api, &
    stats_init_zm_api, &
    stats_init_zt_api, &
    zmscr01, zmscr02, zmscr03, &
    zmscr04, zmscr05, zmscr06, &
    zmscr07, zmscr08, zmscr09, &
    zmscr10, zmscr11, zmscr12, &
    zmscr13, zmscr14, zmscr15, &
    zmscr16, zmscr17, &
    ztscr01, ztscr02, ztscr03, &
    ztscr04, ztscr05, ztscr06, &
    ztscr07, ztscr08, ztscr09, &
    ztscr10, ztscr11, ztscr12, &
    ztscr13, ztscr14, ztscr15, &
    ztscr16, ztscr17, ztscr18, &
    ztscr19, ztscr20, ztscr21

  public &
    ! Needed to use the configurable CLUBB flags
    clubb_config_flags_type, &
    set_default_clubb_config_flags_api, &
    initialize_clubb_config_flags_type_api, &
    print_clubb_config_flags_api

  public &
    thlm_sponge_damp_settings,      & ! Variable(s)
    rtm_sponge_damp_settings,       &
    uv_sponge_damp_settings,        &
    wp2_sponge_damp_settings,       &
    wp3_sponge_damp_settings,       &
    up2_vp2_sponge_damp_settings,   &
    thlm_sponge_damp_profile,       &
    rtm_sponge_damp_profile,        &
    uv_sponge_damp_profile,         &
    wp2_sponge_damp_profile,        &
    wp3_sponge_damp_profile,        &
    up2_vp2_sponge_damp_profile,    &
    initialize_tau_sponge_damp_api, & ! Procedure(s)
    finalize_tau_sponge_damp_api
    
  public &
   copy_single_pdf_params_to_multi, &
   copy_multi_pdf_params_to_single

  interface zt2zm_api
    module procedure zt2zm_scalar_api, zt2zm_prof_api, zt2zm_2D_api
  end interface

  interface zm2zt_api
    module procedure zm2zt_scalar_api, zm2zt_prof_api, zm2zt_2D_api
  end interface
  
  interface setup_pdf_parameters_api
    module procedure setup_pdf_parameters_api_single_col
    module procedure setup_pdf_parameters_api_multi_col
  end interface
  
  interface advance_clubb_core_api
    module procedure advance_clubb_core_api_single_col
    module procedure advance_clubb_core_api_multi_col
  end interface

contains

  !================================================================================================
  ! advance_clubb_core - Advances the model one timestep.
  !================================================================================================

  subroutine advance_clubb_core_api_single_col( gr, &
    l_implemented, dt, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! intent(in)
    sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! intent(in)
    wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! intent(in)
    rtpthlp_forcing, wm_zm, wm_zt, &                        ! intent(in)
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &            ! intent(in)
    wpsclrp_sfc, wpedsclrp_sfc, &                           ! intent(in)
    upwp_sfc_pert, vpwp_sfc_pert, &                         ! intent(in)
    rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &            ! Intent(in)
    p_in_Pa, rho_zm, rho, exner, &                          ! intent(in)
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! intent(in)
    invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &      ! intent(in)
    rfrzm, radf, &                                          ! intent(in)
#ifdef CLUBBND_CAM
    varmu, &                                                ! intent(in)
#endif
    wphydrometp, wp2hmp, rtphmp, thlphmp, &                 ! intent(in)
    host_dx, host_dy, &                                     ! intent(in)
    clubb_params, nu_vert_res_dep, lmin, &                  ! intent(in)
    clubb_config_flags, &                                   ! intent(in)
    stats_zt, stats_zm, stats_sfc, &                        ! intent(inout)
    um, vm, upwp, vpwp, up2, vp2, up3, vp3, &               ! intent(inout)
    thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! intent(inout)
    sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16    ! intent(inout)
#endif
    sclrp2, sclrp3, sclrprtp, sclrpthlp, &                  ! intent(inout)
    wpsclrp, edsclrm, err_code_api, &                       ! intent(inout)
    rcm, cloud_frac, &                                      ! intent(inout)
    wpthvp, wp2thvp, rtpthvp, thlpthvp, &                   ! intent(inout)
    sclrpthvp, &                                            ! intent(inout)
    wp2rtp, wp2thlp, uprcp, vprcp, rc_coef, wp4, &          ! intent(inout)
    wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &      ! intent(inout)
    um_pert, vm_pert, upwp_pert, vpwp_pert, &               ! intent(inout)
    pdf_params, pdf_params_zm, &                            ! intent(inout)
    pdf_implicit_coefs_terms, &                             ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                  ! intent(inout)
               do_liquid_only_in_clubb, &                   ! intent(in)
#endif
    Kh_zm, Kh_zt, &                                         ! intent(out)
#ifdef CLUBB_CAM
    qclvar, &                                               ! intent(out)
#endif
    thlprcp, wprcp, w_up_in_cloud, &                        ! intent(out)
    rcm_in_layer, cloud_cover, invrs_tau_zm )               ! intent(out)

    use advance_clubb_core_module, only : advance_clubb_core

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type(s)

    use parameters_model, only: &
        sclr_dim, & ! Variable(s)
        edsclr_dim

    use model_flags, only: &
        clubb_config_flags_type
        
    use stats_zm_module, only: &
        stats_init_zm ! Procedure(s)
        
    use stats_zt_module, only: &
        stats_init_zt
        
    use stats_sfc_module, only: &
        stats_init_sfc ! Procedure(s)

    implicit none

    type(stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    type(grid), target, intent(in) :: gr
      !!! Input Variables
    logical, intent(in) ::  &
      l_implemented ! Is this part of a larger host model (T/F) ?

    real( kind = core_rknd ), intent(in) ::  &
      dt  ! Current timestep duration    [s]

    real( kind = core_rknd ), intent(in) ::  &
      fcor,  &          ! Coriolis forcing             [s^-1]
      sfc_elevation     ! Elevation of ground level    [m AMSL]

    integer, intent(in) :: &
      hydromet_dim      ! Total number of hydrometeors          [#]

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  &
      thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm,           & ! w mean wind component on momentum levels  [m/s]
      wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      rfrzm              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet           ! Collection of hydrometeors                [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      radf          ! Buoyancy production at the CL top due to LW radiative cooling [m^2/s^3]

#ifdef CLUBBND_CAM 
    real( kind = core_rknd ), intent(in) :: & 
      varmu 
#endif 

    real( kind = core_rknd ), dimension(gr%nz, hydromet_dim), intent(in) :: &
      wphydrometp, & ! Covariance of w and a hydrometeor   [(m/s) <hm units>]
      wp2hmp,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
      rtphmp,      & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp        ! Covariance of thl and a hydrometeor [K <hm units>]

    real( kind = core_rknd ), intent(in) ::  &
      wpthlp_sfc,   & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc,    & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc,     & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc        ! v'w' at surface          [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) :: &
      sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), intent(in), dimension(sclr_dim) ::  &
      wpsclrp_sfc      ! Scalar flux at surface         [{units vary} m/s]

    real( kind = core_rknd ), intent(in) :: &
      upwp_sfc_pert, & ! pertubed u'w' at surface    [m^2/s^2]
      vpwp_sfc_pert    ! pertubed v'w' at surface    [m^2/s^2]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), intent(in), dimension(edsclr_dim) ::  &
      wpedsclrp_sfc    ! Eddy-Scalar flux at surface    [{units vary} m/s]

    ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref, & ! Initial liquid water potential temperature   [K]
      um_ref,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref,   & ! Initial v wind; Michael Falk                 [m/s]
      ug,       & ! u geostrophic wind                           [m/s]
      vg          ! v geostrophic wind                           [m/s]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), intent(in) :: &
      host_dx,  & ! East-West horizontal grid spacing     [m]
      host_dy     ! North-South horizontal grid spacing   [m]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    real( kind = core_rknd ), intent(in) :: &
      lmin    ! Min. value for the length scale    [m]

    type( clubb_config_flags_type ), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags


    !!! Input/Output Variables
    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  &
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      up3,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2,     & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrp3,    & ! sclr'^3 (thermodynamic levels)       [{units vary}^3]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  &
      p_in_Pa, & ! Air pressure (thermodynamic levels)       [Pa]
      exner      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  &
      rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      wp2thvp,    & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) :: &
      sclrpthvp     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  &
      wp2rtp,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      uprcp,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef,           & ! Coef of X'r_c' in Eq. (34) (t-levs.) [K/(kg/kg)]
      wp4,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wpup2,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      wp2up2,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2,            & ! w'^2 v'^2 (momentum levels)          [m^4/s^4]
      ice_supersat_frac    ! ice cloud fraction (thermo. levels)  [-]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) :: &
      um_pert,   & ! perturbed <u>       [m/s]
      vm_pert,   & ! perturbed <v>       [m/s]
      upwp_pert, & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert    ! perturbed <v'w'>    [m^2/s^2]

    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! PDF parameters (thermodynamic levels)    [units vary]
      pdf_params_zm    ! PDF parameters on momentum levels        [units vary]

    type(implicit_coefs_terms), intent(inout) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

#ifdef GFDL
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

      real( kind = core_rknd ), intent(inout), dimension(gr%nz,edsclr_dim) :: &
      edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  &
      rcm_in_layer, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  &
      wprcp,             & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      w_up_in_cloud,     & ! Average upward velocity within liquid cloud   [m/s]
      invrs_tau_zm         ! One divided by tau on zm levels               [1/s]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      Kh_zt, & ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]
      Kh_zm    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

#ifdef CLUBB_CAM
    real( kind = core_rknd), intent(out), dimension(gr%nz) :: &
      qclvar        ! cloud water variance
#endif

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      thlprcp    ! thl'rc'              [K kg/kg]

    !!! Output Variable 
    integer, intent(inout) :: err_code_api ! Diagnostic, for if some calculation goes amiss.

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), intent(inOUT), dimension(gr%nz, min(1,sclr_dim) , 2) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
    logical, intent(in)                 ::  do_liquid_only_in_clubb
#endif


    ! -------------- Local Variables --------------
    type(stats), dimension(1) :: &
      stats_zt_col, &
      stats_zm_col, &
      stats_sfc_col

    type(grid), target, dimension(1) :: &
      gr_col

    real( kind = core_rknd ), dimension(1) ::  &
      fcor_col,  &          ! Coriolis forcing             [s^-1]
      sfc_elevation_col     ! Elevation of ground level    [m AMSL]

    ! Input Variables
    real( kind = core_rknd ), dimension(1,gr%nz) ::  &
      thlm_forcing_col,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing_col,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing_col,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing_col,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wprtp_forcing_col,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing_col,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing_col,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing_col,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing_col, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm_col,           & ! w mean wind component on momentum levels  [m/s]
      wm_zt_col,           & ! w mean wind component on thermo. levels   [m/s]
      rho_zm_col,          & ! Air density on momentum levels            [kg/m^3]
      rho_col,             & ! Air density on thermodynamic levels       [kg/m^3]
      rho_ds_zm_col,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt_col,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm_col, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt_col, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zm_col,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt_col,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      rfrzm_col              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), dimension(1,gr%nz,hydromet_dim) :: &
      hydromet_col           ! Collection of hydrometeors                [units vary]

    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      radf_col          ! Buoyancy production at the CL top due to LW radiative cooling [m^2/s^3]

#ifdef CLUBBND_CAM 
    real( kind = core_rknd ), dimension(1) :: & 
      varmu_col 
#endif 

    real( kind = core_rknd ), dimension(1,gr%nz, hydromet_dim) :: &
      wphydrometp_col, & ! Covariance of w and a hydrometeor   [(m/s) <hm units>]
      wp2hmp_col,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
      rtphmp_col,      & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp_col        ! Covariance of thl and a hydrometeor [K <hm units>]

    real( kind = core_rknd ), dimension(1) ::  &
      wpthlp_sfc_col,   & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc_col,    & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc_col,     & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc_col        ! v'w' at surface          [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), dimension(1,gr%nz,sclr_dim) :: &
      sclrm_forcing_col    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), dimension(1,sclr_dim) ::  &
      wpsclrp_sfc_col      ! Scalar flux at surface         [{units vary} m/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), dimension(1,gr%nz,edsclr_dim) :: &
      edsclrm_forcing_col  ! Eddy passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), dimension(1) :: &
      upwp_sfc_pert_col, & ! pertubed u'w' at surface    [m^2/s^2]
      vpwp_sfc_pert_col    ! pertubed v'w' at surface    [m^2/s^2]

    real( kind = core_rknd ), dimension(1,edsclr_dim) ::  &
      wpedsclrp_sfc_col    ! Eddy-Scalar flux at surface    [{units vary} m/s]

    ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
    real( kind = core_rknd ), dimension(1,gr%nz) ::  &
      rtm_ref_col,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref_col, & ! Initial liquid water potential temperature   [K]
      um_ref_col,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref_col,   & ! Initial v wind; Michael Falk                 [m/s]
      ug_col,       & ! u geostrophic wind                           [m/s]
      vg_col          ! v geostrophic wind                           [m/s]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), dimension(1) :: &
      host_dx_col,  & ! East-West horizontal grid spacing     [m]
      host_dy_col     ! North-South horizontal grid spacing   [m]

    type(nu_vertical_res_dep), dimension(1) :: &
      nu_vert_res_dep_col    ! Vertical resolution dependent nu values

    real( kind = core_rknd ), dimension(1) :: &
      lmin_col    ! Min. value for the length scale    [m]


    !!! Input/Output Variables
    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), dimension(1,gr%nz) ::  &
      um_col,      & ! u mean wind component (thermodynamic levels)   [m/s]
      upwp_col,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vm_col,      & ! v mean wind component (thermodynamic levels)   [m/s]
      vpwp_col,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2_col,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2_col,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      up3_col,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3_col,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm_col,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp_col,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      thlm_col,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp_col,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2_col,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      rtp3_col,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp2_col,   & ! th_l'^2 (momentum levels)                      [K^2]
      thlp3_col,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      rtpthlp_col, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2_col,     & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3_col        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    ! Passive scalar variables
    real( kind = core_rknd ), dimension(1,gr%nz,sclr_dim) :: &
      sclrm_col,     & ! Passive scalar mean (thermo. levels) [units vary]
      wpsclrp_col,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2_col,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrp3_col,    & ! sclr'^3 (thermodynamic levels)       [{units vary}^3]
      sclrprtp_col,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp_col    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), dimension(1,gr%nz) ::  &
      p_in_Pa_col, & ! Air pressure (thermodynamic levels)       [Pa]
      exner_col      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), dimension(1,gr%nz) ::  &
      rcm_col,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac_col, & ! cloud fraction (thermodynamic levels)          [-]
      wpthvp_col,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      wp2thvp_col,    & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      rtpthvp_col,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp_col      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), dimension(1,gr%nz,sclr_dim) :: &
      sclrpthvp_col     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), dimension(1,gr%nz) ::  &
      wp2rtp_col,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp_col,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      uprcp_col,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp_col,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef_col,           & ! Coef of X'r_c' in Eq. (34) (t-levs.) [K/(kg/kg)]
      wp4_col,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wpup2_col,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2_col,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      wp2up2_col,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2_col,            & ! w'^2 v'^2 (momentum levels)          [m^4/s^4]
      ice_supersat_frac_col    ! ice cloud fraction (thermo. levels)  [-]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      um_pert_col,   & ! perturbed <u>       [m/s]
      vm_pert_col,   & ! perturbed <v>       [m/s]
      upwp_pert_col, & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert_col    ! perturbed <v'w'>    [m^2/s^2]

    type(implicit_coefs_terms), dimension(1) :: &
      pdf_implicit_coefs_terms_col    ! Implicit coefs / explicit terms [units vary]

#ifdef GFDL
    real( kind = core_rknd ), dimension(1,gr%nz,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only_col  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

      real( kind = core_rknd ), dimension(1,gr%nz,edsclr_dim) :: &
      edsclrm_col   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), dimension(1,gr%nz) ::  &
      rcm_in_layer_col, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover_col     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), dimension(1,gr%nz) ::  &
      wprcp_col,             & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      w_up_in_cloud_col,     & ! Average upward velocity within liquid cloud   [m/s]
      invrs_tau_zm_col         ! One divided by tau on zm levels               [1/s]

    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      Kh_zt_col, & ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]
      Kh_zm_col    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

#ifdef CLUBB_CAM
    real( kind = core_rknd), dimension(1,gr%nz) :: &
      qclvar_col        ! cloud water variance
#endif

    real( kind = core_rknd ), dimension(1,gr%nz) :: &
      thlprcp_col    ! thl'rc'              [K kg/kg]
      
#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), dimension(1,gr%nz, min(1,sclr_dim) , 2) :: &
      RH_crit_col  ! critical relative humidity for droplet and ice nucleation
    logical, intent(in)                 ::  do_liquid_only_in_clubb
#endif

    integer :: i


    gr_col(1) = gr
    fcor_col(1) = fcor
    sfc_elevation_col(1) = sfc_elevation
    
    thlm_forcing_col(1,:) = thlm_forcing
    rtm_forcing_col(1,:) = rtm_forcing
    um_forcing_col(1,:) = um_forcing
    vm_forcing_col(1,:) = vm_forcing
    sclrm_forcing_col(1,:,:) = sclrm_forcing
    edsclrm_forcing_col(1,:,:) = edsclrm_forcing
    wprtp_forcing_col(1,:) = wprtp_forcing
    wpthlp_forcing_col(1,:) = wpthlp_forcing
    rtp2_forcing_col(1,:) = rtp2_forcing
    thlp2_forcing_col(1,:) = thlp2_forcing
    rtpthlp_forcing_col(1,:) = rtpthlp_forcing
    wm_zm_col(1,:) = wm_zm
    wm_zt_col(1,:) = wm_zt
    
    wpthlp_sfc_col(1) = wpthlp_sfc
    wprtp_sfc_col(1) = wprtp_sfc
    upwp_sfc_col(1) = upwp_sfc
    vpwp_sfc_col(1) = vpwp_sfc
    
    wpsclrp_sfc_col(1,:) = wpsclrp_sfc
    wpedsclrp_sfc_col(1,:) = wpedsclrp_sfc
    
    upwp_sfc_pert_col(1) = upwp_sfc_pert
    vpwp_sfc_pert_col(1) = vpwp_sfc_pert

    rtm_ref_col(1,:) = rtm_ref
    thlm_ref_col(1,:) = thlm_ref
    um_ref_col(1,:) = um_ref
    vm_ref_col(1,:) = vm_ref
    ug_col(1,:) = ug
    vg_col(1,:) = vg
    
    p_in_Pa_col(1,:) = p_in_Pa
    rho_zm_col(1,:) = rho_zm
    rho_col(1,:) = rho
    exner_col(1,:) = exner
    rho_ds_zm_col(1,:) = rho_ds_zm
    rho_ds_zt_col(1,:) = rho_ds_zt
    invrs_rho_ds_zm_col(1,:) = invrs_rho_ds_zm
    invrs_rho_ds_zt_col(1,:) = invrs_rho_ds_zt
    thv_ds_zm_col(1,:) = thv_ds_zm
    thv_ds_zt_col(1,:) = thv_ds_zt
    rfrzm_col(1,:) = rfrzm
    
    hydromet_col(1,:,:) = hydromet
    
    radf_col(1,:) = radf
#ifdef CLUBBND_CAM
    varmu_col(1) = varmu
#endif
    wphydrometp_col(1,:,:) = wphydrometp
    wp2hmp_col(1,:,:) = wp2hmp
    rtphmp_col(1,:,:) = rtphmp
    thlphmp_col(1,:,:) = thlphmp
    
    host_dx_col(1) = host_dx
    host_dy_col(1) = host_dy
    nu_vert_res_dep_col(1) = nu_vert_res_dep
    lmin_col(1) = lmin
    
    stats_zt_col(1) = stats_zt
    stats_zm_col(1) = stats_zm
    stats_sfc_col(1) = stats_sfc
    
    um_col(1,:) = um
    vm_col(1,:) = vm
    upwp_col(1,:) = upwp
    vpwp_col(1,:) = vpwp
    up2_col(1,:) = up2
    vp2_col(1,:) = vp2
    up3_col(1,:) = up3
    vp3_col(1,:) = vp3
    thlm_col(1,:) = thlm
    rtm_col(1,:) = rtm
    wprtp_col(1,:) = wprtp
    wpthlp_col(1,:) = wpthlp
    wp2_col(1,:) = wp2
    wp3_col(1,:) = wp3
    rtp2_col(1,:) = rtp2
    rtp3_col(1,:) = rtp3
    thlp2_col(1,:) = thlp2
    thlp3_col(1,:) = thlp3
    rtpthlp_col(1,:) = rtpthlp
    
    sclrm_col(1,:,:) = sclrm
#ifdef GFDL
    sclrm_trsport_only_col(1,:,:) = sclrm_trsport_only
#endif
    sclrp2_col(1,:,:) = sclrp2
    sclrp3_col(1,:,:) = sclrp3
    sclrprtp_col(1,:,:) = sclrprtp
    sclrpthlp_col(1,:,:) = sclrpthlp
    wpsclrp_col(1,:,:) = wpsclrp
    edsclrm_col(1,:,:) = edsclrm
    
    rcm_col(1,:) = rcm
    cloud_frac_col(1,:) = cloud_frac
    wpthvp_col(1,:) = wpthvp
    wp2thvp_col(1,:) = wp2thvp
    rtpthvp_col(1,:) = rtpthvp
    thlpthvp_col(1,:) = thlpthvp
    sclrpthvp_col(1,:,:) = sclrpthvp
    wp2rtp_col(1,:) = wp2rtp
    wp2thlp_col(1,:) = wp2thlp
    uprcp_col(1,:) = uprcp
    vprcp_col(1,:) = vprcp
    rc_coef_col(1,:) = rc_coef
    wp4_col(1,:) = wp4
    wpup2_col(1,:) = wpup2
    wpvp2_col(1,:) = wpvp2
    wp2up2_col(1,:) = wp2up2
    wp2vp2_col(1,:) = wp2vp2
    ice_supersat_frac_col(1,:) = ice_supersat_frac
    um_pert_col(1,:) = um_pert
    vm_pert_col(1,:) = vm_pert
    upwp_pert_col(1,:) = upwp_pert
    vpwp_pert_col(1,:) = vpwp_pert
    pdf_implicit_coefs_terms_col(1) = pdf_implicit_coefs_terms
#ifdef GFDL
    RH_crit_col(1,:,:,:) = RH_crit
#endif
    Kh_zm_col(1,:) = Kh_zm
    Kh_zt_col(1,:) = Kh_zt
#ifdef CLUBB_CAM
    qclvar_col(1,:) = qclvar
#endif
    thlprcp_col(1,:) = thlprcp
    wprcp_col(1,:) = wprcp
    w_up_in_cloud_col(1,:) = w_up_in_cloud
    rcm_in_layer_col(1,:) = rcm_in_layer
    cloud_cover_col(1,:) = cloud_cover
    invrs_tau_zm_col(1,:) = invrs_tau_zm

    call advance_clubb_core( gr_col, gr%nz, 1, &
      l_implemented, dt, fcor_col, sfc_elevation_col, hydromet_dim, & ! intent(in)
      thlm_forcing_col, rtm_forcing_col, um_forcing_col, vm_forcing_col, &    ! intent(in)
      sclrm_forcing_col, edsclrm_forcing_col, wprtp_forcing_col, &        ! intent(in)
      wpthlp_forcing_col, rtp2_forcing_col, thlp2_forcing_col, &          ! intent(in)
      rtpthlp_forcing_col, wm_zm_col, wm_zt_col, &                        ! intent(in)
      wpthlp_sfc_col, wprtp_sfc_col, upwp_sfc_col, vpwp_sfc_col, &            ! intent(in)
      wpsclrp_sfc_col, wpedsclrp_sfc_col, &                           ! intent(in)
      upwp_sfc_pert_col, vpwp_sfc_pert_col, &                         ! intent(in)
      rtm_ref_col, thlm_ref_col, um_ref_col, vm_ref_col, ug_col, vg_col, &            ! Intent(in)
      p_in_Pa_col, rho_zm_col, rho_col, exner_col, &                          ! intent(in)
      rho_ds_zm_col, rho_ds_zt_col, invrs_rho_ds_zm_col, &                ! intent(in)
      invrs_rho_ds_zt_col, thv_ds_zm_col, thv_ds_zt_col, hydromet_col, &      ! intent(in)
      rfrzm_col, radf_col, &                                          ! intent(in)
#ifdef CLUBBND_CAM
      varmu_col, &
#endif
      wphydrometp_col, wp2hmp_col, rtphmp_col, thlphmp_col, &                 ! intent(in)
      host_dx_col, host_dy_col, &                                     ! intent(in)
      clubb_params, nu_vert_res_dep_col, lmin_col, &                  ! intent(in)
      clubb_config_flags, &                                   ! intent(in)
      stats_zt_col, stats_zm_col, stats_sfc_col, &                        ! intent(inout)
      um_col, vm_col, upwp_col, vpwp_col, up2_col, vp2_col, up3_col, vp3_col, &               ! intent(inout)
      thlm_col, rtm_col, wprtp_col, wpthlp_col, &                             ! intent(inout)
      wp2_col, wp3_col, rtp2_col, rtp3_col, thlp2_col, thlp3_col, rtpthlp_col, &          ! intent(inout)
      sclrm_col,   &
#ifdef GFDL
      sclrm_trsport_only_col,  &  ! h1g, 2010-06-16      ! intent(inout)
#endif
      sclrp2_col, sclrp3_col, sclrprtp_col, sclrpthlp_col, &                  ! intent(inout)
      wpsclrp_col, edsclrm_col, &                                     ! intent(inout)
      rcm_col, cloud_frac_col, &                                      ! intent(inout)
      wpthvp_col, wp2thvp_col, rtpthvp_col, thlpthvp_col, &                   ! intent(inout)
      sclrpthvp_col, &                                            ! intent(inout)
      wp2rtp_col, wp2thlp_col, uprcp_col, vprcp_col, rc_coef_col, wp4_col, & ! intent(inout)
      wpup2_col, wpvp2_col, wp2up2_col, wp2vp2_col, ice_supersat_frac_col, & ! intent(inout)
      um_pert_col, vm_pert_col, upwp_pert_col, vpwp_pert_col, &            ! intent(inout)
      pdf_params, pdf_params_zm, &                            ! intent(inout)
      pdf_implicit_coefs_terms_col, &                             ! intent(inout)
#ifdef GFDL
               RH_crit_col, & !h1g, 2010-06-16                    ! intent(inout)
               do_liquid_only_in_clubb, &                     ! intent(in)
#endif
      Kh_zm_col, Kh_zt_col, &                                         ! intent(out)
#ifdef CLUBB_CAM
               qclvar_col, &                                      ! intent(out)
#endif
      thlprcp_col, wprcp_col, w_up_in_cloud_col, &                ! intent(out)
      rcm_in_layer_col, cloud_cover_col, invrs_tau_zm_col, &      ! intent(out)
      err_code_api )                                          ! intent(out)
    
    
    ! The following does not work for stats 
    !     stats_zt = stats_zt_col(1)
    !     stats_zm = stats_zm_col(1) 
    !     stats_sfc = stats_sfc_col(1)
    ! because of some mysterious pointer issue. However, the only thing that 
    ! updates in stats is the field values, so we can copy only those instead.
    if ( l_stats ) then 
      stats_zm%accum_field_values = stats_zm_col(1)%accum_field_values
      stats_zm%accum_num_samples = stats_zm_col(1)%accum_num_samples
      
      stats_zt%accum_field_values = stats_zt_col(1)%accum_field_values
      stats_zt%accum_num_samples = stats_zt_col(1)%accum_num_samples
      
      stats_sfc%accum_field_values = stats_sfc_col(1)%accum_field_values
      stats_sfc%accum_num_samples = stats_sfc_col(1)%accum_num_samples
    end if
      
      
    um = um_col(1,:)
    upwp = upwp_col(1,:)
    vm = vm_col(1,:)
    vpwp = vpwp_col(1,:)
    up2 = up2_col(1,:)
    vp2 = vp2_col(1,:)
    up3 = up3_col(1,:)
    vp3 = vp3_col(1,:)
    rtm = rtm_col(1,:)
    wprtp = wprtp_col(1,:)
    thlm = thlm_col(1,:)
    wpthlp = wpthlp_col(1,:)
    rtp2 = rtp2_col(1,:)
    rtp3 = rtp3_col(1,:)
    thlp2 = thlp2_col(1,:)
    thlp3 = thlp3_col(1,:)
    rtpthlp = rtpthlp_col(1,:)
    wp2 = wp2_col(1,:)
    wp3 = wp3_col(1,:)
    sclrm = sclrm_col(1,:,:)
    wpsclrp = wpsclrp_col(1,:,:)
    sclrp2 = sclrp2_col(1,:,:)
    sclrp3 = sclrp3_col(1,:,:)
    sclrprtp = sclrprtp_col(1,:,:)
    sclrpthlp = sclrpthlp_col(1,:,:)
    p_in_Pa = p_in_Pa_col(1,:)
    exner = exner_col(1,:)
    rcm = rcm_col(1,:)
    cloud_frac = cloud_frac_col(1,:)
    wpthvp = wpthvp_col(1,:)
    wp2thvp = wp2thvp_col(1,:)
    rtpthvp = rtpthvp_col(1,:)
    thlpthvp = thlpthvp_col(1,:)
    sclrpthvp = sclrpthvp_col(1,:,:)
    wp2rtp = wp2rtp_col(1,:)
    wp2thlp = wp2thlp_col(1,:)
    uprcp = uprcp_col(1,:)
    vprcp = vprcp_col(1,:)
    rc_coef = rc_coef_col(1,:)
    wp4 = wp4_col(1,:)
    wpup2 = wpup2_col(1,:)
    wpvp2 = wpvp2_col(1,:)
    wp2vp2 = wp2vp2_col(1,:)
    wp2up2 = wp2up2_col(1,:)
    ice_supersat_frac = ice_supersat_frac_col(1,:)
    pdf_implicit_coefs_terms = pdf_implicit_coefs_terms_col(1)
#ifdef GFDL
    sclrm_trsport_only = sclrm_trsport_only_col(1,:,:)
#endif
    edsclrm = edsclrm_col(1,:,:)
    rcm_in_layer = rcm_in_layer_col(1,:)
    cloud_cover  = cloud_cover_col(1,:)
    wprcp = wprcp_col(1,:)
    w_up_in_cloud = w_up_in_cloud_col(1,:)
    invrs_tau_zm = invrs_tau_zm_col(1,:)
    Kh_zt = Kh_zt_col(1,:)
    Kh_zm = Kh_zm_col(1,:)
#ifdef CLUBB_CAM
    qclvar = qclvar_col(1,:)
#endif
    thlprcp = thlprcp_col(1,:)
#ifdef GFDL
    RH_crit = RH_crit_col(1,:,:,:)
#endif

  end subroutine advance_clubb_core_api_single_col
  
  subroutine advance_clubb_core_api_multi_col( gr, nz, ngrdcol, &
    l_implemented, dt, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! intent(in)
    sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! intent(in)
    wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! intent(in)
    rtpthlp_forcing, wm_zm, wm_zt, &                        ! intent(in)
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &            ! intent(in)
    wpsclrp_sfc, wpedsclrp_sfc, &                           ! intent(in)
    upwp_sfc_pert, vpwp_sfc_pert, &                         ! intent(in)
    rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &            ! Intent(in)
    p_in_Pa, rho_zm, rho, exner, &                          ! intent(in)
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! intent(in)
    invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &      ! intent(in)
    rfrzm, radf, &                                          ! intent(in)
#ifdef CLUBBND_CAM
    varmu, &                                                ! intent(in)
#endif
    wphydrometp, wp2hmp, rtphmp, thlphmp, &                 ! intent(in)
    host_dx, host_dy, &                                     ! intent(in)
    clubb_params, nu_vert_res_dep, lmin, &                  ! intent(in)
    clubb_config_flags, &                                   ! intent(in)
    stats_zt, stats_zm, stats_sfc, &                        ! intent(inout)
    um, vm, upwp, vpwp, up2, vp2, up3, vp3, &               ! intent(inout)
    thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! intent(inout)
    sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16    ! intent(inout)
#endif
    sclrp2, sclrp3, sclrprtp, sclrpthlp, &                  ! intent(inout)
    wpsclrp, edsclrm, err_code_api, &                       ! intent(inout)
    rcm, cloud_frac, &                                      ! intent(inout)
    wpthvp, wp2thvp, rtpthvp, thlpthvp, &                   ! intent(inout)
    sclrpthvp, &                                            ! intent(inout)
    wp2rtp, wp2thlp, uprcp, vprcp, rc_coef, wp4, &          ! intent(inout)
    wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &      ! intent(inout)
    um_pert, vm_pert, upwp_pert, vpwp_pert, &               ! intent(inout)
    pdf_params, pdf_params_zm, &                            ! intent(inout)
    pdf_implicit_coefs_terms, &                             ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                  ! intent(inout)
               do_liquid_only_in_clubb, &                   ! intent(in)
#endif
    Kh_zm, Kh_zt, &                                         ! intent(out)
#ifdef CLUBB_CAM
    qclvar, &                                               ! intent(out)
#endif
    thlprcp, wprcp, w_up_in_cloud, &                        ! intent(out)
    rcm_in_layer, cloud_cover, invrs_tau_zm )               ! intent(out)

    use advance_clubb_core_module, only : advance_clubb_core

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type(s)

    use parameters_model, only: &
        sclr_dim, & ! Variable(s)
        edsclr_dim

    use model_flags, only: &
        clubb_config_flags_type


    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol   ! Number of grid columns

    type(stats), target, intent(inout), dimension(ngrdcol) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    type(grid), target, intent(in), dimension(ngrdcol) :: gr
    
      !!! Input Variables
    logical, intent(in) ::  &
      l_implemented ! Is this part of a larger host model (T/F) ?

    real( kind = core_rknd ), intent(in) ::  &
      dt  ! Current timestep duration    [s]
      
    real( kind = core_rknd ), intent(in), dimension(ngrdcol) ::  &
      fcor, &           ! Coriolis forcing             [s^-1]
      sfc_elevation     ! Elevation of ground level    [m AMSL]

    integer, intent(in) :: &
      hydromet_dim      ! Total number of hydrometeors          [#]

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) ::  &
      thlm_forcing,    & ! theta_l forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! u wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! v wind forcing (thermodynamic levels)     [m/s/s]
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing,    & ! <r_t'^2> forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! <th_l'^2> forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm,           & ! w mean wind component on momentum levels  [m/s]
      wm_zt,           & ! w mean wind component on thermo. levels   [m/s]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ momentum levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs.  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo. levs.  [K]
      rfrzm              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      hydromet           ! Collection of hydrometeors                [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      radf          ! Buoyancy production at the CL top due to LW radiative cooling [m^2/s^3]

#ifdef CLUBBND_CAM 
    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: & 
      varmu 
#endif 

    real( kind = core_rknd ), dimension(ngrdcol,nz, hydromet_dim), intent(in) :: &
      wphydrometp, & ! Covariance of w and a hydrometeor   [(m/s) <hm units>]
      wp2hmp,      & ! Third moment: <w'^2> * <hydro.'>    [(m/s)^2 <hm units>]
      rtphmp,      & ! Covariance of rt and a hydrometeor  [(kg/kg) <hm units>]
      thlphmp        ! Covariance of thl and a hydrometeor [K <hm units>]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol) ::  &
      wpthlp_sfc,   & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc,    & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc,     & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc        ! v'w' at surface          [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz,sclr_dim) :: &
      sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,sclr_dim) ::  &
      wpsclrp_sfc      ! Scalar flux at surface         [{units vary} m/s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
      upwp_sfc_pert, & ! pertubed u'w' at surface    [m^2/s^2]
      vpwp_sfc_pert    ! pertubed v'w' at surface    [m^2/s^2]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,edsclr_dim) ::  &
      wpedsclrp_sfc    ! Eddy-Scalar flux at surface    [{units vary} m/s]

    ! Reference profiles (used for nudging, sponge damping, and Coriolis effect)
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) ::  &
      rtm_ref,  & ! Initial total water mixing ratio             [kg/kg]
      thlm_ref, & ! Initial liquid water potential temperature   [K]
      um_ref,   & ! Initial u wind; Michael Falk                 [m/s]
      vm_ref,   & ! Initial v wind; Michael Falk                 [m/s]
      ug,       & ! u geostrophic wind                           [m/s]
      vg          ! v geostrophic wind                           [m/s]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
      host_dx,  & ! East-West horizontal grid spacing     [m]
      host_dy     ! North-South horizontal grid spacing   [m]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(nu_vertical_res_dep), intent(in), dimension(ngrdcol) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    real( kind = core_rknd ), intent(in), dimension(ngrdcol) :: &
      lmin    ! Min. value for the length scale    [m]

    type( clubb_config_flags_type ), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags


    !!! Input/Output Variables
    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) ::  &
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      up3,     & ! u'^3 (thermodynamic levels)                    [m^3/s^3]
      vp3,     & ! v'^3 (thermodynamic levels)                    [m^3/s^3]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp,  & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      rtpthlp, & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2,     & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrp3,    & ! sclr'^3 (thermodynamic levels)       [{units vary}^3]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) ::  &
      p_in_Pa, & ! Air pressure (thermodynamic levels)       [Pa]
      exner      ! Exner function (thermodynamic levels)     [-]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) ::  &
      rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      wp2thvp,    & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz,sclr_dim) :: &
      sclrpthvp     ! < sclr' th_v' > (momentum levels)   [units vary]

    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) ::  &
      wp2rtp,            & ! w'^2 rt' (thermodynamic levels)      [m^2/s^2 kg/kg]
      wp2thlp,           & ! w'^2 thl' (thermodynamic levels)     [m^2/s^2 K]
      uprcp,             & ! < u' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      vprcp,             & ! < v' r_c' > (momentum levels)        [(m/s)(kg/kg)]
      rc_coef,           & ! Coef of X'r_c' in Eq. (34) (t-levs.) [K/(kg/kg)]
      wp4,               & ! w'^4 (momentum levels)               [m^4/s^4]
      wpup2,             & ! w'u'^2 (thermodynamic levels)        [m^3/s^3]
      wpvp2,             & ! w'v'^2 (thermodynamic levels)        [m^3/s^3]
      wp2up2,            & ! w'^2 u'^2 (momentum levels)          [m^4/s^4]
      wp2vp2,            & ! w'^2 v'^2 (momentum levels)          [m^4/s^4]
      ice_supersat_frac    ! ice cloud fraction (thermo. levels)  [-]

    ! Variables used to track perturbed version of winds.
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz) :: &
      um_pert,   & ! perturbed <u>       [m/s]
      vm_pert,   & ! perturbed <v>       [m/s]
      upwp_pert, & ! perturbed <u'w'>    [m^2/s^2]
      vpwp_pert    ! perturbed <v'w'>    [m^2/s^2]

    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! PDF parameters (thermodynamic levels)    [units vary]
      pdf_params_zm    ! PDF parameters on momentum levels        [units vary]

    type(implicit_coefs_terms), intent(inout), dimension(ngrdcol) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

#ifdef GFDL
    real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

      real( kind = core_rknd ), intent(inout), dimension(ngrdcol,nz,edsclr_dim) :: &
      edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nz) ::  &
      rcm_in_layer, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nz) ::  &
      wprcp,             & ! w'r_c' (momentum levels)              [(kg/kg) m/s]
      w_up_in_cloud,     & ! Average upward velocity within liquid cloud   [m/s]
      invrs_tau_zm         ! One divided by tau on zm levels               [1/s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      Kh_zt, & ! Eddy diffusivity coefficient on thermodynamic levels   [m^2/s]
      Kh_zm    ! Eddy diffusivity coefficient on momentum levels        [m^2/s]

#ifdef CLUBB_CAM
    real( kind = core_rknd), intent(out), dimension(ngrdcol,nz) :: &
      qclvar        ! cloud water variance
#endif

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      thlprcp    ! thl'rc'              [K kg/kg]

    !!! Output Variable 
    integer, intent(inout) :: err_code_api ! Diagnostic, for if some calculation goes amiss.

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), intent(inOUT), dimension(ngrdcol,nz, min(1,sclr_dim) , 2) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
    logical, intent(in)                 ::  do_liquid_only_in_clubb
#endif

    call advance_clubb_core( gr, nz, ngrdcol, &
      l_implemented, dt, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
      thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! intent(in)
      sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! intent(in)
      wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! intent(in)
      rtpthlp_forcing, wm_zm, wm_zt, &                        ! intent(in)
      wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &            ! intent(in)
      wpsclrp_sfc, wpedsclrp_sfc, &                           ! intent(in)
      upwp_sfc_pert, vpwp_sfc_pert, &                         ! intent(in)
      rtm_ref, thlm_ref, um_ref, vm_ref, ug, vg, &            ! Intent(in)
      p_in_Pa, rho_zm, rho, exner, &                          ! intent(in)
      rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! intent(in)
      invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &      ! intent(in)
      rfrzm, radf, &                                          ! intent(in)
#ifdef CLUBBND_CAM
      varmu, &
#endif
      wphydrometp, wp2hmp, rtphmp, thlphmp, &                 ! intent(in)
      host_dx, host_dy, &                                     ! intent(in)
      clubb_params, nu_vert_res_dep, lmin, &                  ! intent(in)
      clubb_config_flags, &                                   ! intent(in)
      stats_zt, stats_zm, stats_sfc, &                        ! intent(inout)
      um, vm, upwp, vpwp, up2, vp2, up3, vp3, &               ! intent(inout)
      thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
      wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! intent(inout)
      sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16      ! intent(inout)
#endif
      sclrp2, sclrp3, sclrprtp, sclrpthlp, &                  ! intent(inout)
      wpsclrp, edsclrm, &                                     ! intent(inout)
      rcm, cloud_frac, &                                      ! intent(inout)
      wpthvp, wp2thvp, rtpthvp, thlpthvp, &                   ! intent(inout)
      sclrpthvp, &                                            ! intent(inout)
      wp2rtp, wp2thlp, uprcp, vprcp, rc_coef, wp4, &          ! intent(inout)
      wpup2, wpvp2, wp2up2, wp2vp2, ice_supersat_frac, &      ! intent(inout)
      um_pert, vm_pert, upwp_pert, vpwp_pert, &               ! intent(inout)
      pdf_params, pdf_params_zm, &                            ! intent(inout)
      pdf_implicit_coefs_terms, &                             ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                    ! intent(inout)
               do_liquid_only_in_clubb, &                     ! intent(in)
#endif
      Kh_zm, Kh_zt, &                                         ! intent(out)
#ifdef CLUBB_CAM
               qclvar, &                                      ! intent(out)
#endif
      thlprcp, wprcp, w_up_in_cloud, &                        ! intent(out)
      rcm_in_layer, cloud_cover, invrs_tau_zm, &              ! intent(out)
      err_code_api )                                          ! intent(out)

  end subroutine advance_clubb_core_api_multi_col

  !================================================================================================
  ! setup_clubb_core - Sets up the model for execution.
  !================================================================================================

  subroutine setup_clubb_core_api( &
    nzmax, T0_in, ts_nudge_in,                          & ! intent(in)
    hydromet_dim_in, sclr_dim_in,                       & ! intent(in)
    sclr_tol_in, edsclr_dim_in, params,                 & ! intent(in)
    l_host_applies_sfc_fluxes,                          & ! intent(in)
    saturation_formula,                                 & ! intent(in)
    l_input_fields,                                     & ! intent(in)
#ifdef GFDL
    I_sat_sphum,                                        & ! intent(in)  h1g, 2010-06-16
#endif
    l_implemented, grid_type, deltaz, zm_init, zm_top,  & ! intent(in)
    momentum_heights, thermodynamic_heights,            & ! intent(in)
    sfc_elevation,                                      & ! intent(in)
    iiPDF_type,                                         & ! intent(in)
    ipdf_call_placement,                                & ! intent(in)
    l_predict_upwp_vpwp,                                & ! intent(in)
    l_min_xp2_from_corr_wx,                             & ! intent(in)
    l_prescribed_avg_deltaz,                            & ! intent(in)
    l_damp_wp2_using_em,                                & ! intent(in)
    l_stability_correct_tau_zm,                         & ! intent(in)
    l_enable_relaxed_clipping,                          & ! intent(in)
    l_diag_Lscale_from_tau,                             & ! intent(in)
#ifdef GFDL
    cloud_frac_min ,                                    & ! intent(in)  h1g, 2010-06-16
#endif
    gr, lmin, nu_vert_res_dep, err_code_api )             ! intent(out) 

    use advance_clubb_core_module, only : setup_clubb_core

    use parameter_indices, only:  &
        nparams ! Variable(s)
      
    use model_flags, only: &
        clubb_config_flags_type  ! Type

! TODO: This should be called from the api, but all the host models appear to call
!       it directly or not at all.
!   use model_flags, only: &
!     setup_model_flags    ! Subroutine

      implicit none

      type(grid), target, intent(inout) :: gr

    ! Input Variables

    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]

    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation  ! Elevation of ground level    [m AMSL]

    logical, intent(in) :: l_implemented   ! (T/F) CLUBB implemented in host model?

    ! If CLUBB is running on it's own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB model is running by itself, and is using an
    ! evenly-spaced grid (grid_type = 1), it needs the vertical
    ! grid spacing, momentum-level starting altitude, and maximum
    ! altitude as input.
    real( kind = core_rknd ), intent(in) :: &
      deltaz,   & ! Change in altitude per level           [m]
      zm_init,  & ! Initial grid altitude (momentum level) [m]
      zm_top      ! Maximum grid altitude (momentum level) [m]

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(nzmax) :: &
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    ! Model parameters
    real( kind = core_rknd ), intent(in) ::  &
      T0_in, ts_nudge_in

    integer, intent(in) :: &
      hydromet_dim_in,  & ! Number of hydrometeor species
      sclr_dim_in,      & ! Number of passive scalars
      edsclr_dim_in       ! Number of eddy-diff. passive scalars

    real( kind = core_rknd ), intent(in), dimension(sclr_dim_in) :: &
      sclr_tol_in    ! Thresholds for passive scalars

    real( kind = core_rknd ), intent(in), dimension(nparams) :: &
      params  ! Including C1, nu1, nu2, etc.

    ! Flags
    logical, intent(in) ::  &
      l_host_applies_sfc_fluxes ! Whether to apply for the surface flux

    character(len=*), intent(in) :: &
      saturation_formula ! Approximation for saturation vapor pressure

    logical, intent(in) ::  &
      l_input_fields    ! Flag for whether LES input fields are used

    integer, intent(in) :: &
      iiPDF_type,          & ! Selected option for the two-component normal
                             ! (double Gaussian) PDF type to use for the w,
                             ! rt, and theta-l (or w, chi, and eta) portion of
                             ! CLUBB's multivariate, two-component PDF.
      ipdf_call_placement    ! Selected option for the placement of the call to
                             ! CLUBB's PDF.

    logical, intent(in) :: &
      l_predict_upwp_vpwp,         & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                     ! alongside the advancement of <rt>, <w'rt'>, <thl>, <wpthlp>,
                                     ! <sclr>, and <w'sclr'> in subroutine advance_xm_wpxp.
                                     ! Otherwise, <u'w'> and <v'w'> are still approximated by eddy
                                     ! diffusivity when <u> and <v> are advanced in subroutine
                                     ! advance_windm_edsclrm.
      l_min_xp2_from_corr_wx,     & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                    ! thlp2) on keeping the overall correlation of w and x within
                                    ! the limits of -max_mag_correlation_flux to
                                    ! max_mag_correlation_flux.
      l_prescribed_avg_deltaz,    & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
      l_damp_wp2_using_em,        &
      l_stability_correct_tau_zm, &
      l_enable_relaxed_clipping,  & ! Flag to relax clipping on wpxp in
                                    ! xm_wpxp_clipping_and_stats
      l_diag_Lscale_from_tau        ! First diagnose dissipation time tau, and
                                    ! then diagnose the mixing length scale as
                                    ! Lscale = tau * tke

#ifdef GFDL
      logical, intent(in) :: &  ! h1g, 2010-06-16 begin mod
         I_sat_sphum

      real( kind = core_rknd ), intent(in) :: &
         cloud_frac_min         ! h1g, 2010-06-16 end mod
#endif

    ! Output variables 
    real( kind = core_rknd ), intent(out) :: &
      lmin    ! Min. value for the length scale    [m]

    type(nu_vertical_res_dep), intent(out) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(out) :: & 
      err_code_api   ! Diagnostic for a problem with the setup 

    call setup_clubb_core &
      ( nzmax, T0_in, ts_nudge_in,                          & ! intent(in)
      hydromet_dim_in, sclr_dim_in,                         & ! intent(in)
      sclr_tol_in, edsclr_dim_in, params,                   & ! intent(in)
      l_host_applies_sfc_fluxes,                            & ! intent(in)
      saturation_formula,                                   & ! intent(in)
      l_input_fields,                                       & ! intent(in)
#ifdef GFDL
      I_sat_sphum,                                          & ! intent(in)  h1g, 2010-06-16
#endif
      l_implemented, grid_type, deltaz, zm_init, zm_top,    & ! intent(in)
      momentum_heights, thermodynamic_heights,              & ! intent(in)
      sfc_elevation,                                        & ! intent(in)
      iiPDF_type,                                           & ! intent(in)
      ipdf_call_placement,                                  & ! intent(in)
      l_predict_upwp_vpwp,                                  & ! intent(in)
      l_min_xp2_from_corr_wx,                               & ! intent(in)
      l_prescribed_avg_deltaz,                              & ! intent(in)
      l_damp_wp2_using_em,                                  & ! intent(in)
      l_stability_correct_tau_zm,                           & ! intent(in)
      l_enable_relaxed_clipping,                            & ! intent(in)
      l_diag_Lscale_from_tau,                               & ! intent(in)
#ifdef GFDL
      cloud_frac_min,                                       & ! intent(in)  h1g, 2010-06-16
#endif
      gr, lmin, nu_vert_res_dep, err_code_api )               ! intent(out)

  end subroutine setup_clubb_core_api

  !================================================================================================
  ! cleanup_clubb_core_api - Frees memory used by the model.
  !================================================================================================

  subroutine cleanup_clubb_core_api( gr )

    use advance_clubb_core_module, only : cleanup_clubb_core

    

    implicit none

    type(grid), target, intent(inout) :: gr
    call cleanup_clubb_core( gr ) ! intent(inout)

  end subroutine cleanup_clubb_core_api

  !================================================================================================
  ! gregorian2julian_day - Computes the number of days since 1 January 4713 BC.
  !================================================================================================

  integer function gregorian2julian_day_api( &
    day, month, year )

    use calendar, only : gregorian2julian_day

    implicit none

    ! Input Variables
    integer, intent(in) ::  &
      day,        & ! Gregorian Calendar Day for given Month        [dd]
      month,      & ! Gregorian Calendar Month for given Year       [mm]
      year          ! Gregorian Calendar Year                       [yyyy]

    gregorian2julian_day_api = gregorian2julian_day( &
      day, month, year )
  end function gregorian2julian_day_api

  !================================================================================================
  ! compute_current_date - Computes the current date and the seconds since that date.
  !================================================================================================

  subroutine compute_current_date_api( &
    previous_day, previous_month, &
    previous_year,  &
    seconds_since_previous_date, &
    current_day, current_month, &
    current_year, &
    seconds_since_current_date )

    use calendar, only : compute_current_date

    implicit none

    ! Previous date
    integer, intent(in) :: &
      previous_day,    & ! Day of the month      [dd]
      previous_month,  & ! Month of the year     [mm]
      previous_year      ! Year                  [yyyy]

    real(kind=time_precision), intent(in) :: &
      seconds_since_previous_date ! [s]

    ! Output Variable(s)

    ! Current date
    integer, intent(out) :: &
      current_day,     & ! Day of the month      [dd]
      current_month,   & ! Month of the year     [mm]
      current_year       ! Year                  [yyyy]

    real(kind=time_precision), intent(out) :: &
      seconds_since_current_date

    call compute_current_date( &
      previous_day, previous_month, & ! intent(in)
      previous_year,  & ! intent(in)
      seconds_since_previous_date, & ! intent(in)
      current_day, current_month, & ! intent(out)
      current_year, & ! intent(out)
      seconds_since_current_date ) ! intent(out)
  end subroutine compute_current_date_api

  !================================================================================================
  ! leap_year - Determines if the given year is a leap year.
  !================================================================================================

  logical function leap_year_api( &
    year )

    use calendar, only : leap_year

    implicit none

    ! External
    intrinsic :: mod

    ! Input Variable(s)
    integer, intent(in) :: year ! Gregorian Calendar Year [yyyy]

    leap_year_api = leap_year( &
      year )
  end function leap_year_api

  !================================================================================================
  ! setup_corr_varnce_array - Creates a correlation array with x'^2/xm^2 variables on the diagonal
  !================================================================================================

  subroutine setup_corr_varnce_array_api( &
    input_file_cloud, input_file_below, iunit, &
    l_fix_w_chi_eta_correlations )

    use corr_varnce_module, only : setup_corr_varnce_array

    implicit none

    ! External
    intrinsic :: max, epsilon, trim

    ! Input Variables
    integer, intent(in) :: &
      iunit ! The file unit

    character(len=*), intent(in) :: &
      input_file_cloud, &  ! Path to the in cloud correlation file
      input_file_below     ! Path to the out of cloud correlation file

    logical, intent(in) :: &
      l_fix_w_chi_eta_correlations ! Use a fixed correlation for s and t Mellor(chi/eta)

    call setup_corr_varnce_array( &
      input_file_cloud, input_file_below, iunit, & ! intent(in)
      l_fix_w_chi_eta_correlations ) ! intent(in)

  end subroutine setup_corr_varnce_array_api

  !================================================================================================
  ! set_clubb_debug_level - Controls the importance of error messages sent to the console.
  !================================================================================================

  subroutine set_clubb_debug_level_api( &
    level )

    use error_code, only: &
        set_clubb_debug_level ! Procedure

    implicit none

    ! Input variable
    integer, intent(in) :: level ! The debug level being checked against the current setting

    call set_clubb_debug_level( &
      level ) ! intent(in)
  end subroutine set_clubb_debug_level_api

  !================================================================================================
  ! clubb_at_least_debug_level - Checks to see if clubb has been set to a specified debug level.
  !================================================================================================

  logical function clubb_at_least_debug_level_api( &
    level )
    
    use error_code, only: &
        clubb_at_least_debug_level ! Procedure

    implicit none

    ! Input variable
    integer, intent(in) :: level   ! The debug level being checked against the current setting

    clubb_at_least_debug_level_api = clubb_at_least_debug_level( &
      level )
  end function clubb_at_least_debug_level_api

  !================================================================================================
  ! fill_holes_driver - Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
  !================================================================================================

  subroutine fill_holes_driver_api( gr, &
    nz, dt, hydromet_dim,        & ! Intent(in)
    l_fill_holes_hm,             & ! Intent(in)
    rho_ds_zm, rho_ds_zt, exner, & ! Intent(in)
    stats_zt, &                    ! intent(inout)
    thlm_mc, rvm_mc, hydromet )    ! Intent(inout)

    use fill_holes, only : fill_holes_driver

    

    implicit none

    type(stats), target, intent(inout) :: &
      stats_zt

    type(grid), target, intent(in) :: gr

    intrinsic :: trim

    ! Input Variables
    integer, intent(in) :: hydromet_dim, nz

    logical, intent(in) :: l_fill_holes_hm

    real( kind = core_rknd ), intent(in) ::  &
      dt           ! Timestep         [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rho_ds_zm, & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt    ! Dry, static density on thermo. levels    [kg/m^3]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      exner  ! Exner function                                       [-]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz, hydromet_dim), intent(inout) :: &
      hydromet

    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      rvm_mc,  & ! Microphysics contributions to vapor water            [kg/kg/s]
      thlm_mc    ! Microphysics contributions to liquid potential temp. [K/s]

    call fill_holes_driver( gr,    & ! intent(in)
      nz, dt, hydromet_dim,        & ! Intent(in)
      l_fill_holes_hm,             & ! Intent(in)
      rho_ds_zm, rho_ds_zt, exner, & ! Intent(in)
      stats_zt,                    & ! intent(inout)
      thlm_mc, rvm_mc, hydromet )    ! Intent(inout)
  end subroutine fill_holes_driver_api

  !================================================================================================
  ! fill_holes_vertical - clips values of 'field' that are below 'threshold' as much as possible.
  !================================================================================================

  subroutine fill_holes_vertical_api( gr, &
    num_pts, threshold, field_grid, &
    rho_ds, rho_ds_zm, &
    field )

    use fill_holes, only : fill_holes_vertical

     ! Type

    implicit none

    type(grid), target, intent(in) :: gr

    ! Input variables
    integer, intent(in) :: &
      num_pts  ! The number of points on either side of the hole;
               ! Mass is drawn from these points to fill the hole.  []

    real( kind = core_rknd ), intent(in) :: &
      threshold  ! A threshold (e.g. w_tol*w_tol) below which field must not
                 ! fall                           [Units vary; same as field]

    character(len=2), intent(in) :: &
      field_grid ! The grid of the field, either stats_zt or stats_zm

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      rho_ds,    & ! Dry, static density on thermodynamic levels    [kg/m^3]
      rho_ds_zm    ! Dry, static density on momentum levels         [kg/m^3]

    ! Input/Output variable
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      field  ! The field (e.g. wp2) that contains holes [Units same as threshold]

    call fill_holes_vertical( gr, & ! intent(in)
      num_pts, threshold, field_grid, & ! intent(in)
      rho_ds, rho_ds_zm, & ! intent(in)
      field ) ! intent(inout)
  end subroutine fill_holes_vertical_api

  !=============================================================================
  ! fill_holes_hydromet - fills holes in a hydrometeor using mass from another
  ! hydrometeor that has the same phase.
  !=============================================================================
  subroutine fill_holes_hydromet_api( nz, hydromet_dim, hydromet, & ! Intent(in)
                                      hydromet_filled ) ! Intent(out)

    use fill_holes, only: &
        fill_holes_hydromet ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim, & ! Number of hydrometeor fields
      nz              ! Number of vertical grid levels

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet    ! Mean of hydrometeor fields    [units vary] 

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(out) :: &
      hydromet_filled ! Mean of hydrometeor fields after hole filling [un. vary]


    call fill_holes_hydromet( nz, hydromet_dim, hydromet, & ! Intent(in)
                              hydromet_filled ) ! Intent(out)


  end subroutine fill_holes_hydromet_api

  !================================================================================================
  ! vertical_integral - Computes the vertical integral.
  !================================================================================================

  function vertical_integral_api( &
    total_idx, rho_ds, &
    field, dz )

    use fill_holes, only : vertical_integral

    implicit none

    ! Input variables
    integer, intent(in) :: &
      total_idx  ! The total numer of indices within the range of averaging

    real( kind = core_rknd ), dimension(total_idx), intent(in) ::  &
      rho_ds,  & ! Dry, static density                   [kg/m^3]
      field,   & ! The field to be vertically averaged   [Units vary]
      dz         ! Level thickness                       [1/m]
    ! Note:  The rho_ds and field points need to be arranged from
    !        lowest to highest in altitude, with rho_ds(1) and
    !        field(1) actually their respective values at level k = begin_idx.

    real( kind = core_rknd ) :: &
      vertical_integral_api ! Integral in the numerator (see description)

    vertical_integral_api = vertical_integral( &
      total_idx, rho_ds, &
      field, dz )
  end function vertical_integral_api

  !================================================================================================
  ! setup_grid_heights - Sets the heights and interpolation weights of the column.
  !================================================================================================

  subroutine setup_grid_heights_api( &
    l_implemented, grid_type,  &
    deltaz, zm_init, momentum_heights,  &
    gr, thermodynamic_heights )

    use grid_class, only: & 
        grid, & ! Type
        setup_grid_heights
    
    use error_code, only : &
        clubb_fatal_error       ! Constant

    implicit none
   
    type(grid), target, intent(inout) :: gr

    ! Input Variables

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented

    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB model is running by itself, and is using an evenly-spaced
    ! grid (grid_type = 1), it needs the vertical grid spacing and
    ! momentum-level starting altitude as input.
    real( kind = core_rknd ), intent(in) ::  &
      deltaz,   & ! Vertical grid spacing                  [m]
      zm_init     ! Initial grid altitude (momentum level) [m]


    ! If the CLUBB parameterization is implemented in a host model, it needs to
    ! use the host model's momentum level altitudes and thermodynamic level
    ! altitudes.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on thermodynamic levels (grid_type = 2), it needs to use the
    ! thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on momentum levels (grid_type = 3), it needs to use the momentum
    ! level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  &
      momentum_heights,   & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

    call setup_grid_heights( &
      l_implemented, grid_type,  & ! intent(in)
      deltaz, zm_init, momentum_heights,  & ! intent(in)
      thermodynamic_heights, & ! intent(in)
      gr ) ! intent(inout)

    if ( err_code == clubb_fatal_error ) error stop

  end subroutine setup_grid_heights_api
  
  !================================================================================================
  ! setup_grid - This subroutine sets up the CLUBB vertical grid.
  !================================================================================================
  
  subroutine setup_grid_api( nzmax, sfc_elevation, l_implemented, &
                             grid_type, deltaz, zm_init, zm_top, &
                             momentum_heights, thermodynamic_heights, &
                             gr, begin_height, end_height )
                            
    use grid_class, only: & 
        grid, & ! Type
        setup_grid

    implicit none

    type(grid), target, intent(inout) :: gr

    ! Input Variables
    integer, intent(in) ::  & 
      nzmax  ! Number of vertical levels in grid      [#]

    real( kind = core_rknd ), intent(in) ::  &
      sfc_elevation  ! Elevation of ground level    [m AMSL]
      
    logical, intent(in) :: l_implemented
    
    integer, intent(in) :: grid_type
    
    real( kind = core_rknd ), intent(in) ::  & 
      deltaz,   & ! Vertical grid spacing                  [m]
      zm_init,  & ! Initial grid altitude (momentum level) [m]
      zm_top      ! Maximum grid altitude (momentum level) [m]
      
    real( kind = core_rknd ), intent(in), dimension(nzmax) ::  & 
      momentum_heights,   & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights ! Thermodynamic level altitudes (input) [m]

    integer, intent(out) :: &
      begin_height, &  ! Lower bound for *_heights arrays [-]
      end_height       ! Upper bound for *_heights arrays [-]


    call setup_grid( nzmax, sfc_elevation, l_implemented,     & ! intent(in)
                     grid_type, deltaz, zm_init, zm_top,      & ! intent(in)
                     momentum_heights, thermodynamic_heights, & ! intent(in)
                     gr, begin_height, end_height             ) ! intent(out)

  end subroutine setup_grid_api

  !================================================================================================
  ! lin_interpolate_two_points - Computes a linear interpolation of the value of a variable.
  !================================================================================================

  function lin_interpolate_two_points_api( &
    height_int, height_high, height_low, &
    var_high, var_low )

    use interpolation, only : lin_interpolate_two_points

    implicit none

    real( kind = core_rknd ), intent(in) :: &
      height_int,  & ! Height to be interpolated to     [m]
      height_high, & ! Height above the interpolation   [m]
      height_low,  & ! Height below the interpolation   [m]
      var_high,    & ! Variable above the interpolation [units vary]
      var_low        ! Variable below the interpolation [units vary]

    ! Output Variables
    real( kind = core_rknd ) :: lin_interpolate_two_points_api

    lin_interpolate_two_points_api = lin_interpolate_two_points( &
      height_int, height_high, height_low, &
      var_high, var_low )

  end function lin_interpolate_two_points_api

  !================================================================================================
  ! lin_interpolate_on_grid - Linear interpolation for 25 June 1996 altocumulus case.
  !================================================================================================

  subroutine lin_interpolate_on_grid_api( &
    nparam, xlist, tlist, xvalue, tvalue )

    use interpolation, only : lin_interpolate_on_grid

    implicit none

    ! Input Variables
    integer, intent(in) :: nparam ! Number of parameters in xlist and tlist

    ! Input/Output Variables
    real( kind = core_rknd ), intent(inout), dimension(nparam) ::  &
      xlist,  & ! List of x-values (independent variable)
      tlist     ! List of t-values (dependent variable)

    real( kind = core_rknd ), intent(in) ::  &
      xvalue  ! x-value at which to interpolate

    real( kind = core_rknd ), intent(inout) ::  &
      tvalue  ! t-value solved by interpolation

    call lin_interpolate_on_grid( &
      nparam, xlist, tlist, xvalue, & ! intent(in)
      tvalue ) ! intent(inout)

  end subroutine lin_interpolate_on_grid_api

  !================================================================================================
  ! read_parameters - Read a namelist containing the model parameters.
  !================================================================================================

  subroutine read_parameters_api( iunit, filename, &
                                  C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
                                  C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
                                  C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
                                  C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
                                  C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
                                  C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
                                  c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
                                  c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
                                  slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
                                  coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
                                  gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
                                  omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
                                  lambda0_stability_coef, mult_coef, taumin, taumax, &
                                  Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
                                  Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
                                  thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
                                  Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
                                  altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
                                  C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, & 
                                  C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
                                  C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
                                  C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
                                  Cx_min, Cx_max, Richardson_num_min, &
                                  Richardson_num_max, a3_coef_min, &
                                  params )

    use parameters_tunable, only : read_parameters

    use parameter_indices, only:  &
        nparams ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Input/Output variables
    real( kind = core_rknd ), intent(inout) :: & 
      C1, C1b, C1c, C2rt, C2thl, C2rtthl, & 
      C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
      C7, C7b, C7c, C8, C8b, C10, & 
      C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
      C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, & 
      C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
      c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8,  & 
      c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
      slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
      coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
      gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
      omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
      lambda0_stability_coef, mult_coef, taumin, taumax, Lscale_mu_coef, &
      Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, c_K10h, &
      thlp2_rad_coef, thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
      Skw_max_mag, xp3_coef_base, xp3_coef_slope, altitude_threshold, &
      rtp2_clip_coef, C_invrs_tau_bkgnd, C_invrs_tau_sfc, &
      C_invrs_tau_shear, C_invrs_tau_N2, C_invrs_tau_N2_wp2, &
      C_invrs_tau_N2_xp2, C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
      C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, a3_coef_min

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(nparams) :: params

    call read_parameters( iunit, filename, & ! intent(in)
                          C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
                          C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
                          C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
                          C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
                          C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
                          C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
                          c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
                          c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
                          slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
                          coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
                          gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
                          omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
                          lambda0_stability_coef, mult_coef, taumin, taumax, &
                          Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
                          Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
                          thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
                          Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
                          altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
                          C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, & 
                          C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
                          C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
                          C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
                          Cx_min, Cx_max, Richardson_num_min, &
                          Richardson_num_max, a3_coef_min, &
                          params ) ! intent(out)

  end subroutine read_parameters_api

  !================================================================================================
  ! setup_parameters - Sets up model parameters.
  !================================================================================================

  subroutine setup_parameters_api &
           ( deltaz, params, nzmax, &
             grid_type, momentum_heights, thermodynamic_heights, &
             l_prescribed_avg_deltaz, &
             lmin, nu_vert_res_dep, err_code_api )

    use parameters_tunable, only: &
        setup_parameters

    use parameter_indices, only:  &
        nparams ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      deltaz  ! Change per height level        [m]

    real( kind = core_rknd ), intent(in), dimension(nparams) :: &
      params  ! Tuneable model parameters      [-]

    ! Grid definition
    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]

    ! If CLUBB is running on its own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(nzmax) :: &
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    logical, intent(in) :: &
      l_prescribed_avg_deltaz ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz

    ! Output Variables 
    real( kind = core_rknd ), intent(out) :: &
      lmin    ! Min. value for the length scale    [m]

    type(nu_vertical_res_dep), intent(out) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(out) ::  & 	 	      
      err_code_api ! Error condition 

    call setup_parameters & 
            ( deltaz, params, nzmax, &
              grid_type, momentum_heights, thermodynamic_heights, &
              l_prescribed_avg_deltaz, &
              lmin, nu_vert_res_dep, err_code_api )

  end subroutine setup_parameters_api

  !================================================================================================
  ! adj_low_res_nu - Adjusts values of background eddy diffusivity based on vertical grid spacing.
  !================================================================================================

  subroutine adj_low_res_nu_api( nzmax, grid_type, deltaz,  & ! Intent(in)
                                 momentum_heights, thermodynamic_heights, & ! Intent(in)
                                 l_prescribed_avg_deltaz, mult_coef, &  ! Intent(in)
                                 nu1, nu2, nu6, nu8, nu9, nu10, nu_hm, &  ! Intent(in)
                                 nu_vert_res_dep )  ! Intent(out)

    use parameters_tunable, only : adj_low_res_nu

    implicit none

    ! Input Variables

    ! Grid definition
    integer, intent(in) :: nzmax  ! Vertical grid levels            [#]

    ! If CLUBB is running on it's own, this option determines
    ! if it is using:
    ! 1) an evenly-spaced grid,
    ! 2) a stretched (unevenly-spaced) grid entered on the
    !    thermodynamic grid levels (with momentum levels set
    !    halfway between thermodynamic levels), or
    ! 3) a stretched (unevenly-spaced) grid entered on the
    !    momentum grid levels (with thermodynamic levels set
    !    halfway between momentum levels).
    integer, intent(in) :: grid_type

    real( kind = core_rknd ), intent(in) ::  &
      deltaz  ! Change per height level        [m]

    ! If the CLUBB parameterization is implemented in a host model,
    ! it needs to use the host model's momentum level altitudes
    ! and thermodynamic level altitudes.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on thermodynamic levels (grid_type = 2),
    ! it needs to use the thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a
    ! stretched grid entered on momentum levels (grid_type = 3),
    ! it needs to use the momentum level altitudes as input.
    real( kind = core_rknd ), intent(in), dimension(nzmax) :: &
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    logical, intent(in) :: &
      l_prescribed_avg_deltaz ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz

    real( kind = core_rknd ), intent(in) :: &
      mult_coef, & ! CLUBB tunable parameter mult_coef
      nu1,       & ! CLUBB tunable parameter nu1
      nu2,       & ! CLUBB tunable parameter nu2
      nu6,       & ! CLUBB tunable parameter nu6
      nu8,       & ! CLUBB tunable parameter nu8
      nu9,       & ! CLUBB tunable parameter nu9
      nu10,      & ! CLUBB tunable parameter nu10
      nu_hm        ! CLUBB tunable parameter nu_hm

    ! Output Variables
    type(nu_vertical_res_dep), intent(out) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    call adj_low_res_nu( nzmax, grid_type, deltaz,  & ! Intent(in)
                         momentum_heights, thermodynamic_heights, & ! Intent(in)
                         l_prescribed_avg_deltaz, mult_coef, &  ! Intent(in)
                         nu1, nu2, nu6, nu8, nu9, nu10, nu_hm, &  ! Intent(in)
                         nu_vert_res_dep )  ! Intent(out)

  end subroutine adj_low_res_nu_api

! The CLUBB_CAM preprocessor directives are being commented out because this
! code is now also used for WRF-CLUBB.
!#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
  !================================================================================================
  ! pack_pdf_params - Returns a two dimensional real array with all values.
  !================================================================================================

  subroutine pack_pdf_params_api( pdf_params, nz, r_param_array, &
                                  k_start, k_end )

    use pdf_parameter_module, only : pack_pdf_params

    !use statements

    implicit none

    integer, intent(in) :: nz ! Num Vert Model Levs
    
    ! Input a pdf_parameter array with nz instances of pdf_parameter
    type (pdf_parameter), intent(inout) :: pdf_params

    ! Output a two dimensional real array with all values
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(inout) :: &
      r_param_array
      
    integer, optional, intent(in) :: k_start, k_end
      
    if( present( k_start ) .and. present( k_end ) ) then
        call pack_pdf_params( pdf_params, nz, & ! intent(in)
                              r_param_array, & ! intent(out)
                              k_start, k_end ) ! intent(in/optional)
    else 
        call pack_pdf_params( pdf_params, nz, & ! intent(in)
                              r_param_array ) ! intent(out)
    end if

  end subroutine pack_pdf_params_api

  !================================================================================================
  ! unpack_pdf_params - Returns a pdf_parameter array with nz instances of pdf_parameter.
  !================================================================================================

  subroutine unpack_pdf_params_api( r_param_array, nz, pdf_params, &
                                    k_start, k_end )

    use pdf_parameter_module, only : unpack_pdf_params

    implicit none
    
    integer, intent(in) :: nz ! Num Vert Model Levs
    
    ! Input a two dimensional real array with pdf values
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(in) :: &
      r_param_array

    ! Output a pdf_parameter array with nz instances of pdf_parameter
    type (pdf_parameter), intent(inout) :: pdf_params
    
    integer, optional, intent(in) :: k_start, k_end

    if( present( k_start ) .and. present( k_end ) ) then
        call unpack_pdf_params( r_param_array, nz, & ! intent(in)
                                pdf_params, & ! intent(inout)
                                k_start, k_end ) ! intent(in/optional)
    else  
        call unpack_pdf_params( r_param_array, nz, & ! intent(in)
                                pdf_params ) ! intent(inout)
    end if
    
    
  end subroutine unpack_pdf_params_api
  
!#endif
  !================================================================================================
  ! init_pdf_params - allocates arrays for pdf_params
  !================================================================================================
  subroutine init_pdf_params_api( nz, ngrdcol, pdf_params )
  
    use pdf_parameter_module, only : init_pdf_params
    
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      nz, &   ! Number of vertical grid levels    [-]
      ngrdcol ! Number of grid columns            [-]

    ! Output Variable(s)
    type(pdf_parameter), intent(out) :: &
      pdf_params    ! PDF parameters            [units vary]
    
    call init_pdf_params( nz, ngrdcol, & ! intent(in)
                          pdf_params ) ! intent(out)
    
  end subroutine init_pdf_params_api
  
  !================================================================================================
  ! copy_single_pdf_params_to_multi - copies values of a single column version of pdf_params
  !   to a multiple column version for a specified column.
  !
  ! NOTE: THIS SUBROUTINE IS INTENDED TO BE TEMPORARY AND SHOULD BECOME UNNECESSARY ONCE 
  !       CLUBB IS ABLE TO OPERATE OVER MULTIPLE COLUMNS.
  !       See https://github.com/larson-group/cam/issues/129#issuecomment-827944454
  !================================================================================================
  subroutine copy_single_pdf_params_to_multi( pdf_params_single, icol, &
                                              pdf_params_multi )
    
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      icol   ! Column number to copy to
      
    type(pdf_parameter), intent(in) :: &
      pdf_params_single  ! PDF parameters            [units vary]

    ! Output Variable(s)
    type(pdf_parameter), intent(inout) :: &
      pdf_params_multi  ! PDF parameters            [units vary]
      
    pdf_params_multi%w_1(icol,:)                  = pdf_params_single%w_1(1,:) 
    pdf_params_multi%w_2(icol,:)                  = pdf_params_single%w_2(1,:) 
    pdf_params_multi%varnce_w_1(icol,:)           = pdf_params_single%varnce_w_1(1,:) 
    pdf_params_multi%varnce_w_2(icol,:)           = pdf_params_single%varnce_w_2(1,:) 
    pdf_params_multi%rt_1(icol,:)                 = pdf_params_single%rt_1(1,:) 
    pdf_params_multi%rt_2(icol,:)                 = pdf_params_single%rt_2(1,:) 
    pdf_params_multi%varnce_rt_1(icol,:)          = pdf_params_single%varnce_rt_1(1,:) 
    pdf_params_multi%varnce_rt_2(icol,:)          = pdf_params_single%varnce_rt_2(1,:) 
    pdf_params_multi%thl_1(icol,:)                = pdf_params_single%thl_1(1,:) 
    pdf_params_multi%thl_2(icol,:)                = pdf_params_single%thl_2(1,:) 
    pdf_params_multi%varnce_thl_1(icol,:)         = pdf_params_single%varnce_thl_1(1,:) 
    pdf_params_multi%varnce_thl_2(icol,:)         = pdf_params_single%varnce_thl_2(1,:) 
    pdf_params_multi%corr_w_rt_1(icol,:)          = pdf_params_single%corr_w_rt_1(1,:) 
    pdf_params_multi%corr_w_rt_2(icol,:)          = pdf_params_single%corr_w_rt_2(1,:) 
    pdf_params_multi%corr_w_thl_1(icol,:)         = pdf_params_single%corr_w_thl_1(1,:) 
    pdf_params_multi%corr_w_thl_2(icol,:)         = pdf_params_single%corr_w_thl_2(1,:) 
    pdf_params_multi%corr_rt_thl_1(icol,:)        = pdf_params_single%corr_rt_thl_1(1,:) 
    pdf_params_multi%corr_rt_thl_2(icol,:)        = pdf_params_single%corr_rt_thl_2(1,:) 
    pdf_params_multi%alpha_thl(icol,:)            = pdf_params_single%alpha_thl(1,:) 
    pdf_params_multi%alpha_rt(icol,:)             = pdf_params_single%alpha_rt(1,:) 
    pdf_params_multi%crt_1(icol,:)                = pdf_params_single%crt_1(1,:) 
    pdf_params_multi%crt_2(icol,:)                = pdf_params_single%crt_2(1,:) 
    pdf_params_multi%cthl_1(icol,:)               = pdf_params_single%cthl_1(1,:) 
    pdf_params_multi%cthl_2(icol,:)               = pdf_params_single%cthl_2(1,:) 
    pdf_params_multi%chi_1(icol,:)                = pdf_params_single%chi_1(1,:) 
    pdf_params_multi%chi_2(icol,:)                = pdf_params_single%chi_2(1,:) 
    pdf_params_multi%stdev_chi_1(icol,:)          = pdf_params_single%stdev_chi_1(1,:) 
    pdf_params_multi%stdev_chi_2(icol,:)          = pdf_params_single%stdev_chi_2(1,:) 
    pdf_params_multi%stdev_eta_1(icol,:)          = pdf_params_single%stdev_eta_1(1,:) 
    pdf_params_multi%stdev_eta_2(icol,:)          = pdf_params_single%stdev_eta_2(1,:) 
    pdf_params_multi%covar_chi_eta_1(icol,:)      = pdf_params_single%covar_chi_eta_1(1,:) 
    pdf_params_multi%covar_chi_eta_2(icol,:)      = pdf_params_single%covar_chi_eta_2(1,:) 
    pdf_params_multi%corr_w_chi_1(icol,:)         = pdf_params_single%corr_w_chi_1(1,:) 
    pdf_params_multi%corr_w_chi_2(icol,:)         = pdf_params_single%corr_w_chi_2(1,:) 
    pdf_params_multi%corr_w_eta_1(icol,:)         = pdf_params_single%corr_w_eta_1(1,:) 
    pdf_params_multi%corr_w_eta_2(icol,:)         = pdf_params_single%corr_w_eta_2(1,:) 
    pdf_params_multi%corr_chi_eta_1(icol,:)       = pdf_params_single%corr_chi_eta_1(1,:) 
    pdf_params_multi%corr_chi_eta_2(icol,:)       = pdf_params_single%corr_chi_eta_2(1,:) 
    pdf_params_multi%rsatl_1(icol,:)              = pdf_params_single%rsatl_1(1,:) 
    pdf_params_multi%rsatl_2(icol,:)              = pdf_params_single%rsatl_2(1,:) 
    pdf_params_multi%rc_1(icol,:)                 = pdf_params_single%rc_1(1,:) 
    pdf_params_multi%rc_2(icol,:)                 = pdf_params_single%rc_2(1,:) 
    pdf_params_multi%cloud_frac_1(icol,:)         = pdf_params_single%cloud_frac_1(1,:) 
    pdf_params_multi%cloud_frac_2(icol,:)         = pdf_params_single%cloud_frac_2(1,:) 
    pdf_params_multi%mixt_frac(icol,:)            = pdf_params_single%mixt_frac(1,:) 
    pdf_params_multi%ice_supersat_frac_1(icol,:)  = pdf_params_single%ice_supersat_frac_1(1,:) 
    pdf_params_multi%ice_supersat_frac_2(icol,:)  = pdf_params_single%ice_supersat_frac_2(1,:) 
    
  end subroutine copy_single_pdf_params_to_multi
  
  !================================================================================================
  ! copy_multi_pdf_params_to_single - copies values of a multiple column version of pdf_params
  !   at a specified column to a single column version.
  !
  ! NOTE: THIS SUBROUTINE IS INTENDED TO BE TEMPORARY AND SHOULD BECOME UNNECESSARY ONCE 
  !       CLUBB IS ABLE TO OPERATE OVER MULTIPLE COLUMNS.
  !       See https://github.com/larson-group/cam/issues/129#issuecomment-827944454
  !================================================================================================
  subroutine copy_multi_pdf_params_to_single( pdf_params_multi, icol, &
                                              pdf_params_single )
    
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: &
      icol   ! Column number to copy to
      
    type(pdf_parameter), intent(in) :: &
      pdf_params_multi  ! PDF parameters            [units vary]

    ! Output Variable(s)
    type(pdf_parameter), intent(inout) :: &
      pdf_params_single   ! PDF parameters            [units vary]
      
    pdf_params_single%w_1(1,:)                  = pdf_params_multi%w_1(icol,:) 
    pdf_params_single%w_2(1,:)                  = pdf_params_multi%w_2(icol,:) 
    pdf_params_single%varnce_w_1(1,:)           = pdf_params_multi%varnce_w_1(icol,:) 
    pdf_params_single%varnce_w_2(1,:)           = pdf_params_multi%varnce_w_2(icol,:) 
    pdf_params_single%rt_1(1,:)                 = pdf_params_multi%rt_1(icol,:) 
    pdf_params_single%rt_2(1,:)                 = pdf_params_multi%rt_2(icol,:) 
    pdf_params_single%varnce_rt_1(1,:)          = pdf_params_multi%varnce_rt_1(icol,:) 
    pdf_params_single%varnce_rt_2(1,:)          = pdf_params_multi%varnce_rt_2(icol,:) 
    pdf_params_single%thl_1(1,:)                = pdf_params_multi%thl_1(icol,:) 
    pdf_params_single%thl_2(1,:)                = pdf_params_multi%thl_2(icol,:) 
    pdf_params_single%varnce_thl_1(1,:)         = pdf_params_multi%varnce_thl_1(icol,:) 
    pdf_params_single%varnce_thl_2(1,:)         = pdf_params_multi%varnce_thl_2(icol,:) 
    pdf_params_single%corr_w_rt_1(1,:)          = pdf_params_multi%corr_w_rt_1(icol,:) 
    pdf_params_single%corr_w_rt_2(1,:)          = pdf_params_multi%corr_w_rt_2(icol,:) 
    pdf_params_single%corr_w_thl_1(1,:)         = pdf_params_multi%corr_w_thl_1(icol,:) 
    pdf_params_single%corr_w_thl_2(1,:)         = pdf_params_multi%corr_w_thl_2(icol,:) 
    pdf_params_single%corr_rt_thl_1(1,:)        = pdf_params_multi%corr_rt_thl_1(icol,:) 
    pdf_params_single%corr_rt_thl_2(1,:)        = pdf_params_multi%corr_rt_thl_2(icol,:) 
    pdf_params_single%alpha_thl(1,:)            = pdf_params_multi%alpha_thl(icol,:) 
    pdf_params_single%alpha_rt(1,:)             = pdf_params_multi%alpha_rt(icol,:) 
    pdf_params_single%crt_1(1,:)                = pdf_params_multi%crt_1(icol,:) 
    pdf_params_single%crt_2(1,:)                = pdf_params_multi%crt_2(icol,:) 
    pdf_params_single%cthl_1(1,:)               = pdf_params_multi%cthl_1(icol,:) 
    pdf_params_single%cthl_2(1,:)               = pdf_params_multi%cthl_2(icol,:) 
    pdf_params_single%chi_1(1,:)                = pdf_params_multi%chi_1(icol,:) 
    pdf_params_single%chi_2(1,:)                = pdf_params_multi%chi_2(icol,:) 
    pdf_params_single%stdev_chi_1(1,:)          = pdf_params_multi%stdev_chi_1(icol,:) 
    pdf_params_single%stdev_chi_2(1,:)          = pdf_params_multi%stdev_chi_2(icol,:) 
    pdf_params_single%stdev_eta_1(1,:)          = pdf_params_multi%stdev_eta_1(icol,:) 
    pdf_params_single%stdev_eta_2(1,:)          = pdf_params_multi%stdev_eta_2(icol,:) 
    pdf_params_single%covar_chi_eta_1(1,:)      = pdf_params_multi%covar_chi_eta_1(icol,:) 
    pdf_params_single%covar_chi_eta_2(1,:)      = pdf_params_multi%covar_chi_eta_2(icol,:) 
    pdf_params_single%corr_w_chi_1(1,:)         = pdf_params_multi%corr_w_chi_1(icol,:) 
    pdf_params_single%corr_w_chi_2(1,:)         = pdf_params_multi%corr_w_chi_2(icol,:) 
    pdf_params_single%corr_w_eta_1(1,:)         = pdf_params_multi%corr_w_eta_1(icol,:) 
    pdf_params_single%corr_w_eta_2(1,:)         = pdf_params_multi%corr_w_eta_2(icol,:) 
    pdf_params_single%corr_chi_eta_1(1,:)       = pdf_params_multi%corr_chi_eta_1(icol,:) 
    pdf_params_single%corr_chi_eta_2(1,:)       = pdf_params_multi%corr_chi_eta_2(icol,:) 
    pdf_params_single%rsatl_1(1,:)              = pdf_params_multi%rsatl_1(icol,:) 
    pdf_params_single%rsatl_2(1,:)              = pdf_params_multi%rsatl_2(icol,:) 
    pdf_params_single%rc_1(1,:)                 = pdf_params_multi%rc_1(icol,:) 
    pdf_params_single%rc_2(1,:)                 = pdf_params_multi%rc_2(icol,:) 
    pdf_params_single%cloud_frac_1(1,:)         = pdf_params_multi%cloud_frac_1(icol,:) 
    pdf_params_single%cloud_frac_2(1,:)         = pdf_params_multi%cloud_frac_2(icol,:) 
    pdf_params_single%mixt_frac(1,:)            = pdf_params_multi%mixt_frac(icol,:) 
    pdf_params_single%ice_supersat_frac_1(1,:)  = pdf_params_multi%ice_supersat_frac_1(icol,:) 
    pdf_params_single%ice_supersat_frac_2(1,:)  = pdf_params_multi%ice_supersat_frac_2(icol,:) 
    
  end subroutine copy_multi_pdf_params_to_single

  
  !================================================================================================
  ! init_precip_fracs - allocates arrays for precip_fracs
  !================================================================================================
  subroutine init_precip_fracs_api( nz, ngrdcol, &
                                    precip_fracs )

    use hydromet_pdf_parameter_module, only : init_precip_fracs

    implicit none
    
    ! Input Variable(s)
    integer, intent(in) :: &
      nz,     & ! Number of vertical grid levels    [-]
      ngrdcol   ! Number of grid columns            [-]

    ! Output Variable
    type(precipitation_fractions), intent(out) :: &
      precip_fracs    ! Hydrometeor PDF parameters      [units vary]

    call init_precip_fracs( nz, ngrdcol, & ! intent(in)
                            precip_fracs ) ! intent(out)

    return

  end subroutine init_precip_fracs_api
  
  !================================================================================================
  ! init_pdf_implicit_coefs_terms - allocates arrays for the PDF implicit
  ! coefficient and explicit terms.
  !================================================================================================
  subroutine init_pdf_implicit_coefs_terms_api( nz, sclr_dim, &
                                                pdf_implicit_coefs_terms )

    use pdf_parameter_module, only: &
        init_pdf_implicit_coefs_terms    ! Procedure(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,       & ! Number of vertical grid levels    [-]
      sclr_dim    ! Number of scalar variables        [-]

    ! Output Variable
    type(implicit_coefs_terms), intent(out) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    call init_pdf_implicit_coefs_terms( nz, sclr_dim, & ! intent(in)
                                        pdf_implicit_coefs_terms ) ! intent(out)

  end subroutine init_pdf_implicit_coefs_terms_api

  !================================================================================================
  ! setup_pdf_parameters
  !================================================================================================

  subroutine setup_pdf_parameters_api_single_col( gr, & ! intent(in)
    nz, pdf_dim, dt, &                      ! Intent(in)
    Nc_in_cloud, rcm, cloud_frac, Kh_zm, &      ! Intent(in)
    ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
    corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
    pdf_params, l_stats_samp, &                 ! Intent(in)
    clubb_params, &                             ! Intent(in)
    iiPDF_type, &                               ! Intent(in)
    l_use_precip_frac, &                        ! Intent(in)
    l_predict_upwp_vpwp, &                      ! Intent(in)
    l_diagnose_correlations, &                  ! Intent(in)
    l_calc_w_corr, &                            ! Intent(in)
    l_const_Nc_in_cloud, &                      ! Intent(in)
    l_fix_w_chi_eta_correlations, &             ! Intent(in)
    stats_zt, stats_zm, stats_sfc, &            ! intent(inout)
    hydrometp2, &                               ! Intent(inout)
    mu_x_1_n, mu_x_2_n, &                       ! Intent(out)
    sigma_x_1_n, sigma_x_2_n, &                 ! Intent(out)
    corr_array_1_n, corr_array_2_n, &           ! Intent(out)
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
    precip_fracs,  &                            ! Intent(out)
    hydromet_pdf_params )                       ! Intent(out)

    use setup_clubb_pdf_params, only : setup_pdf_parameters

    use parameter_indices, only: &
        nparams    ! Variable(s)

    use error_code, only : &
        err_code, &         ! Error Indicator
        clubb_fatal_error   ! Constant

    

    implicit none

    type(stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    type(grid), target, intent(in) :: gr

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      pdf_dim   ! Number of variables in the correlation array

    real( kind = core_rknd ), intent(in) ::  &
      dt    ! Model timestep                                           [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      rcm,               & ! Mean cloud water mixing ratio, < r_c >    [kg/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      Kh_zm,             & ! Eddy diffusivity coef. on momentum levels [m^2/s]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), &
      intent(in) :: &
      corr_array_n_cloud, & ! Prescribed norm. space corr. array in cloud    [-]
      corr_array_n_below    ! Prescribed norm. space corr. array below cloud [-]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations    ! Use a fixed correlation for s and t Mellor(chi/eta)

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(inout) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    ! Output Variables
    real( kind = core_rknd ), dimension(nz,pdf_dim), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(nz,pdf_dim,pdf_dim), &
      intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space):  PDF vars. (comp. 1)   [-]
      corr_array_2_n    ! Corr. array (normal space):  PDF vars. (comp. 2)   [-]

    real( kind = core_rknd ), dimension(nz,pdf_dim,pdf_dim), &
      intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    ! This is only an output, but it contains allocated arrays, so we need to treat it as inout
    type(precipitation_fractions), intent(inout) :: &
      precip_fracs           ! Precipitation fractions      [-]

    type(hydromet_pdf_parameter), dimension(nz), optional, intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]
      
    ! -------------- Local Variables --------------
    
    real( kind = core_rknd ), dimension(1,nz) :: &
      Nc_in_cloud_col,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      rcm_col,               & ! Mean cloud water mixing ratio, < r_c >    [kg/kg]
      cloud_frac_col,        & ! Cloud fraction                            [-]
      Kh_zm_col,             & ! Eddy diffusivity coef. on momentum levels [m^2/s]
      ice_supersat_frac_col    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(1,nz,hydromet_dim) :: &
      hydromet_col,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp_col    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]
      
    ! Input/Output Variables
    real( kind = core_rknd ), dimension(1,nz,hydromet_dim) :: &
      hydrometp2_col    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    ! Output Variables
    real( kind = core_rknd ), dimension(1,nz,pdf_dim) :: &
      mu_x_1_n_col,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n_col,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n_col, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n_col    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(1,nz,pdf_dim,pdf_dim) :: &
      corr_array_1_n_col, & ! Corr. array (normal space):  PDF vars. (comp. 1)   [-]
      corr_array_2_n_col    ! Corr. array (normal space):  PDF vars. (comp. 2)   [-]

    real( kind = core_rknd ), dimension(1,nz,pdf_dim,pdf_dim) :: &
      corr_cholesky_mtx_1_col, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2_col    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    type(hydromet_pdf_parameter), dimension(1,nz) :: &
      hydromet_pdf_params_col    ! Hydrometeor PDF parameters        [units vary]
      
    type(stats), dimension(1) :: &
      stats_zt_col, &
      stats_zm_col, &
      stats_sfc_col

    type(grid), target, dimension(1) :: &
      gr_col


    Nc_in_cloud_col(1,:) = Nc_in_cloud
    rcm_col(1,:) = rcm
    cloud_frac_col(1,:) = cloud_frac
    Kh_zm_col(1,:) = Kh_zm
    ice_supersat_frac_col(1,:) = ice_supersat_frac
    
    hydromet_col(1,:,:) = hydromet
    wphydrometp_col(1,:,:) = wphydrometp
    
    gr_col(1) = gr
    
    stats_zt_col(1) = stats_zt
    stats_zm_col(1) = stats_zm
    stats_sfc_col(1) = stats_sfc

    call setup_pdf_parameters( gr_col, &                          ! intent(in)
      nz, 1, pdf_dim, dt, &                                   ! Intent(in)
      Nc_in_cloud_col, rcm_col, cloud_frac_col, Kh_zm_col, &  ! Intent(in)
      ice_supersat_frac_col, hydromet_col, wphydrometp_col, & ! Intent(in)
      corr_array_n_cloud, corr_array_n_below, &               ! Intent(in)
      pdf_params, l_stats_samp, &                             ! Intent(in)
      clubb_params, &                                         ! Intent(in)
      iiPDF_type, &                                           ! Intent(in)
      l_use_precip_frac, &                                    ! Intent(in)
      l_predict_upwp_vpwp, &                                  ! Intent(in)
      l_diagnose_correlations, &                              ! Intent(in)
      l_calc_w_corr, &                                        ! Intent(in)
      l_const_Nc_in_cloud, &                                  ! Intent(in)
      l_fix_w_chi_eta_correlations, &                         ! Intent(in)
      stats_zt_col, stats_zm_col, stats_sfc_col, &            ! Intent(inout)
      hydrometp2_col, &                                       ! Intent(inout)
      mu_x_1_n_col, mu_x_2_n_col, &                           ! Intent(out)
      sigma_x_1_n_col, sigma_x_2_n_col, &                     ! Intent(out)
      corr_array_1_n_col, corr_array_2_n_col, &               ! Intent(out)
      corr_cholesky_mtx_1_col, corr_cholesky_mtx_2_col, &     ! Intent(out)
      precip_fracs, &                                         ! Intent(inout)
      hydromet_pdf_params_col )                               ! Optional(out)

    if ( err_code == clubb_fatal_error ) error stop
    
    ! The following does not work for stats 
    !     stats_zt = stats_zt_col(1)
    !     stats_zm = stats_zm_col(1) 
    !     stats_sfc = stats_sfc_col(1)
    ! because of some mysterious pointer issue. However, the only thing that 
    ! updates in stats is the field values, so we can copy only those instead.
    if ( l_stats ) then 
      stats_zm%accum_field_values = stats_zm_col(1)%accum_field_values
      stats_zm%accum_num_samples = stats_zm_col(1)%accum_num_samples
      
      stats_zt%accum_field_values = stats_zt_col(1)%accum_field_values
      stats_zt%accum_num_samples = stats_zt_col(1)%accum_num_samples
      
      stats_sfc%accum_field_values = stats_sfc_col(1)%accum_field_values
      stats_sfc%accum_num_samples = stats_sfc_col(1)%accum_num_samples
    end if
    
    hydrometp2 = hydrometp2_col(1,:,:)
    mu_x_1_n = mu_x_1_n_col(1,:,:)
    mu_x_2_n = mu_x_2_n_col(1,:,:)
    sigma_x_1_n = sigma_x_1_n_col(1,:,:)
    sigma_x_2_n = sigma_x_2_n_col(1,:,:)
    corr_array_1_n = corr_array_1_n_col(1,:,:,:)
    corr_array_2_n = corr_array_2_n_col(1,:,:,:)
    corr_cholesky_mtx_1 = corr_cholesky_mtx_1_col(1,:,:,:)
    corr_cholesky_mtx_2 = corr_cholesky_mtx_2_col(1,:,:,:)
    if ( present(hydromet_pdf_params) ) then
      hydromet_pdf_params = hydromet_pdf_params_col(1,:)
    end if

  end subroutine setup_pdf_parameters_api_single_col
!===========================================================================! 
  subroutine setup_pdf_parameters_api_multi_col( gr, &
    nz, ngrdcol, pdf_dim, dt, &                 ! Intent(in)
    Nc_in_cloud, rcm, cloud_frac, Kh_zm, &      ! Intent(in)
    ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
    corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
    pdf_params, l_stats_samp, &                 ! Intent(in)
    clubb_params, &                             ! Intent(in)
    iiPDF_type, &                               ! Intent(in)
    l_use_precip_frac, &                        ! Intent(in)
    l_predict_upwp_vpwp, &                      ! Intent(in)
    l_diagnose_correlations, &                  ! Intent(in)
    l_calc_w_corr, &                            ! Intent(in)
    l_const_Nc_in_cloud, &                      ! Intent(in)
    l_fix_w_chi_eta_correlations, &             ! Intent(in)
    stats_zt, stats_zm, stats_sfc, &            ! intent(inout)
    hydrometp2, &                               ! Intent(inout)
    mu_x_1_n, mu_x_2_n, &                       ! Intent(out)
    sigma_x_1_n, sigma_x_2_n, &                 ! Intent(out)
    corr_array_1_n, corr_array_2_n, &           ! Intent(out)
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
    precip_fracs, &                             ! Intent(out)
    hydromet_pdf_params )                       ! Intent(out)

    use setup_clubb_pdf_params, only : setup_pdf_parameters

    use parameter_indices, only: &
        nparams    ! Variable(s)

    use error_code, only : &
        err_code, &         ! Error Indicator
        clubb_fatal_error   ! Constant

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      pdf_dim,     & ! Number of variables in the correlation array
      ngrdcol        ! Number of grid columns
      
    type(grid), target, dimension(ngrdcol), intent(in) :: gr

    real( kind = core_rknd ), intent(in) ::  &
      dt    ! Model timestep                                           [s]

    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      rcm,               & ! Mean cloud water mixing ratio, < r_c >    [kg/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      Kh_zm,             & ! Eddy diffusivity coef. on momentum levels [m^2/s]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(in) :: &
      hydromet,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(pdf_dim,pdf_dim), &
      intent(in) :: &
      corr_array_n_cloud, & ! Prescribed norm. space corr. array in cloud    [-]
      corr_array_n_below    ! Prescribed norm. space corr. array below cloud [-]

    type(pdf_parameter), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    integer, intent(in) :: &
      iiPDF_type    ! Selected option for the two-component normal (double
                    ! Gaussian) PDF type to use for the w, rt, and theta-l (or
                    ! w, chi, and eta) portion of CLUBB's multivariate,
                    ! two-component PDF.

    logical, intent(in) :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations    ! Use a fixed correlation for s and t Mellor(chi/eta)

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz,hydromet_dim), intent(inout) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    type(stats), target, dimension(ngrdcol), intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim), &
      intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space):  PDF vars. (comp. 1)   [-]
      corr_array_2_n    ! Corr. array (normal space):  PDF vars. (comp. 2)   [-]

    real( kind = core_rknd ), dimension(ngrdcol,nz,pdf_dim,pdf_dim), &
      intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    type(hydromet_pdf_parameter), dimension(ngrdcol,nz), optional, intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]
      
    ! This is only an output, but it contains allocated arrays, so we need to treat it as inout
    type(precipitation_fractions), intent(inout) :: &
      precip_fracs           ! Precipitation fractions      [-]

    call setup_pdf_parameters( gr, &              ! intent(in)
      nz, ngrdcol, pdf_dim, dt, &                 ! Intent(in)
      Nc_in_cloud, rcm, cloud_frac, Kh_zm, &      ! Intent(in)
      ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
      corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
      pdf_params, l_stats_samp, &                 ! Intent(in)
      clubb_params, &                             ! Intent(in)
      iiPDF_type, &                               ! Intent(in)
      l_use_precip_frac, &                        ! Intent(in)
      l_predict_upwp_vpwp, &                      ! Intent(in)
      l_diagnose_correlations, &                  ! Intent(in)
      l_calc_w_corr, &                            ! Intent(in)
      l_const_Nc_in_cloud, &                      ! Intent(in)
      l_fix_w_chi_eta_correlations, &             ! Intent(in)
      stats_zt, stats_zm, stats_sfc, &            ! intent(inout)
      hydrometp2, &                               ! Intent(inout)
      mu_x_1_n, mu_x_2_n, &                       ! Intent(out)
      sigma_x_1_n, sigma_x_2_n, &                 ! Intent(out)
      corr_array_1_n, corr_array_2_n, &           ! Intent(out)
      corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
      precip_fracs, &                             ! Intent(inout)
      hydromet_pdf_params )                       ! Optional(out)

    if ( err_code == clubb_fatal_error ) error stop

  end subroutine setup_pdf_parameters_api_multi_col

  !================================================================================================
  ! stats_init - Initializes the statistics saving functionality of the CLUBB model.
  !================================================================================================

  subroutine stats_init_api( &
    iunit, fname_prefix, fdir, l_stats_in, &
    stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
    nzmax, nlon, nlat, gzt, gzm, nnrad_zt, &
    grad_zt, nnrad_zm, grad_zm, day, month, year, &
    lon_vals, lat_vals, time_current, delt, l_silhs_out_in, &
    stats_zt, stats_zm, stats_sfc, &
    stats_lh_zt, stats_lh_sfc, &
    stats_rad_zt, stats_rad_zm )

    use stats_clubb_utilities, only : stats_init

    implicit none

    type(stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    ! Input Variables
    integer, intent(in) :: iunit  ! File unit for fnamelist

    character(len=*), intent(in) ::  &
      fname_prefix, & ! Start of the stats filenames
      fdir            ! Directory to output to

    logical, intent(in) :: &
      l_stats_in      ! Stats on? T/F

    character(len=*), intent(in) :: &
      stats_fmt_in    ! Format of the stats file output

    real( kind = core_rknd ), intent(in) ::  &
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    character(len=*), intent(in) :: &
      fnamelist          ! Filename holding the &statsnl

    integer, intent(in) :: &
      nlon, & ! Number of points in the X direction [-]
      nlat, & ! Number of points in the Y direction [-]
      nzmax   ! Grid points in the vertical         [-]

    real( kind = core_rknd ), intent(in), dimension(nzmax) ::  &
      gzt, gzm  ! Thermodynamic and momentum levels           [m]

    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zt) :: grad_zt ! Radiation levels [m]

    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real( kind = core_rknd ), intent(in), dimension(nnrad_zm) :: grad_zm ! Radiation levels [m]

    integer, intent(in) :: day, month, year  ! Time of year

    real( kind = core_rknd ), dimension(nlon), intent(in) ::  &
      lon_vals  ! Longitude values [Degrees E]

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  &
      lat_vals  ! Latitude values  [Degrees N]

    real( kind = time_precision ), intent(in) ::  &
      time_current ! Model time                         [s]

    real( kind = core_rknd ), intent(in) ::  &
      delt         ! Timestep (dt_main in CLUBB)         [s]

    logical, intent(in) :: &
      l_silhs_out_in  ! Whether to output SILHS files (stats_lh_zt,stats_lh_sfc) [dimensionless]

    call stats_init( &
      iunit, fname_prefix, fdir, l_stats_in, & ! intent(in)
      stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, & ! intent(in)
      nzmax, nlon, nlat, gzt, gzm, nnrad_zt, & ! intent(in)
      grad_zt, nnrad_zm, grad_zm, day, month, year, & ! intent(in)
      lon_vals, lat_vals, time_current, delt, l_silhs_out_in, & ! intent(in)
      stats_zt, stats_zm, stats_sfc, & ! intent(inout)
      stats_lh_zt, stats_lh_sfc, & ! intent(inout)
      stats_rad_zt, stats_rad_zm ) ! intent(inout)

    if ( err_code == clubb_fatal_error ) error stop
    
  end subroutine stats_init_api

  !================================================================================================
  ! stats_begin_timestep - Sets flags determining specific timestep info.
  !================================================================================================

  subroutine stats_begin_timestep_api( &
    itime, stats_nsamp, stats_nout )


    use stats_clubb_utilities, only : stats_begin_timestep

    implicit none

    ! External
    intrinsic :: mod

    ! Input Variable(s)
    integer, intent(in) ::  &
      itime,        & ! Elapsed model time       [timestep]
      stats_nsamp,  & ! Stats sampling interval  [timestep]
      stats_nout      ! Stats output interval    [timestep]

    call stats_begin_timestep( &
      itime, stats_nsamp, stats_nout ) ! intent(in)
  end subroutine stats_begin_timestep_api

  !================================================================================================
  ! stats_end_timestep - Calls statistics to be written to the output format.
  !================================================================================================

  subroutine stats_end_timestep_api( clubb_params, &
                                     stats_zt, stats_zm, stats_sfc, &
                                     stats_lh_zt, stats_lh_sfc, &
                                     stats_rad_zt, stats_rad_zm &
#ifdef NETCDF
                                     , l_uv_nudge, &
                                     l_tke_aniso, &
                                     l_standard_term_ta &
#endif
                                      )

    use stats_clubb_utilities, only : stats_end_timestep

    use parameter_indices, only: nparams


    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

    type(stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

#ifdef NETCDF
    logical, intent(in) :: &
      l_uv_nudge,         & ! For wind speed nudging.
      l_tke_aniso,        & ! For anisotropic turbulent kinetic energy, i.e.
                            ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_standard_term_ta    ! Use the standard discretization for the turbulent advection terms.
                            ! Setting to .false. means that a_1 and a_3 are pulled outside of the
                            ! derivative in advance_wp2_wp3_module.F90 and in
                            ! advance_xp2_xpyp_module.F90.
#endif

    call stats_end_timestep( clubb_params, &              ! intent(in)
                             stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                             stats_lh_zt, stats_lh_sfc, & ! intent(inout)
                             stats_rad_zt, stats_rad_zm & ! intent(inout)
#ifdef NETCDF
                             , l_uv_nudge, & ! Intent(in)
                             l_tke_aniso, & ! Intent(in)
                             l_standard_term_ta & ! Intent(in)
#endif
                              )

    if ( err_code == clubb_fatal_error ) error stop

  end subroutine stats_end_timestep_api

  !================================================================================================
  ! stats_accumulate_hydromet - Computes stats related the hydrometeors.
  !================================================================================================

  subroutine stats_accumulate_hydromet_api( gr, &
                                            hydromet, rho_ds_zt, &
                                            stats_zt, stats_sfc )

    use stats_clubb_utilities, only : stats_accumulate_hydromet

    

    implicit none

    type(stats), target, intent(inout) :: &
      stats_zt, &
      stats_sfc

    type(grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet ! All hydrometeors except for rcm        [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho_ds_zt ! Dry, static density (thermo. levs.)      [kg/m^3]

    call stats_accumulate_hydromet( gr, & ! intent(in)
      hydromet, rho_ds_zt, & ! intent(in)
      stats_zt, stats_sfc ) ! intent(inout)
  end subroutine stats_accumulate_hydromet_api

  !================================================================================================
  ! stats_finalize - Close NetCDF files and deallocate scratch space and stats file structures.
  !================================================================================================

  subroutine stats_finalize_api ( stats_zt, stats_zm, stats_sfc, &
                                  stats_lh_zt, stats_lh_sfc, &
                                  stats_rad_zt, stats_rad_zm )

    use stats_clubb_utilities, only : stats_finalize

    implicit none

    type(stats), target, intent(inout) :: &
      stats_zt, &
      stats_zm, &
      stats_sfc, &
      stats_lh_zt, &
      stats_lh_sfc, &
      stats_rad_zt, &
      stats_rad_zm

    call stats_finalize ( stats_zt, stats_zm, stats_sfc, & ! intent(inout)
                          stats_lh_zt, stats_lh_sfc, & ! intent(inout)
                          stats_rad_zt, stats_rad_zm ) ! intent(inout)

  end subroutine stats_finalize_api

  !================================================================================================
  ! stats_init_rad_zm - Initializes array indices for rad_zm variables.
  !================================================================================================

  subroutine stats_init_rad_zm_api( &
                                    vars_rad_zm, l_error, &
                                    stats_rad_zm )

    use stats_rad_zm_module, only : stats_init_rad_zm, nvarmax_rad_zm

    implicit none

    type(stats), target, intent(inout) :: &
      stats_rad_zm

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zm), intent(in) :: vars_rad_zm

    ! Input/Output Variable
    logical, intent(inout) :: l_error

    call stats_init_rad_zm( &
      vars_rad_zm, & ! intent(in)
      l_error, & ! intent(inout)
      stats_rad_zm ) ! intent(inout)
  end subroutine stats_init_rad_zm_api

  !================================================================================================
  ! stats_init_rad_zt - Initializes array indices for zt.
  !================================================================================================

  subroutine stats_init_rad_zt_api( &
                                    vars_rad_zt, l_error, &
                                    stats_rad_zt )

    use stats_rad_zt_module, only : stats_init_rad_zt, nvarmax_rad_zt

    implicit none

    type(stats), target, intent(inout) :: &
      stats_rad_zt

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zt), intent(in) :: vars_rad_zt

    ! Input/Output Variable
    logical, intent(inout) :: l_error

    call stats_init_rad_zt( &
      vars_rad_zt, & ! intent(in)
      l_error, & ! intent(inout)
      stats_rad_zt ) ! intent(inout)

  end subroutine stats_init_rad_zt_api

  !================================================================================================
  ! stats_init_zm - Initializes array indices for zm.
  !================================================================================================

  subroutine stats_init_zm_api( &
                                vars_zm, l_error, &
                                stats_zm )

    use stats_zm_module, only : stats_init_zm, nvarmax_zm

    implicit none

    type(stats), target, intent(inout) :: &
      stats_zm

    ! Input Variable
    character(len= * ), dimension(nvarmax_zm), intent(in) :: vars_zm ! zm variable names

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    call stats_init_zm( &
      vars_zm, & ! intent(in)
      l_error, & ! intent(inout)
      stats_zm ) ! intent(inout)

  end subroutine stats_init_zm_api

  !================================================================================================
  ! stats_init_zt - Initializes array indices for zt.
  !================================================================================================

  subroutine stats_init_zt_api( &
                                vars_zt, l_error, &
                                stats_zt )

    use stats_zt_module, only : stats_init_zt, nvarmax_zt

    implicit none

    type(stats), target, intent(inout) :: &
      stats_zt

    ! Input Variable
    character(len= * ), dimension(nvarmax_zt), intent(in) :: vars_zt

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    call stats_init_zt( &
      vars_zt, & ! intent(in)
      l_error, & ! intent(inout)
      stats_zt ) ! intent(inout)

  end subroutine stats_init_zt_api

  !================================================================================================
  ! stats_init_sfc - Initializes array indices for sfc.
  !================================================================================================

  subroutine stats_init_sfc_api( &
                                 vars_sfc, l_error, &
                                 stats_sfc )

    use stats_sfc_module, only : stats_init_sfc, nvarmax_sfc

    implicit none

    type(stats), target, intent(inout) :: &
      stats_sfc

    ! Input Variable
    character(len= * ), dimension(nvarmax_sfc), intent(in) :: vars_sfc

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    call stats_init_sfc( &
      vars_sfc, & ! intent(in)
      l_error, & ! intent(inout)
      stats_sfc ) ! intent(inout)

  end subroutine stats_init_sfc_api

  !================================================================================================
  ! thlm2T_in_K - Calculates absolute temperature from liquid water potential temperature.
  !================================================================================================

  elemental function thlm2T_in_K_api( &
    thlm, exner, rcm )  &
    result( T_in_K )

    use T_in_K_module, only : thlm2T_in_K

    implicit none

    ! Input
    real( kind = core_rknd ), intent(in) :: &
      thlm,   & ! Liquid potential temperature  [K]
      exner,  & ! Exner function                [-]
      rcm       ! Liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ) :: &
      T_in_K ! Result temperature [K]

    T_in_K = thlm2T_in_K( &
      thlm, exner, rcm )

  end function thlm2T_in_K_api

  !================================================================================================
  ! T_in_K2thlm - Calculates liquid water potential temperature from absolute temperature
  !================================================================================================

  elemental function T_in_K2thlm_api( &
    T_in_K, exner, rcm )  &
    result( thlm )

    use T_in_K_module, only : T_in_K2thlm

    implicit none

    ! Input
    real( kind = core_rknd ), intent(in) :: &
      T_in_K, &! Result temperature [K]
      exner,  & ! Exner function                [-]
      rcm       ! Liquid water mixing ratio     [kg/kg]

    real( kind = core_rknd ) :: &
      thlm    ! Liquid potential temperature  [K]

    thlm = T_in_K2thlm( &
      T_in_K, exner, rcm )

  end function T_in_K2thlm_api

  !================================================================================================
  ! calculate_spurious_source - Checks whether there is conservation within the column.
  !================================================================================================
  function calculate_spurious_source_api ( &
    integral_after, integral_before, &
    flux_top, flux_sfc, &
    integral_forcing, dt ) result( spurious_source )

    use numerical_check, only : calculate_spurious_source

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      integral_after, &   ! Vertically-integrated quantity after dt time  [units vary]
      integral_before, &  ! Vertically-integrated quantity before dt time [units vary]
      flux_top, &         ! Total flux at the top of the domain           [units vary]
      flux_sfc, &         ! Total flux at the bottom of the domain        [units vary]
      integral_forcing, & ! Vertically-integrated forcing                 [units vary]
      dt                  ! Timestep size                                 [s]

    ! Return Variable
    real( kind = core_rknd ) :: spurious_source ! [units vary]

    spurious_source = calculate_spurious_source( &
      integral_after, integral_before, &
      flux_top, flux_sfc, &
      integral_forcing, dt )

  end function calculate_spurious_source_api

  !================================================================================================
  ! zm2zt_scalar - Interpolates a variable from zm to zt grid at one height level
  !================================================================================================
  function zm2zt_scalar_api( gr, azm, k )

    use grid_class, only: & 
        grid, & ! Type
        zm2zt

    implicit none

    type(grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    integer, intent(in) :: &
      k      ! Vertical level index

    ! Return Variable
    real( kind = core_rknd ) :: &
      zm2zt_scalar_api   ! Variable when interp. to thermo. levels

    zm2zt_scalar_api = zm2zt( gr, azm, k )

  end function zm2zt_scalar_api

  !================================================================================================
  ! zt2zm_scalar - Interpolates a variable from zt to zm grid at one height level
  !================================================================================================
  function zt2zm_scalar_api( gr, azt, k )

    use grid_class, only: &
        grid, & ! Type
        zt2zm

    implicit none

    type(grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    integer, intent(in) :: &
      k      ! Vertical level index

    ! Return Variable
    real( kind = core_rknd ) :: &
      zt2zm_scalar_api   ! Variable when interp. to momentum levels

    zt2zm_scalar_api = zt2zm( gr, azt, k )

  end function zt2zm_scalar_api

  !================================================================================================
  ! zt2zm_prof - Interpolates a variable (profile) from zt to zm grid
  !================================================================================================
  function zt2zm_prof_api( gr, azt )

    use grid_class, only: &
        grid, & ! Type
        zt2zm

    implicit none

    type(grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      zt2zm_prof_api   ! Variable when interp. to momentum levels

    zt2zm_prof_api = zt2zm( gr, azt )

  end function zt2zm_prof_api

  !================================================================================================
  ! zm2zt_prof - Interpolates a variable (profile) from zm to zt grid
  !================================================================================================
  function zm2zt_prof_api( gr, azm )

    use grid_class, only: &
        grid, & ! Type
        zm2zt

    implicit none

    type(grid), target, intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      zm2zt_prof_api   ! Variable when interp. to thermo. levels

    zm2zt_prof_api = zm2zt( gr, azm )

  end function zm2zt_prof_api
  
  !================================================================================================
  ! zm2zt_2D - Interpolates a variable (profile) from zm to zt grid for all columns
  !================================================================================================
  function zm2zt_2D_api( nz, ngrdcol, gr, azm )

    use grid_class, only: &
        grid, & ! Type
        zm2zt

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type(grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      zm2zt_2D_api   ! Variable when interp. to thermo. levels

    zm2zt_2D_api = zm2zt( nz, ngrdcol, gr, azm )

  end function zm2zt_2D_api
  
  !================================================================================================
  ! zt2zm_2D - Interpolates a variable (profile) from zt to zm grid for all columns
  !================================================================================================
  function zt2zm_2D_api( nz, ngrdcol, gr, azm )

    use grid_class, only: &
        grid, & ! Type
        zt2zm

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    type(grid), target, dimension(ngrdcol), intent(in) :: gr

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      zt2zm_2D_api   ! Variable when interp. to thermo. levels

    zt2zm_2D_api = zt2zm( nz, ngrdcol, gr, azm )

  end function zt2zm_2D_api

  !================================================================================================
  ! calculate_thlp2_rad - Computes the contribution of radiative cooling to thlp2
  !================================================================================================
  pure subroutine calculate_thlp2_rad_api &
                  ( nz, rcm_zm, thlprcp, radht_zm, & ! Intent(in)
                    clubb_params,                  & ! Intent(in)
                    thlp2_forcing )                  ! Intent(inout)

    use clubb_precision, only: &
        core_rknd                     ! Constant(s)

    use advance_clubb_core_module, only: &
        calculate_thlp2_rad

    use parameter_indices, only: &
        nparams

    implicit none

  ! Input Variables
    integer, intent(in) :: &
      nz                    ! Number of vertical levels                      [-]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      rcm_zm, &             ! Cloud water mixing ratio on momentum grid      [kg/kg]
      thlprcp, &            ! thl'rc'                                        [K kg/kg]
      radht_zm              ! SW + LW heating rate (on momentum grid)        [K/s]

    real( kind = core_rknd ), dimension(nparams), intent(in) :: &
      clubb_params    ! Array of CLUBB's tunable parameters    [units vary]

  ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      thlp2_forcing         ! <th_l'^2> forcing (momentum levels)            [K^2/s]
  !----------------------------------------------------------------------

    call calculate_thlp2_rad( nz, rcm_zm, thlprcp, radht_zm, & ! intent(in)
                              clubb_params,                  & ! intent(in)
                              thlp2_forcing )                  ! intent(inout)

    return
  end subroutine calculate_thlp2_rad_api

  !================================================================================================
  ! update_xp2_mc - Calculates the effects of rain evaporation on rtp2 and thlp2
  !================================================================================================
  subroutine update_xp2_mc_api( gr, nz, dt, cloud_frac, rcm, rvm, thlm,        &
                            wm, exner, rrm_evap, pdf_params,        &
                            rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc,    &
                            rtpthlp_mc )

    use advance_xp2_xpyp_module, only: &
        update_xp2_mc

     ! Type

    implicit none

    type(grid), target, intent(in) :: gr

    !input parameters
    integer, intent(in) :: nz ! Points in the Vertical        [-]

    real( kind = core_rknd ), intent(in) :: dt ! Model timestep        [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      cloud_frac, &       !Cloud fraction                        [-]
      rcm, &              !Cloud water mixing ratio              [kg/kg]
      rvm, &              !Vapor water mixing ratio              [kg/kg]
      thlm, &             !Liquid potential temperature          [K]
      wm, &               !Mean vertical velocity                [m/s]
      exner, &            !Exner function                        [-]
      rrm_evap         !Evaporation of rain                   [kg/kg/s]
                          !It is expected that this variable is negative, as
                          !that is the convention in Morrison microphysics

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters

    !input/output variables
    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      rtp2_mc, &    !Tendency of <rt'^2> due to evaporation   [(kg/kg)^2/s]
      thlp2_mc, &   !Tendency of <thl'^2> due to evaporation  [K^2/s]
      wprtp_mc, &   !Tendency of <w'rt'> due to evaporation   [m*(kg/kg)/s^2]
      wpthlp_mc, &  !Tendency of <w'thl'> due to evaporation  [m*K/s^2] 
      rtpthlp_mc    !Tendency of <rt'thl'> due to evaporation [K*(kg/kg)/s]

    call update_xp2_mc( gr, nz, dt, cloud_frac, rcm, rvm, thlm,        & ! intent(in)
                        wm, exner, rrm_evap, pdf_params,        & ! intent(in)
                        rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc,    & ! intent(inout)
                        rtpthlp_mc ) ! intent(inout)
    return
  end subroutine update_xp2_mc_api

  !================================================================================================
  ! sat_mixrat_liq - computes the saturation mixing ratio of liquid water
  !================================================================================================
  elemental real( kind = core_rknd ) function sat_mixrat_liq_api( p_in_Pa, T_in_K )

    use saturation, only: sat_mixrat_liq

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) ::  & 
      p_in_Pa,  & ! Pressure    [Pa]
      T_in_K      ! Temperature [K]

    sat_mixrat_liq_api = sat_mixrat_liq( p_in_Pa, T_in_K )
    return
  end function sat_mixrat_liq_api

    
  !================================================================================================
  ! subroutine init_pdf_hydromet_arrays_api
  ! 
  ! DESCRIPTION: 
  !     This subroutine intializes the hydromet arrays(iirr, iiNr, etc.) to the values specified by
  !     the input arguements, this determines which hyrometeors are to be used by the microphysics
  !     scheme. It also sets up the corresponding pdf and hydromet arrays, and calculates the 
  !     subgrid variance ratio for each hydrometeor.
  ! 
  ! OPTIONAL FUNCTIONALITY:
  !     The subgrid variance ratio for each hydrometeor is calculated based on the grid spacing 
  !     defined by the host model. The calculation is a linear equation defined by a slope and
  !     intercept, each of which may or may not be passed in to this subroutine. If the slope
  !     and/or intercept are not passed in through the arguement list the default values, which 
  !     are set in the corresponding type definitions, will be used. Otherwise the values
  !     specified by the aruements will be used.
  ! 
  ! NOTES: 
  !     'hmp2_ip_on_hmm2_ip_slope_in' is of type 'hmp2_ip_on_hmm2_ip_slope_type' and
  !     'hmp2_ip_on_hmm2_ip_intrcpt_in' is of type 'hmp2_ip_on_hmm2_ip_intrcpt_in', both of which 
  !     are deinfed in corr_vrnce_module.F90, and made public through this API.
  ! 
  !     If full control over the hydrometeor variance ratios is desired, pass in slopes that are
  !     initialized to 0.0, this causes the ratios to no longer depend on the grid spacing. Then
  !     pass in the intercepts set to the values of the desired ratios.
  ! 
  ! ARGUEMENTS:
  !     host_dx (real) - Horizontal grid spacings
  !     host_dy (real)
  ! 
  !     hydromet_dim (integer) - Number of enabled hydrometeors
  ! 
  !         Each of these is an index value corresponding to a hydrometeor,
  !         used to index the hydrometeor arrays. Each index has to be unqiue
  !         for each different hyrometeor that is enabled. Setting one of these
  !         indices to -1 disables that hydrometeor
  !     iirr_in (integer) - Index of rain water mixing ratio
  !     iiri_in (integer) - Index of rain drop concentration
  !     iirs_in (integer) - Index of ice mixing ratio
  !     iirg_in (integer) - Index of ice crystal concentration
  !     iiNr_in (integer) - Index of snow mixing ratio
  !     iiNi_in (integer) - Index of snow flake concentration
  !     iiNs_in (integer) - Index of graupel mixing ratio
  !     iiNg_in (integer) - Index of graupel concentration
  ! 
  !     hmp2_ip_on_hmm2_ip_slope_in (hmp2_ip_on_hmm2_ip_slope_type) - Custom slope values
  !     hmp2_ip_on_hmm2_ip_intrcpt_in (hmp2_ip_on_hmm2_ip_intrcpt_type) - Custom intercept values
  ! 
  !================================================================================================
  subroutine init_pdf_hydromet_arrays_api( host_dx, host_dy, hydromet_dim_in,   & ! intent(in)
                                           iirr_in, iiri_in, iirs_in, iirg_in,  & ! intent(in)
                                           iiNr_in, iiNi_in, iiNs_in, iiNg_in,  & ! intent(in)
                                           hmp2_ip_on_hmm2_ip_slope_in,         & ! optional(in)
                                           hmp2_ip_on_hmm2_ip_intrcpt_in        ) ! optional(in)

    use array_index, only: &
        iirr, & ! Indicies for the hydromet arrays
        iiNr, &
        iirs, &
        iiri, &
        iirg, &
        iiNs, & 
        iiNi, &
        iiNg

    use corr_varnce_module, only: &
        init_pdf_indices,                   & ! Procedures
        init_hydromet_arrays,               &
        hmp2_ip_on_hmm2_ip_slope_type,      & ! Types
        hmp2_ip_on_hmm2_ip_intrcpt_type,    &
        hmp2_ip_on_hmm2_ip                    ! Array of hydromet ratios

    use parameters_model, only: &
        hydromet_dim

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim_in,  & ! Total number of hydrometeor species.
      iirr_in,          & ! Index of rain water mixing ratio
      iiNr_in,          & ! Index of rain drop concentration
      iiri_in,          & ! Index of ice mixing ratio
      iiNi_in,          & ! Index of ice crystal concentration
      iirs_in,          & ! Index of snow mixing ratio
      iiNs_in,          & ! Index of snow flake concentration
      iirg_in,          & ! Index of graupel mixing ratio
      iiNg_in             ! Index of graupel concentration

    real( kind = core_rknd ), intent(in) :: &
      host_dx, host_dy  ! Horizontal grid spacing, defined by host model [m]


    ! Optional Input Variables

    ! Used to overwrite default values of slope and intercept
    type( hmp2_ip_on_hmm2_ip_slope_type ), optional, intent(in) :: &
        hmp2_ip_on_hmm2_ip_slope_in     ! Custom slopes to overwrite defaults [1/m]
      
    type( hmp2_ip_on_hmm2_ip_intrcpt_type ), optional, intent(in) :: &
        hmp2_ip_on_hmm2_ip_intrcpt_in   ! Custom intercepts to overwrite defaults [-]


    ! Local Variables

    ! Slope and intercept are initialized with default values
    type( hmp2_ip_on_hmm2_ip_slope_type ) :: &
        hmp2_ip_on_hmm2_ip_slope        ! Slopes used to calculated hydromet variance [1/m]
      
    type( hmp2_ip_on_hmm2_ip_intrcpt_type ) :: &
        hmp2_ip_on_hmm2_ip_intrcpt      ! Intercepts used to calculated hydromet variance [-]

    !----------------------- Begin Code -----------------------------

    ! If slope and intercept are present in call, then overwrite default values
    if ( present( hmp2_ip_on_hmm2_ip_slope_in ) ) then
        hmp2_ip_on_hmm2_ip_slope = hmp2_ip_on_hmm2_ip_slope_in
    end if

    if ( present( hmp2_ip_on_hmm2_ip_intrcpt_in ) ) then
        hmp2_ip_on_hmm2_ip_intrcpt = hmp2_ip_on_hmm2_ip_intrcpt_in
    end if


    ! Initialize the hydromet indices and hydromet_dim
    hydromet_dim = hydromet_dim_in
    iirr = iirr_in
    iiri = iiri_in
    iirs = iirs_in
    iirg = iirg_in
    iiNr = iiNr_in
    iiNi = iiNi_in
    iiNs = iiNs_in
    iiNg = iiNg_in

    ! Calculate the subgrid variances of the hydrometeors
    allocate( hmp2_ip_on_hmm2_ip(hydromet_dim) )

    if ( iirr > 0 ) then
       hmp2_ip_on_hmm2_ip(iirr) = hmp2_ip_on_hmm2_ip_intrcpt%rr + &
                                  hmp2_ip_on_hmm2_ip_slope%rr * max( host_dx, host_dy )
    endif

    if ( iirs > 0 ) then
       hmp2_ip_on_hmm2_ip(iirs) = hmp2_ip_on_hmm2_ip_intrcpt%rs + &
                                  hmp2_ip_on_hmm2_ip_slope%rs * max( host_dx, host_dy )
    endif

    if ( iiri > 0 ) then
       hmp2_ip_on_hmm2_ip(iiri) = hmp2_ip_on_hmm2_ip_intrcpt%ri + &
                                  hmp2_ip_on_hmm2_ip_slope%ri * max( host_dx, host_dy )
    endif

    if ( iirg > 0 ) then
       hmp2_ip_on_hmm2_ip(iirg) = hmp2_ip_on_hmm2_ip_intrcpt%rg + &
                                  hmp2_ip_on_hmm2_ip_slope%rg * max( host_dx, host_dy )
    endif

    if ( iiNr > 0 ) then
       hmp2_ip_on_hmm2_ip(iiNr) = hmp2_ip_on_hmm2_ip_intrcpt%Nr + &
                                  hmp2_ip_on_hmm2_ip_slope%Nr * max( host_dx, host_dy )
    endif

    if ( iiNs > 0 ) then
       hmp2_ip_on_hmm2_ip(iiNs) = hmp2_ip_on_hmm2_ip_intrcpt%Ns + &
                                  hmp2_ip_on_hmm2_ip_slope%Ns * max( host_dx, host_dy )
    endif

    if ( iiNi > 0 ) then
       hmp2_ip_on_hmm2_ip(iiNi) = hmp2_ip_on_hmm2_ip_intrcpt%Ni + &
                                  hmp2_ip_on_hmm2_ip_slope%Ni * max( host_dx, host_dy )
    endif

    if ( iiNg > 0 ) then
       hmp2_ip_on_hmm2_ip(iiNg) = hmp2_ip_on_hmm2_ip_intrcpt%Ng + &
                                  hmp2_ip_on_hmm2_ip_slope%Ng * max( host_dx, host_dy )
    endif

    ! Hydromet arrays are Initialized based on the hydromet indices
    call init_hydromet_arrays( hydromet_dim, iirr, iiNr,    & ! intent(in)
                               iiri, iiNi, iirs, iiNs,      & ! intent(in)
                               iirg, iiNg )                   ! intent(in)

    
    ! Initialize the PDF indices based on the hydromet indices
    call init_pdf_indices( hydromet_dim,iirr, iiNr, & ! intent(in)
                           iiri, iiNi, iirs, iiNs,  & ! intent(in)
                           iirg, iiNg )               ! intent(in)

    return

  end subroutine init_pdf_hydromet_arrays_api

  !================================================================================================
  ! set_default_parameters: Sets all CLUBB tunable parameters to a default setting
  !================================================================================================
  subroutine set_default_parameters_api( &
               C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
               C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
               C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
               C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
               C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
               C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
               c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
               slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
               coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
               gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
               omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
               lambda0_stability_coef, mult_coef, taumin, taumax, &
               Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
               Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
               thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
               Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
               altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
               C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, &
               C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
               C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
               C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
               Cx_min, Cx_max, Richardson_num_min, &
               Richardson_num_max, a3_coef_min )

    use parameters_tunable, only: &
        set_default_parameters    ! Procedure(s)

    implicit none

    ! Output variables
    real( kind = core_rknd ), intent(out) :: &
      C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
      C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
      C7, C7b, C7c, C8, C8b, C10, &
      C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
      C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
      C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
      c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8,  &
      c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
      slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
      coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
      gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
      omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
      lambda0_stability_coef, mult_coef, taumin, taumax, Lscale_mu_coef, &
      Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, c_K10h, &
      thlp2_rad_coef, thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
      Skw_max_mag, xp3_coef_base, xp3_coef_slope, altitude_threshold, &
      rtp2_clip_coef, C_invrs_tau_bkgnd, C_invrs_tau_sfc, &
      C_invrs_tau_shear, C_invrs_tau_N2, C_invrs_tau_N2_wp2, &
      C_invrs_tau_N2_xp2, C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
      C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, a3_coef_min

    call set_default_parameters( &
               C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
               C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, &
               C6thl, C6thlb, C6thlc, C7, C7b, C7c, C8, C8b, C10, &
               C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
               C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, &
               C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, &
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, &
               c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
               slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
               coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
               gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, &
               omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
               lambda0_stability_coef, mult_coef, taumin, taumax, &
               Lscale_mu_coef, Lscale_pert_coef, alpha_corr, &
               Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, &
               thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
               Skw_max_mag, xp3_coef_base, xp3_coef_slope, &
               altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, &
               C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, &
               C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &
               C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
               C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
               Cx_min, Cx_max, Richardson_num_min, &
               Richardson_num_max, a3_coef_min )

    return

  end subroutine set_default_parameters_api

  !================================================================================================
  ! set_default_clubb_config_flags: Sets all CLUBB flags to a default setting
  !================================================================================================
  subroutine set_default_clubb_config_flags_api( iiPDF_type, & ! Out
                                                 ipdf_call_placement, & ! Out
                                                 l_use_precip_frac, & ! Out
                                                 l_predict_upwp_vpwp, & ! Out
                                                 l_min_wp2_from_corr_wx, & ! Out
                                                 l_min_xp2_from_corr_wx, & ! Out
                                                 l_C2_cloud_frac, & ! Out
                                                 l_diffuse_rtm_and_thlm, & ! Out
                                                 l_stability_correct_Kh_N2_zm, & ! Out
                                                 l_calc_thlp2_rad, & ! Out
                                                 l_upwind_xpyp_ta, & ! Out
                                                 l_upwind_xm_ma, & ! Out
                                                 l_uv_nudge, & ! Out
                                                 l_rtm_nudge, & ! Out
                                                 l_tke_aniso, & ! Out
                                                 l_vert_avg_closure, & ! Out
                                                 l_trapezoidal_rule_zt, & ! Out
                                                 l_trapezoidal_rule_zm, & ! Out
                                                 l_call_pdf_closure_twice, & ! Out
                                                 l_standard_term_ta, & ! Out
                                                 l_partial_upwind_wp3, & ! Out
                                                 l_godunov_upwind_wpxp_ta, & ! Out
                                                 l_godunov_upwind_xpyp_ta, & ! Out
                                                 l_use_cloud_cover, & ! Out
                                                 l_diagnose_correlations, & ! Out
                                                 l_calc_w_corr, & ! Out
                                                 l_const_Nc_in_cloud, & ! Out
                                                 l_fix_w_chi_eta_correlations, & ! Out
                                                 l_stability_correct_tau_zm, & ! Out
                                                 l_damp_wp2_using_em, & ! Out
                                                 l_do_expldiff_rtm_thlm, & ! Out
                                                 l_Lscale_plume_centered, & ! Out
                                                 l_diag_Lscale_from_tau, & ! Out
                                                 l_use_C7_Richardson, & ! Out
                                                 l_use_C11_Richardson, & ! Out
                                                 l_use_shear_Richardson, & ! Out
                                                 l_brunt_vaisala_freq_moist, & ! Out
                                                 l_use_thvm_in_bv_freq, & ! Out
                                                 l_rcm_supersat_adj, & ! Out
                                                 l_damp_wp3_Skw_squared, & ! Out
                                                 l_prescribed_avg_deltaz, & ! Out
                                                 l_lmm_stepping, & ! Out
                                                 l_e3sm_config, & ! Out
                                                 l_vary_convect_depth, & ! Out
                                                 l_use_tke_in_wp3_pr_turb_term, & ! Out
                                                 l_use_tke_in_wp2_wp3_K_dfsn, & ! Out
                                                 l_smooth_Heaviside_tau_wpxp, & ! Out
                                                 l_enable_relaxed_clipping, & ! Out
                                                 l_linearize_pbl_winds ) ! Out

    use model_flags, only: &
        set_default_clubb_config_flags  ! Procedure

    implicit none

    ! Output variables
    integer, intent(out) :: &
      iiPDF_type,          & ! Selected option for the two-component normal
                             ! (double Gaussian) PDF type to use for the w, rt,
                             ! and theta-l (or w, chi, and eta) portion of
                             ! CLUBB's multivariate, two-component PDF.
      ipdf_call_placement    ! Selected option for the placement of the call to
                             ! CLUBB's PDF.

    logical, intent(out) :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_min_wp2_from_corr_wx,       & ! Flag to base the threshold minimum value of wp2 on keeping
                                      ! the overall correlation of w and x (w and rt, as well as w
                                      ! and theta-l) within the limits of -max_mag_correlation_flux
                                      ! to max_mag_correlation_flux.
      l_min_xp2_from_corr_wx,       & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                      ! thlp2) on keeping the overall correlation of w and x within
                                      ! the limits of -max_mag_correlation_flux to
                                      ! max_mag_correlation_flux.
      l_C2_cloud_frac,              & ! Flag to use cloud fraction to adjust the value of the
                                      ! turbulent dissipation coefficient, C2.
      l_diffuse_rtm_and_thlm,       & ! Diffuses rtm and thlm
      l_stability_correct_Kh_N2_zm, & ! Divides Kh_N2_zm by a stability factor
      l_calc_thlp2_rad,             & ! Include the contribution of radiation to thlp2
      l_upwind_xpyp_ta,             & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, &
                                      ! sclrpthlp.
      l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtm, thlm, sclrm, um and vm.
      l_uv_nudge,                   & ! For wind speed nudging.
      l_rtm_nudge,                  & ! For rtm nudging
      l_tke_aniso,                  & ! For anisotropic turbulent kinetic energy, i.e.
                                      ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_vert_avg_closure,           & ! Use 2 calls to pdf_closure and the trapezoidal rule to
                                      ! compute the varibles that are output from high order
                                      ! closure
      l_trapezoidal_rule_zt,        & ! If true, the trapezoidal rule is called for the
                                      ! thermodynamic-level variables output from pdf_closure.
      l_trapezoidal_rule_zm,        & ! If true, the trapezoidal rule is called for three
                                      ! momentum-level variables - wpthvp, thlpthvp, and rtpthvp -
                                      ! output from pdf_closure.
      l_call_pdf_closure_twice,     & ! This logical flag determines whether or not to call
                                      ! subroutine pdf_closure twice.  If true, pdf_closure is
                                      ! called first on thermodynamic levels and then on momentum
                                      ! levels so that each variable is computed on its native
                                      ! level.  If false, pdf_closure is only called on
                                      ! thermodynamic levels, and variables which belong on
                                      ! momentum levels are interpolated.
      l_standard_term_ta,           & ! Use the standard discretization for the turbulent advection
                                      ! terms.  Setting to .false. means that a_1 and a_3 are
                                      ! pulled outside of the derivative in
                                      ! advance_wp2_wp3_module.F90 and in
                                      ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3,         & ! Flag to use an "upwind" discretization rather
                                      ! than a centered discretization for the portion
                                      ! of the wp3 turbulent advection term for ADG1
                                      ! that is linearized in terms of wp3<t+1>.
                                      ! (Requires ADG1 PDF and l_standard_term_ta).
      l_godunov_upwind_wpxp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. 
                                      ! It affects  wpxp only.
      l_godunov_upwind_xpyp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. It affects
                                      ! xpyp only.
      l_use_cloud_cover,            & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac
                                      ! and rcm to help increase cloudiness at coarser grid
                                      ! resolutions.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations, & ! Use a fixed correlation for s and t Mellor(chi/eta)
      l_stability_correct_tau_zm,   & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability
                                      ! correction
      l_damp_wp2_using_em,          & ! In wp2 equation, use a dissipation formula of
                                      ! -(2/3)*em/tau_zm, as in Bougeault (1981)
      l_do_expldiff_rtm_thlm,       & ! Diffuse rtm and thlm explicitly
      l_Lscale_plume_centered,      & ! Alternate that uses the PDF to compute the perturbed values
      l_diag_Lscale_from_tau,       & ! First diagnose dissipation time tau, and then diagnose the
                                      ! mixing length scale as Lscale = tau * tke
      l_use_C7_Richardson,          & ! Parameterize C7 based on Richardson number
      l_use_C11_Richardson,         & ! Parameterize C11 and C16 based on Richardson number
      l_use_shear_Richardson,       & ! Use shear in the calculation of Richardson number
      l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                      ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq,        & ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_rcm_supersat_adj,           & ! Add excess supersaturated vapor to cloud water
      l_damp_wp3_Skw_squared,       & ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_prescribed_avg_deltaz,      & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
      l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
      l_e3sm_config,                & ! Run model with E3SM settings
      l_vary_convect_depth,         & ! Flag used to calculate convective velocity using
                                      ! a variable estimate of layer depth based on the depth
                                      ! over which wpthlp is positive near the ground when true
                                      ! More information can be found by
                                      ! Looking at issue #905 on the clubb repo
      l_use_tke_in_wp3_pr_turb_term,& ! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn,  & ! Use TKE in eddy diffusion for wp2 and wp3
      l_smooth_Heaviside_tau_wpxp,  & ! Use smoothed Heaviside 'Peskin' function
                                      ! in the calculation of H_invrs_tau_wpxp_N2
                                      ! in src/CLUBB_core/mixing_length.F90
      l_enable_relaxed_clipping,    & ! Flag to relax clipping on wpxp
                                      ! in xm_wpxp_clipping_and_stats
      l_linearize_pbl_winds           ! Code to linearize PBL winds

    call set_default_clubb_config_flags( iiPDF_type, & ! Out
                                         ipdf_call_placement, & ! Out
                                         l_use_precip_frac, & ! Out
                                         l_predict_upwp_vpwp, & ! Out
                                         l_min_wp2_from_corr_wx, & ! Out
                                         l_min_xp2_from_corr_wx, & ! Out
                                         l_C2_cloud_frac, & ! Out
                                         l_diffuse_rtm_and_thlm, & ! Out
                                         l_stability_correct_Kh_N2_zm, & ! Out
                                         l_calc_thlp2_rad, & ! Out
                                         l_upwind_xpyp_ta, & ! Out
                                         l_upwind_xm_ma, & ! Out
                                         l_uv_nudge, & ! Out
                                         l_rtm_nudge, & ! Out
                                         l_tke_aniso, & ! Out
                                         l_vert_avg_closure, & ! Out
                                         l_trapezoidal_rule_zt, & ! Out
                                         l_trapezoidal_rule_zm, & ! Out
                                         l_call_pdf_closure_twice, & ! Out
                                         l_standard_term_ta, & ! Out
                                         l_partial_upwind_wp3, & ! Out
                                         l_godunov_upwind_wpxp_ta, & ! Out
                                         l_godunov_upwind_xpyp_ta, & ! Out
                                         l_use_cloud_cover, & ! Out
                                         l_diagnose_correlations, & ! Out
                                         l_calc_w_corr, & ! Out
                                         l_const_Nc_in_cloud, & ! Out
                                         l_fix_w_chi_eta_correlations, & ! Out
                                         l_stability_correct_tau_zm, & ! Out
                                         l_damp_wp2_using_em, & ! Out
                                         l_do_expldiff_rtm_thlm, & ! Out
                                         l_Lscale_plume_centered, & ! Out
                                         l_diag_Lscale_from_tau, & ! Out
                                         l_use_C7_Richardson, & ! Out
                                         l_use_C11_Richardson, & ! Out
                                         l_use_shear_Richardson, & ! Out
                                         l_brunt_vaisala_freq_moist, & ! Out
                                         l_use_thvm_in_bv_freq, & ! Out
                                         l_rcm_supersat_adj, & ! Out
                                         l_damp_wp3_Skw_squared, & ! Out
                                         l_prescribed_avg_deltaz, & ! Out
                                         l_lmm_stepping, & ! Out
                                         l_e3sm_config, & ! Out
                                         l_vary_convect_depth, & ! Out
                                         l_use_tke_in_wp3_pr_turb_term, & ! Out
                                         l_use_tke_in_wp2_wp3_K_dfsn, & ! Out
                                         l_smooth_Heaviside_tau_wpxp, & ! Out
                                         l_enable_relaxed_clipping, & ! Out
                                         l_linearize_pbl_winds ) ! Out

  end subroutine set_default_clubb_config_flags_api

  !================================================================================================
  ! initialize_clubb_config_flags_type: Initialize the clubb_config_flags_type
  !================================================================================================
  subroutine initialize_clubb_config_flags_type_api( iiPDF_type, & ! In
                                                     ipdf_call_placement, & ! In
                                                     l_use_precip_frac, & ! In
                                                     l_predict_upwp_vpwp, & ! In
                                                     l_min_wp2_from_corr_wx, & ! In
                                                     l_min_xp2_from_corr_wx, & ! In
                                                     l_C2_cloud_frac, & ! In
                                                     l_diffuse_rtm_and_thlm, & ! In
                                                     l_stability_correct_Kh_N2_zm, & ! In
                                                     l_calc_thlp2_rad, & ! In
                                                     l_upwind_xpyp_ta, & ! In
                                                     l_upwind_xm_ma, & ! In
                                                     l_uv_nudge, & ! In
                                                     l_rtm_nudge, & ! In
                                                     l_tke_aniso, & ! In
                                                     l_vert_avg_closure, & ! In
                                                     l_trapezoidal_rule_zt, & ! In
                                                     l_trapezoidal_rule_zm, & ! In
                                                     l_call_pdf_closure_twice, & ! In
                                                     l_standard_term_ta, & ! In
                                                     l_partial_upwind_wp3, & ! In
                                                     l_godunov_upwind_wpxp_ta, & ! In
                                                     l_godunov_upwind_xpyp_ta, & ! In
                                                     l_use_cloud_cover, & ! In
                                                     l_diagnose_correlations, & ! In
                                                     l_calc_w_corr, & ! In
                                                     l_const_Nc_in_cloud, & ! In
                                                     l_fix_w_chi_eta_correlations, & ! In
                                                     l_stability_correct_tau_zm, & ! In
                                                     l_damp_wp2_using_em, & ! In
                                                     l_do_expldiff_rtm_thlm, & ! In
                                                     l_Lscale_plume_centered, & ! In
                                                     l_diag_Lscale_from_tau, & ! In
                                                     l_use_C7_Richardson, & ! In
                                                     l_use_C11_Richardson, & ! In
                                                     l_use_shear_Richardson, & ! In
                                                     l_brunt_vaisala_freq_moist, & ! In
                                                     l_use_thvm_in_bv_freq, & ! In
                                                     l_rcm_supersat_adj, & ! In
                                                     l_damp_wp3_Skw_squared, & ! In
                                                     l_prescribed_avg_deltaz, & ! In
                                                     l_lmm_stepping, & ! In
                                                     l_e3sm_config, & ! In
                                                     l_vary_convect_depth, & ! In
                                                     l_use_tke_in_wp3_pr_turb_term, & ! In
                                                     l_use_tke_in_wp2_wp3_K_dfsn, & ! In
                                                     l_smooth_Heaviside_tau_wpxp, & ! In
                                                     l_enable_relaxed_clipping, & ! In
                                                     l_linearize_pbl_winds, & ! In
                                                     clubb_config_flags ) ! Out

    use model_flags, only: &
        clubb_config_flags_type, &          ! Type
        initialize_clubb_config_flags_type  ! Procedure

    implicit none

    ! Input variables
    integer, intent(in) :: &
      iiPDF_type,          & ! Selected option for the two-component normal
                             ! (double Gaussian) PDF type to use for the w, rt,
                             ! and theta-l (or w, chi, and eta) portion of
                             ! CLUBB's multivariate, two-component PDF.
      ipdf_call_placement    ! Selected option for the placement of the call to
                             ! CLUBB's PDF.

    logical, intent(in) :: &
      l_use_precip_frac,            & ! Flag to use precipitation fraction in KK microphysics. The
                                      ! precipitation fraction is automatically set to 1 when this
                                      ! flag is turned off.
      l_predict_upwp_vpwp,          & ! Flag to predict <u'w'> and <v'w'> along with <u> and <v>
                                      ! alongside the advancement of <rt>, <w'rt'>, <thl>,
                                      ! <wpthlp>, <sclr>, and <w'sclr'> in subroutine
                                      ! advance_xm_wpxp.  Otherwise, <u'w'> and <v'w'> are still
                                      ! approximated by eddy diffusivity when <u> and <v> are
                                      ! advanced in subroutine advance_windm_edsclrm.
      l_min_wp2_from_corr_wx,       & ! Flag to base the threshold minimum value of wp2 on keeping
                                      ! the overall correlation of w and x (w and rt, as well as w
                                      ! and theta-l) within the limits of -max_mag_correlation_flux
                                      ! to max_mag_correlation_flux.
      l_min_xp2_from_corr_wx,       & ! Flag to base the threshold minimum value of xp2 (rtp2 and
                                      ! thlp2) on keeping the overall correlation of w and x within
                                      ! the limits of -max_mag_correlation_flux to
                                      ! max_mag_correlation_flux.
      l_C2_cloud_frac,              & ! Flag to use cloud fraction to adjust the value of the
                                      ! turbulent dissipation coefficient, C2.
      l_diffuse_rtm_and_thlm,       & ! Diffuses rtm and thlm
      l_stability_correct_Kh_N2_zm, & ! Divides Kh_N2_zm by a stability factor
      l_calc_thlp2_rad,             & ! Include the contribution of radiation to thlp2
      l_upwind_xpyp_ta,             & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtp2, thlp2, up2, vp2, sclrp2, rtpthlp, sclrprtp, &
                                      ! sclrpthlp.
      l_upwind_xm_ma,               & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent or mean advection terms. It
                                      ! affects rtm, thlm, sclrm, um and vm.
      l_uv_nudge,                   & ! For wind speed nudging.
      l_rtm_nudge,                  & ! For rtm nudging
      l_tke_aniso,                  & ! For anisotropic turbulent kinetic energy, i.e.
                                      ! TKE = 1/2 (u'^2 + v'^2 + w'^2)
      l_vert_avg_closure,           & ! Use 2 calls to pdf_closure and the trapezoidal rule to
                                      ! compute the varibles that are output from high order
                                      ! closure
      l_trapezoidal_rule_zt,        & ! If true, the trapezoidal rule is called for the
                                      ! thermodynamic-level variables output from pdf_closure.
      l_trapezoidal_rule_zm,        & ! If true, the trapezoidal rule is called for three
                                      ! momentum-level variables - wpthvp, thlpthvp, and rtpthvp -
                                      ! output from pdf_closure.
      l_call_pdf_closure_twice,     & ! This logical flag determines whether or not to call
                                      ! subroutine pdf_closure twice.  If true, pdf_closure is
                                      ! called first on thermodynamic levels and then on momentum
                                      ! levels so that each variable is computed on its native
                                      ! level.  If false, pdf_closure is only called on
                                      ! thermodynamic levels, and variables which belong on
                                      ! momentum levels are interpolated.
      l_standard_term_ta,           & ! Use the standard discretization for the turbulent advection
                                      ! terms.  Setting to .false. means that a_1 and a_3 are
                                      ! pulled outside of the derivative in
                                      ! advance_wp2_wp3_module.F90 and in
                                      ! advance_xp2_xpyp_module.F90.
      l_partial_upwind_wp3,         & ! Flag to use an "upwind" discretization rather
                                      ! than a centered discretization for the portion
                                      ! of the wp3 turbulent advection term for ADG1
                                      ! that is linearized in terms of wp3<t+1>.
                                      ! (Requires ADG1 PDF and l_standard_term_ta).
      l_godunov_upwind_wpxp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered
                                      ! differencing for turbulent advection terms. 
                                      ! It affects  wpxp only.
      l_godunov_upwind_xpyp_ta,     & ! This flag determines whether we want to use an upwind
                                      ! differencing approximation rather than a centered 
                                      ! differencing for turbulent advection terms. It affects
                                      ! xpyp only.
      l_use_cloud_cover,            & ! Use cloud_cover and rcm_in_layer to help boost cloud_frac
                                      ! and rcm to help increase cloudiness at coarser grid
                                      ! resolutions.
      l_diagnose_correlations,      & ! Diagnose correlations instead of using fixed ones
      l_calc_w_corr,                & ! Calculate the correlations between w and the hydrometeors
      l_const_Nc_in_cloud,          & ! Use a constant cloud droplet conc. within cloud (K&K)
      l_fix_w_chi_eta_correlations, & ! Use a fixed correlation for s and t Mellor(chi/eta)
      l_stability_correct_tau_zm,   & ! Use tau_N2_zm instead of tau_zm in wpxp_pr1 stability
                                      ! correction
      l_damp_wp2_using_em,          & ! In wp2 equation, use a dissipation formula of
                                      ! -(2/3)*em/tau_zm, as in Bougeault (1981)
      l_do_expldiff_rtm_thlm,       & ! Diffuse rtm and thlm explicitly
      l_Lscale_plume_centered,      & ! Alternate that uses the PDF to compute the perturbed values
      l_diag_Lscale_from_tau,       & ! First diagnose dissipation time tau, and then diagnose the
                                      ! mixing length scale as Lscale = tau * tke
      l_use_C7_Richardson,          & ! Parameterize C7 based on Richardson number
      l_use_C11_Richardson,         & ! Parameterize C11 and C16 based on Richardson number
      l_use_shear_Richardson,       & ! Use shear in the calculation of Richardson number
      l_brunt_vaisala_freq_moist,   & ! Use a different formula for the Brunt-Vaisala frequency in
                                      ! saturated atmospheres (from Durran and Klemp, 1982)
      l_use_thvm_in_bv_freq,        & ! Use thvm in the calculation of Brunt-Vaisala frequency
      l_rcm_supersat_adj,           & ! Add excess supersaturated vapor to cloud water
      l_damp_wp3_Skw_squared,       & ! Set damping on wp3 to use Skw^2 rather than Skw^4
      l_prescribed_avg_deltaz,      & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
      l_lmm_stepping,               & ! Apply Linear Multistep Method (LMM) Stepping
      l_e3sm_config,                & ! Run model with E3SM settings
      l_vary_convect_depth,         & ! Flag used to calculate convective velocity using
                                      ! a variable estimate of layer depth based on the depth
                                      ! over which wpthlp is positive near the ground when true
                                      ! More information can be found by
                                      ! Looking at issue #905 on the clubb repo
      l_use_tke_in_wp3_pr_turb_term,& ! Use TKE formulation for wp3 pr_turb term
      l_use_tke_in_wp2_wp3_K_dfsn,  & ! Use TKE in eddy diffusion for wp2 and wp3
      l_smooth_Heaviside_tau_wpxp,  & ! Use smoothed Heaviside 'Peskin' function
                                      ! in the calculation of H_invrs_tau_wpxp_N2
                                      ! in src/CLUBB_core/mixing_length.F90
      l_enable_relaxed_clipping,    & ! Flag to relax clipping on wpxp
                                      ! in xm_wpxp_clipping_and_stats
      l_linearize_pbl_winds           ! Code to linearize PBL winds

    ! Output variables
    type(clubb_config_flags_type), intent(out) :: &
      clubb_config_flags            ! Derived type holding all configurable CLUBB flags

    call initialize_clubb_config_flags_type( iiPDF_type, & ! In
                                             ipdf_call_placement, & ! In 
                                             l_use_precip_frac, & ! In
                                             l_predict_upwp_vpwp, & ! In
                                             l_min_wp2_from_corr_wx, & ! In
                                             l_min_xp2_from_corr_wx, & ! In
                                             l_C2_cloud_frac, & ! In
                                             l_diffuse_rtm_and_thlm, & ! In
                                             l_stability_correct_Kh_N2_zm, & ! In
                                             l_calc_thlp2_rad, & ! In
                                             l_upwind_xpyp_ta, & ! In
                                             l_upwind_xm_ma, & ! In
                                             l_uv_nudge, & ! In
                                             l_rtm_nudge, & ! In
                                             l_tke_aniso, & ! In
                                             l_vert_avg_closure, & ! In
                                             l_trapezoidal_rule_zt, & ! In
                                             l_trapezoidal_rule_zm, & ! In
                                             l_call_pdf_closure_twice, & ! In
                                             l_standard_term_ta, & ! In
                                             l_partial_upwind_wp3, & ! In
                                             l_godunov_upwind_wpxp_ta, & ! In
                                             l_godunov_upwind_xpyp_ta, & ! In
                                             l_use_cloud_cover, & ! In
                                             l_diagnose_correlations, & ! In
                                             l_calc_w_corr, & ! In
                                             l_const_Nc_in_cloud, & ! In
                                             l_fix_w_chi_eta_correlations, & ! In
                                             l_stability_correct_tau_zm, & ! In
                                             l_damp_wp2_using_em, & ! In
                                             l_do_expldiff_rtm_thlm, & ! In
                                             l_Lscale_plume_centered, & ! In
                                             l_diag_Lscale_from_tau, & ! In
                                             l_use_C7_Richardson, & ! In
                                             l_use_C11_Richardson, & ! In
                                             l_use_shear_Richardson, & ! In
                                             l_brunt_vaisala_freq_moist, & ! In
                                             l_use_thvm_in_bv_freq, & ! In
                                             l_rcm_supersat_adj, & ! In
                                             l_damp_wp3_Skw_squared, & ! In
                                             l_prescribed_avg_deltaz, & ! In
                                             l_lmm_stepping, & ! In
                                             l_e3sm_config, & ! In
                                             l_vary_convect_depth, & ! In
                                             l_use_tke_in_wp3_pr_turb_term, & ! In
                                             l_use_tke_in_wp2_wp3_K_dfsn, & ! In
                                             l_smooth_Heaviside_tau_wpxp, & ! In
                                             l_enable_relaxed_clipping, & ! In
                                             l_linearize_pbl_winds, & ! In
                                             clubb_config_flags ) ! Out

  end subroutine initialize_clubb_config_flags_type_api

  !================================================================================================
  ! print_clubb_config_flags: Prints the clubb_config_flags
  !================================================================================================
  subroutine print_clubb_config_flags_api( iunit, clubb_config_flags ) ! In

    use model_flags, only: &
        clubb_config_flags_type, &          ! Type
        print_clubb_config_flags            ! Procedure

    implicit none

    ! Input variables
    integer, intent(in) :: &
      iunit ! The file to write to

    type(clubb_config_flags_type), intent(in) :: &
      clubb_config_flags ! Derived type holding all configurable CLUBB flags

    call print_clubb_config_flags( iunit, clubb_config_flags ) ! In

  end subroutine print_clubb_config_flags_api

  !================================================================================================
  ! initialize_tau_sponge_damp
  !================================================================================================
  subroutine initialize_tau_sponge_damp_api( gr, dt, z, settings, damping_profile )

    use sponge_layer_damping, only: &
        sponge_damp_settings,       & ! Variable(s)
        sponge_damp_profile,        &
        initialize_tau_sponge_damp    ! Procedure(s)

    

    implicit none

    type(grid), target, intent(in) :: gr

    ! Input Variable(s)
    real( kind = core_rknd ), intent(in) :: &
      dt    ! Model Timestep    [s]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      z    ! Height of model grid levels    [m]

    type(sponge_damp_settings), intent(in) :: &
      settings

    ! Output Variable(s)
    type(sponge_damp_profile), intent(out) :: &
      damping_profile

    call initialize_tau_sponge_damp( gr, dt, z, settings, & ! intent(in)
                                     damping_profile ) ! intent(inout)

  end subroutine initialize_tau_sponge_damp_api

  !================================================================================================
  ! finalize_tau_sponge_damp
  !================================================================================================
  subroutine finalize_tau_sponge_damp_api( damping_profile )

    use sponge_layer_damping, only: &
        sponge_damp_profile,      & ! Variable(s)
        finalize_tau_sponge_damp    ! Procedure(s)

    implicit none

    ! Input/Output Variable(s)
    type(sponge_damp_profile), intent(inout) :: &
      damping_profile ! Information for damping the profile

    call finalize_tau_sponge_damp( damping_profile ) ! intent(inout)

  end subroutine finalize_tau_sponge_damp_api
    
end module clubb_api_module
