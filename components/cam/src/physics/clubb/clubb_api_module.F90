!--------------------------------------------------------------------------------------------------
! $Id: clubb_api_module.F90 8220 2016-07-21 18:48:32Z raut@uwm.edu $
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
    iiNgm, & ! Hydrometeor array index for graupel concentration, Ng
    iiNim, & ! Hydrometeor array index for ice concentration, Ni
    iiNrm, & ! Hydrometeor array index for rain drop concentration, Nr
    iiNsm, & ! Hydrometeor array index for snow concentration, Ns
    iirgm, & ! Hydrometeor array index for graupel mixing ratio, rg
    iirim, & ! Hydrometeor array index for ice mixing ratio, ri
    iirrm, & ! Hydrometeor array index for rain water mixing ratio, rr
    iirsm, & ! Hydrometeor array index for snow mixing ratio, rs
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
      d_variables,        &
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
      hmp2_ip_on_hmm2_ip_ratios_type

  use error_code, only : &
    clubb_no_error ! Enum representing that no errors have occurred in CLUBB

  use grid_class, only : &
    gr

  use hydromet_pdf_parameter_module, only : &
    hydromet_pdf_parameter

  use model_flags, only : &
    l_use_boussinesq, & ! Use Boussinesq form of predictive equations (default is Anelastic).
    l_diagnose_correlations, & ! Diagnose correlations instead of using fixed ones
    l_calc_w_corr, & ! Calculate the correlations between w and the hydrometeors
    l_use_cloud_cover, & ! helps to increase cloudiness at coarser grid resolutions.
    l_use_precip_frac, & ! Flag to use precipitation fraction in KK microphysics.
    l_tke_aniso, & ! For anisotropic turbulent kinetic energy
    l_fix_chi_eta_correlations, & ! Use a fixed correlation for s and t Mellor(chi/eta)
    l_const_Nc_in_cloud, &   ! Use a constant cloud droplet conc. within cloud (K&K)
    l_diffuse_rtm_and_thlm, &
    l_stability_correct_Kh_N2_zm, &
    l_stability_correct_tau_zm, &
    l_do_expldiff_rtm_thlm, &
    l_Lscale_plume_centered, &
    l_use_ice_latent, &
    l_use_C7_Richardson, &
    l_use_C11_Richardson, &
    l_brunt_vaisala_freq_moist, &
    l_use_thvm_in_bv_freq, &
    l_rcm_supersat_adj

  use parameters_model, only : &
    hydromet_dim    ! Number of hydrometeor species

  use parameters_tunable, only : &
    l_prescribed_avg_deltaz, & ! used in adj_low_res_nu. If .true., avg_deltaz = deltaz
    mu

  use parameter_indices, only:  &
    nparams, & ! Variable(s)
    iC1, iC1b, iC1c, iC2, iC2b, iC2c,  &
    iC2rt, iC2thl, iC2rtthl, iC4, iC5, &
    iC6rt, iC6rtb, iC6rtc, iC6thl, iC6thlb, iC6thlc, &
    iC7, iC7b, iC7c, iC8, iC8b, iC10, iC11, iC11b, iC11c, &
    iC12, iC13, iC14, iC15, iC6rt_Lscale0, iC6thl_Lscale0, &
    iC7_Lscale0, iwpxp_L_thresh, ic_K, ic_K1, inu1, ic_K2, inu2, &
    ic_K6, inu6, ic_K8, inu8, ic_K9, inu9, inu10, ic_K_hm, ic_K_hmb, iK_hm_min_coef, &
    inu_hm, ibeta, igamma_coef, igamma_coefb, igamma_coefc, ilmin_coef, &
    iomicron, izeta_vrnce_rat, iupsilon_precip_frac_rat, &
    ilambda0_stability_coef, imult_coef, itaumin, itaumax, imu, iLscale_mu_coef, &
    iLscale_pert_coef, ialpha_corr, iSkw_denom_coef, ic_K10, ic_K10h, ithlp2_rad_coef, &
    ithlp2_rad_cloud_frac_thresh, iup2_vp2_factor

  use pdf_parameter_module, only : &
#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
    num_pdf_params, &
#endif
    pdf_parameter

  use stat_file_module, only : &
    clubb_i, &    ! Used to output multiple columns
    clubb_j       ! The indices must not exceed nlon (for i) or nlat (for j).

  use stats_rad_zm_module, only : &
    nvarmax_rad_zm ! Maximum variables allowed

  use stats_rad_zt_module, only : &
    nvarmax_rad_zt  ! Maximum variables allowed

  use stats_variables, only : &
    stats_zt, & ! zt grid
    stats_zm, & ! zm grid
    stats_rad_zt, & ! rad_zt grid
    stats_rad_zm, & ! rad_zm grid
    stats_sfc, &
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

  use variables_diagnostic_module, only : &
    Lscale, & ! Mixing lengths
    wp2_zt, & ! w'^2 on thermo. grid     [m^2/s^2]
    wphydrometp ! Covariance of w and hydrometeor (momentum levels) [(m/s)un]

  implicit none

  private

  public &
    ! To Implement CLUBB:
    read_parameters_api, &
    setup_clubb_core_api, &
        ! CLUBB can be set more specifically using these flags:
        l_use_boussinesq, &
        l_diagnose_correlations, &
        l_calc_w_corr, &
        l_use_cloud_cover, &
        l_use_precip_frac, &
        l_tke_aniso, &
        l_fix_chi_eta_correlations, &
        l_const_Nc_in_cloud, &
        l_diffuse_rtm_and_thlm, &
        l_stability_correct_Kh_N2_zm, &
        l_stability_correct_tau_zm, &
        l_do_expldiff_rtm_thlm, &
        l_Lscale_plume_centered, &
        l_use_ice_latent, &
        l_use_C7_Richardson, &
        l_use_C11_Richardson, &
        l_brunt_vaisala_freq_moist, &
        l_use_thvm_in_bv_freq, &
        l_rcm_supersat_adj, &
        ! The parameters of CLUBB can be retrieved and tuned using these indices:
        iC1, iC1b, iC1c, iC2, iC2b, iC2c,  &
        iC2rt, iC2thl, iC2rtthl, iC4, iC5, &
        iC6rt, iC6rtb, iC6rtc, iC6thl, iC6thlb, iC6thlc, &
        iC7, iC7b, iC7c, iC8, iC8b, iC10, iC11, iC11b, iC11c, &
        iC12, iC13, iC14, iC15, iC6rt_Lscale0, iC6thl_Lscale0, &
        iC7_Lscale0, iwpxp_L_thresh, ic_K, ic_K1, inu1, ic_K2, inu2, &
        ic_K6, inu6, ic_K8, inu8, ic_K9, inu9, inu10, ic_K_hm, ic_K_hmb, iK_hm_min_coef, &
        inu_hm, ibeta, igamma_coef, igamma_coefb, igamma_coefc, ilmin_coef, &
        iomicron, izeta_vrnce_rat, iupsilon_precip_frac_rat, &
        ilambda0_stability_coef, imult_coef, itaumin, itaumax, imu, iLscale_mu_coef, &
        iLscale_pert_coef, ialpha_corr, iSkw_denom_coef, ic_K10, ic_K10h, ithlp2_rad_coef, &
        ithlp2_rad_cloud_frac_thresh, iup2_vp2_factor



  public &
    advance_clubb_core_api, &
        pdf_parameter, &
        ! A hydromet array is required, and these variables are required for a hydromet array:
        hydromet_list, &
        hydromet_tol, &
        hydromet_dim, &
        iiNgm, &
        iiNim, &
        iiNrm, &
        iiNsm, &
        iirgm, &
        iirim, &
        iirrm, &
        iirsm, &
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
    setup_pdf_indices_api, &
    setup_corr_varnce_array_api, &
    setup_pdf_parameters_api, &
    hydromet_pdf_parameter, &
    ! lh_subcolumn_generator - SILHS API
    genrand_init_api, & ! if you are doing restarts)
    genrand_state, &
    genrand_srepr, &
    genrand_intg, &
    ! To use the results, you will need these variables:
    corr_array_n_cloud, &
    corr_array_n_below, &
    d_variables,        &
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
    hmp2_ip_on_hmm2_ip_ratios_type

  public &
    ! To Interact With CLUBB's Grid:
    gr, &
    ! For Varying Grids
    setup_grid_heights_api    ! if heights vary with time

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
    calculate_thlp2_rad_api, mu, update_xp2_mc_api, sat_mixrat_liq_api

  public :: &
    ! To Convert Between Common CLUBB-related quantities:
    lin_interpolate_two_points_api, & ! OR
    lin_interpolate_on_grid_api, &
    T_in_K2thlm_api, &
    thlm2T_in_K_api, &
    zt2zm_api, &
    zm2zt_api

  public &
    ! To Check For and Handle CLUBB's Errors:
    calculate_spurious_source_api, &
    clubb_at_least_debug_level_api, &
    clubb_no_error, &
    fatal_error_api, &
    fill_holes_driver_api, & ! OR
    fill_holes_vertical_api, &
    report_error_api, &
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
#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
    pack_pdf_params_api, &
    unpack_pdf_params_api, &
    num_pdf_params, &
#endif
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
    l_prescribed_avg_deltaz, &
    leap_year_api, &
    Lscale, &
    nvarmax_rad_zm, &
    nvarmax_rad_zt, &
    nvarmax_sfc, &
    nvarmax_zm, &
    nvarmax_zt, &
    stats_rad_zm, &
    stats_rad_zt
    public &
    nparams, &
    setup_parameters_api, &
    stats_sfc, &
    stat_nknd, &
    stat_rknd, &
    stats_accumulate_hydromet_api, &
    stats_init_rad_zm_api, &
    stats_init_rad_zt_api, &
    stats_init_sfc_api, &
    stats_init_zm_api, &
    stats_init_zt_api, &
    wp2_zt, &
    wphydrometp, &
    stats_zm, &
    zmscr01, zmscr02, zmscr03, &
    zmscr04, zmscr05, zmscr06, &
    zmscr07, zmscr08, zmscr09, &
    zmscr10, zmscr11, zmscr12, &
    zmscr13, zmscr14, zmscr15, &
    zmscr16, zmscr17, &
    stats_zt, &
    ztscr01, ztscr02, ztscr03, &
    ztscr04, ztscr05, ztscr06, &
    ztscr07, ztscr08, ztscr09, &
    ztscr10, ztscr11, ztscr12, &
    ztscr13, ztscr14, ztscr15, &
    ztscr16, ztscr17, ztscr18, &
    ztscr19, ztscr20, ztscr21

  interface zt2zm_api
    module procedure zt2zm_scalar_api, zt2zm_prof_api
  end interface

  interface zm2zt_api
    module procedure zm2zt_scalar_api, zm2zt_prof_api
  end interface

contains

  !================================================================================================
  ! advance_clubb_core - Advances the model one timestep.
  !================================================================================================

  subroutine advance_clubb_core_api( &
    l_implemented, dt, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
    thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! intent(in)
    sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! intent(in)
    wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! intent(in)
    rtpthlp_forcing, wm_zm, wm_zt, &                        ! intent(in)
    wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &            ! intent(in)
    wpsclrp_sfc, wpedsclrp_sfc, &                           ! intent(in)
    p_in_Pa, rho_zm, rho, exner, &                          ! intent(in)
    rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! intent(in)
    invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &      ! intent(in)
    rfrzm, radf, &                                          ! intent(in)
#ifdef CLUBBND_CAM
    varmu, &                                                ! intent(in)
#endif
    wphydrometp, wp2hmp, rtphmp, thlphmp, &    ! intent(in)
    host_dx, host_dy, &                                     ! intent(in)
    um, vm, upwp, vpwp, up2, vp2, &                         ! intent(inout)
    thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
    wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! intent(inout)
    sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16    ! intent(inout)
#endif
    sclrp2, sclrprtp, sclrpthlp, &                          ! intent(inout)
    wpsclrp, edsclrm, err_code, &                           ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                  ! intent(inout)
               do_liquid_only_in_clubb, &                   ! intent(in)
#endif
    rcm, wprcp, cloud_frac, ice_supersat_frac, &            ! intent(out)
    rcm_in_layer, cloud_cover, &                            ! intent(out)
#if defined(CLUBB_CAM) || defined(GFDL)
    khzm, khzt, &                                           ! intent(out)
#endif
#ifdef CLUBB_CAM
    qclvar, thlprcp_out, &                                  ! intent(out)
#endif
    pdf_params )                                            ! intent(out)

    use advance_clubb_core_module, only : advance_clubb_core

    use parameters_model, only: &
      sclr_dim, & ! Variable(s)
      edsclr_dim

    implicit none
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
      p_in_Pa,         & ! Air pressure (thermodynamic levels)       [Pa]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      exner,           & ! Exner function (thermodynamic levels)     [-]
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

    real( kind = core_rknd ), intent(in),  dimension(sclr_dim) ::  &
      wpsclrp_sfc      ! Scalar flux at surface         [{units vary} m/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), intent(in),  dimension(edsclr_dim) ::  &
      wpedsclrp_sfc    ! Eddy-Scalar flux at surface    [{units vary} m/s]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), intent(in) :: &
      host_dx,  & ! East-West horizontal grid spacing     [m]
      host_dy     ! North-South horizontal grid spacing   [m]


    !!! Input/Output Variables
    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  &
      um,      & ! u mean wind component (thermodynamic levels)   [m/s]
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vm,      & ! v mean wind component (thermodynamic levels)   [m/s]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
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
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

#ifdef GFDL
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

      real( kind = core_rknd ), intent(inout), dimension(gr%nz,edsclr_dim) :: &
      edsclrm   ! Eddy passive scalar mean (thermo. levels)   [units vary]

    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  &
      rcm,          & ! cloud water mixing ratio, r_c (thermo. levels)  [kg/kg]
      rcm_in_layer, & ! rcm in cloud layer                              [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    type(pdf_parameter), dimension(gr%nz), intent(out) :: &
      pdf_params      ! PDF parameters   [units vary]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  &
      wprcp,            & ! w'r_c' (momentum levels)                  [(kg/kg) m/s]
      cloud_frac,       & ! cloud fraction (thermodynamic levels)     [-]
      ice_supersat_frac   ! ice cloud fraction (thermodynamic levels) [-]

#if defined(CLUBB_CAM) || defined(GFDL)
    real( kind = core_rknd ), intent(out), dimension(gr%nz) :: &
      khzt, &       ! eddy diffusivity on thermo levels
      khzm          ! eddy diffusivity on momentum levels
#endif

#ifdef CLUBB_CAM
    real( kind = core_rknd), intent(out), dimension(gr%nz) :: &
      qclvar, &     ! cloud water variance
      thlprcp_out
#endif

      !!! Output Variable
      integer, intent(inout) :: err_code ! Diagnostic, for if some calculation goes amiss.

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), intent(inOUT), dimension(gr%nz, min(1,sclr_dim) , 2) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
    logical, intent(in)                 ::  do_liquid_only_in_clubb
#endif
    call advance_clubb_core( &
      l_implemented, dt, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
      thlm_forcing, rtm_forcing, um_forcing, vm_forcing, &    ! intent(in)
      sclrm_forcing, edsclrm_forcing, wprtp_forcing, &        ! intent(in)
      wpthlp_forcing, rtp2_forcing, thlp2_forcing, &          ! intent(in)
      rtpthlp_forcing, wm_zm, wm_zt, &                        ! intent(in)
      wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &            ! intent(in)
      wpsclrp_sfc, wpedsclrp_sfc, &                           ! intent(in)
      p_in_Pa, rho_zm, rho, exner, &                          ! intent(in)
      rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &                ! intent(in)
      invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &      ! intent(in)
      rfrzm, radf, &                                          ! intent(in)
#ifdef CLUBBND_CAM
      varmu, &
#endif
      wphydrometp, wp2hmp, rtphmp, thlphmp, &                 ! intent(in)
      host_dx, host_dy, &                                     ! intent(in)
      um, vm, upwp, vpwp, up2, vp2, &                         ! intent(inout)
      thlm, rtm, wprtp, wpthlp, &                             ! intent(inout)
      wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &          ! intent(inout)
      sclrm,   &
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16               ! intent(inout)
#endif
      sclrp2, sclrprtp, sclrpthlp, &                          ! intent(inout)
      wpsclrp, edsclrm, err_code, &                           ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                             ! intent(inout)
               do_liquid_only_in_clubb, &                              ! intent(in)
#endif
      rcm, wprcp, cloud_frac, ice_supersat_frac, &            ! intent(out)
      rcm_in_layer, cloud_cover, &                            ! intent(out)
#if defined(CLUBB_CAM) || defined(GFDL)
               khzm, khzt, &                                           ! intent(out)
#endif
#ifdef CLUBB_CAM
               qclvar, thlprcp_out, &                                               ! intent(out)
#endif
      pdf_params )                                            ! intent(out)
  end subroutine advance_clubb_core_api

  !================================================================================================
  ! setup_clubb_core - Sets up the model for execution.
  !================================================================================================

  subroutine setup_clubb_core_api( &
    nzmax, T0_in, ts_nudge_in,              & ! intent(in)
    hydromet_dim_in, sclr_dim_in,           & ! intent(in)
    sclr_tol_in, edsclr_dim_in, params,     & ! intent(in)
    l_host_applies_sfc_fluxes,              & ! intent(in)
    l_uv_nudge, saturation_formula,         & ! intent(in)
#ifdef GFDL
      I_sat_sphum,                                       & ! intent(in)  h1g, 2010-06-16
#endif
    l_implemented, grid_type, deltaz, zm_init, zm_top, & ! intent(in)
    momentum_heights, thermodynamic_heights,           & ! intent(in)
    sfc_elevation,                                     & ! intent(in)
#ifdef GFDL
      cloud_frac_min ,                                   & ! intent(in)  h1g, 2010-06-16
#endif
    err_code )                                           ! intent(out)

    use advance_clubb_core_module, only : setup_clubb_core

    use parameter_indices, only:  &
      nparams ! Variable(s)

! TODO: This should be called from the api, but all the host models appear to call
!       it directly or not at all.
!   use model_flags, only: &
!     setup_model_flags    ! Subroutine

#ifdef MKL
      use csr_matrix_class, only: &
        initialize_csr_class, & ! Subroutine
        intlc_5d_5d_ja_size     ! Variable

#endif

      implicit none

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
      l_uv_nudge,             & ! Wind nudging
      l_host_applies_sfc_fluxes ! Whether to apply for the surface flux

    character(len=*), intent(in) :: &
      saturation_formula ! Approximation for saturation vapor pressure

#ifdef GFDL
      logical, intent(in) :: &  ! h1g, 2010-06-16 begin mod
         I_sat_sphum

      real( kind = core_rknd ), intent(in) :: &
         cloud_frac_min         ! h1g, 2010-06-16 end mod
#endif

      ! Output variables
      integer, intent(out) :: &
      err_code   ! Diagnostic for a problem with the setup

    call setup_clubb_core &
      ( nzmax, T0_in, ts_nudge_in,              & ! intent(in)
      hydromet_dim_in, sclr_dim_in,           & ! intent(in)
      sclr_tol_in, edsclr_dim_in, params,     & ! intent(in)
      l_host_applies_sfc_fluxes,              & ! intent(in)
      l_uv_nudge, saturation_formula,         & ! intent(in)
#ifdef GFDL
      I_sat_sphum,                                       & ! intent(in)  h1g, 2010-06-16
#endif
      l_implemented, grid_type, deltaz, zm_init, zm_top, & ! intent(in)
      momentum_heights, thermodynamic_heights,           & ! intent(in)
      sfc_elevation,                                     & ! intent(in)
#ifdef GFDL
      cloud_frac_min ,                                   & ! intent(in)  h1g, 2010-06-16
#endif
      err_code )                                           ! intent(out)

  end subroutine setup_clubb_core_api

  !================================================================================================
  ! cleanup_clubb_core_api - Frees memory used by the model.
  !================================================================================================

  subroutine cleanup_clubb_core_api( &
    l_implemented )

    use advance_clubb_core_module, only : cleanup_clubb_core

    implicit none

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, intent(in) :: l_implemented   ! (T/F)

    call cleanup_clubb_core( &
    l_implemented )
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
      previous_day, previous_month, &
      previous_year,  &
      seconds_since_previous_date, &
      current_day, current_month, &
      current_year, &
      seconds_since_current_date )
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
    input_file_cloud, input_file_below, iunit )

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

    call setup_corr_varnce_array( &
      input_file_cloud, input_file_below, iunit )

  end subroutine setup_corr_varnce_array_api

  !================================================================================================
  ! setup_pdf_indices - Sets up the iiPDF indices.
  !================================================================================================

  subroutine setup_pdf_indices_api( &
    hydromet_dim, iirrm, iiNrm, &
    iirim, iiNim, iirsm, iiNsm, &
    iirgm, iiNgm )

    use corr_varnce_module, only : setup_pdf_indices

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim    ! Total number of hydrometeor species.

    integer, intent(in) :: &
      iirrm,    & ! Index of rain water mixing ratio
      iiNrm,       & ! Index of rain drop concentration
      iirim,     & ! Index of ice mixing ratio
      iiNim,       & ! Index of ice crystal concentration
      iirsm,    & ! Index of snow mixing ratio
      iiNsm,    & ! Index of snow flake concentration
      iirgm, & ! Index of graupel mixing ratio
      iiNgm    ! Index of graupel concentration

    call setup_pdf_indices( &
      hydromet_dim, iirrm, iiNrm, &
      iirim, iiNim, iirsm, iiNsm, &
      iirgm, iiNgm )
  end subroutine setup_pdf_indices_api

  !================================================================================================
  ! report_error - Reports the meaning of an error code to the console.
  !================================================================================================

  subroutine report_error_api( &
    err_code)

    use error_code, only: &
      report_error  ! Procedure

    implicit none

    ! Input Variable
    integer, intent(in) :: err_code ! Error Code being examined

    call report_error( &
      err_code)
  end subroutine report_error_api

  !================================================================================================
  ! fatal_error - Checks to see if an error code is usually one which causes an exit elsewhere.
  !================================================================================================

  elemental function fatal_error_api( &
    err_code )

    use error_code, only : fatal_error

    implicit none

    ! Input Variable
    integer, intent(in) :: err_code ! Error Code being examined

    ! Output variable
    logical :: fatal_error_api

    fatal_error_api = fatal_error( &
      err_code )
  end function fatal_error_api

  !================================================================================================
  ! set_clubb_debug_level - Controls the importance of error messages sent to the console.
  !================================================================================================

  subroutine set_clubb_debug_level_api( &
    level )

    use error_code, only : set_clubb_debug_level

    implicit none

    ! Input variable
    integer, intent(in) :: level ! The debug level being checked against the current setting

    call set_clubb_debug_level( &
      level )
  end subroutine set_clubb_debug_level_api

  !================================================================================================
  ! clubb_at_least_debug_level - Checks to see if clubb has been set to a specified debug level.
  !================================================================================================

  logical function clubb_at_least_debug_level_api( &
    level )

    use error_code, only : clubb_at_least_debug_level

    implicit none

    ! Input variable
    integer, intent(in) :: level   ! The debug level being checked against the current setting

    clubb_at_least_debug_level_api = clubb_at_least_debug_level( &
      level )
  end function clubb_at_least_debug_level_api

  !================================================================================================
  ! fill_holes_driver - Fills holes between same-phase hydrometeors(i.e. for frozen hydrometeors).
  !================================================================================================

  subroutine fill_holes_driver_api( &
    nz, dt, hydromet_dim,        & ! Intent(in)
    l_fill_holes_hm,             & ! Intent(in)
    rho_ds_zm, rho_ds_zt, exner, & ! Intent(in)
    thlm_mc, rvm_mc, hydromet )    ! Intent(inout)

    use fill_holes, only : fill_holes_driver

    use constants_clubb, only: &
      four_thirds,     &
      rho_ice

    implicit none

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

    call fill_holes_driver( &
      nz, dt, hydromet_dim,        & ! Intent(in)
      l_fill_holes_hm,             & ! Intent(in)
      rho_ds_zm, rho_ds_zt, exner, & ! Intent(in)
      thlm_mc, rvm_mc, hydromet )    ! Intent(inout)
  end subroutine fill_holes_driver_api

  !================================================================================================
  ! fill_holes_vertical - clips values of 'field' that are below 'threshold' as much as possible.
  !================================================================================================

  subroutine fill_holes_vertical_api( &
    num_pts, threshold, field_grid, &
    rho_ds, rho_ds_zm, &
    field )

    use fill_holes, only : fill_holes_vertical

    implicit none

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

    call fill_holes_vertical( &
      num_pts, threshold, field_grid, &
      rho_ds, rho_ds_zm, &
      field )
  end subroutine fill_holes_vertical_api

  !================================================================================================
  ! vertical_integral - Computes the vertical integral.
  !================================================================================================

  function vertical_integral_api( &
    total_idx, rho_ds, &
    field, invrs_dz )

    use fill_holes, only : vertical_integral

    implicit none

    ! Input variables
    integer, intent(in) :: &
      total_idx  ! The total numer of indices within the range of averaging

    real( kind = core_rknd ), dimension(total_idx), intent(in) ::  &
      rho_ds,  & ! Dry, static density                   [kg/m^3]
      field,   & ! The field to be vertically averaged   [Units vary]
      invrs_dz   ! Level thickness                       [1/m]
    ! Note:  The rho_ds and field points need to be arranged from
    !        lowest to highest in altitude, with rho_ds(1) and
    !        field(1) actually their respective values at level k = begin_idx.

    real( kind = core_rknd ) :: &
      vertical_integral_api ! Integral in the numerator (see description)

    vertical_integral_api = vertical_integral( &
      total_idx, rho_ds, &
      field, invrs_dz )
  end function vertical_integral_api

  !================================================================================================
  ! setup_grid_heights - Sets the heights and interpolation weights of the column.
  !================================================================================================

  subroutine setup_grid_heights_api( &
    l_implemented, grid_type,  &
    deltaz, zm_init, momentum_heights,  &
    thermodynamic_heights )

    use grid_class, only : setup_grid_heights

    implicit none

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
      l_implemented, grid_type,  &
      deltaz, zm_init, momentum_heights,  &
      thermodynamic_heights )

  end subroutine setup_grid_heights_api


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
      nparam, xlist, tlist, xvalue, tvalue )

  end subroutine lin_interpolate_on_grid_api

  !================================================================================================
  ! read_parameters - Read a namelist containing the model parameters.
  !================================================================================================

  subroutine read_parameters_api( &
    iunit, filename, params )

    use parameters_tunable, only : read_parameters

    use parameter_indices, only:  &
      nparams ! Variable(s)

    implicit none

    ! Input variables
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(nparams) :: params

    call read_parameters( &
      iunit, filename, params )

  end subroutine read_parameters_api

  !================================================================================================
  ! setup_parameters - Sets up model parameters.
  !================================================================================================

  subroutine setup_parameters_api( &
    deltaz, params, nzmax, &
    grid_type, momentum_heights, thermodynamic_heights, &
    err_code )

    use parameters_tunable, only: &
      setup_parameters

    use constants_clubb, only:  &
      fstderr ! Variable(s)

    use error_code, only:  &
      clubb_var_out_of_bounds ! Variable(s)

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

    ! Output Variables
    integer, intent(out) ::  &
      err_code ! Error condition

    call setup_parameters( &
      deltaz, params, nzmax, &
      grid_type, momentum_heights, thermodynamic_heights, &
      err_code )

  end subroutine setup_parameters_api

  !================================================================================================
  ! adj_low_res_nu - Adjusts values of background eddy diffusivity based on vertical grid spacing.
  !================================================================================================

  subroutine adj_low_res_nu_api( &
    nzmax, grid_type, deltaz, & ! Intent(in)
    momentum_heights, thermodynamic_heights )  ! Intent(in)

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

    call adj_low_res_nu( &
      nzmax, grid_type, deltaz, & ! Intent(in)
      momentum_heights, thermodynamic_heights )  ! Intent(in)
  end subroutine adj_low_res_nu_api

#ifdef CLUBB_CAM /* Code for storing pdf_parameter structs in pbuf as array */
  !================================================================================================
  ! pack_pdf_params - Returns a two dimensional real array with all values.
  !================================================================================================

  subroutine pack_pdf_params_api( &
    pdf_params, nz, r_param_array)

    use pdf_parameter_module, only : pack_pdf_params

    !use statements

    implicit none

    ! Input a pdf_parameter array with nz instances of pdf_parameter
    integer, intent(in) :: nz ! Num Vert Model Levs
    type (pdf_parameter), dimension(nz), intent(in) :: pdf_params

    ! Output a two dimensional real array with all values
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(out) :: &
      r_param_array

    call pack_pdf_params( &
      pdf_params, nz, r_param_array)

  end subroutine pack_pdf_params_api

  !================================================================================================
  ! unpack_pdf_params - Returns a pdf_parameter array with nz instances of pdf_parameter.
  !================================================================================================

  subroutine unpack_pdf_params_api( &
    r_param_array, nz, pdf_params)

    use pdf_parameter_module, only : unpack_pdf_params

    implicit none

    ! Input a two dimensional real array with pdf values
    integer, intent(in) :: nz ! Num Vert Model Levs
    real (kind = core_rknd), dimension(nz,num_pdf_params), intent(in) :: &
      r_param_array

    ! Output a pdf_parameter array with nz instances of pdf_parameter
    type (pdf_parameter), dimension(nz), intent(out) :: pdf_params

    call unpack_pdf_params( &
      r_param_array, nz, pdf_params)
  end subroutine unpack_pdf_params_api
#endif

  !================================================================================================
  ! setup_pdf_parameters
  !================================================================================================

  subroutine setup_pdf_parameters_api( &
    nz, d_variables, dt, &                      ! Intent(in)
    Nc_in_cloud, rcm, cloud_frac, &             ! Intent(in)
    ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
    corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
    pdf_params, l_stats_samp, &                 ! Intent(in)
    hydrometp2, &                               ! Intent(inout)
    mu_x_1_n, mu_x_2_n, &                       ! Intent(out)
    sigma_x_1_n, sigma_x_2_n, &                 ! Intent(out)
    corr_array_1_n, corr_array_2_n, &           ! Intent(out)
    corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
    hydromet_pdf_params )                       ! Intent(out)

    use setup_clubb_pdf_params, only : setup_pdf_parameters

    use constants_clubb, only: &
      one,            & ! Constant(s)
      Ncn_tol,        &
      cloud_frac_min

    use advance_windm_edsclrm_module, only: &
      xpwp_fnc

    use parameters_tunable, only: &
      c_K_hm

    use clip_explicit, only: &
      clip_wphydrometp    ! Variables(s)

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,          & ! Number of model vertical grid levels
      d_variables    ! Number of variables in the correlation array

    real( kind = core_rknd ), intent(in) ::  &
      dt    ! Model timestep                                           [s]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      Nc_in_cloud,       & ! Mean (in-cloud) cloud droplet conc.       [num/kg]
      rcm,               & ! Mean cloud water mixing ratio, < r_c >    [kg/kg]
      cloud_frac,        & ! Cloud fraction                            [-]
      ice_supersat_frac    ! Ice supersaturation fraction              [-]

    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(in) :: &
      hydromet,    & ! Mean of hydrometeor, hm (overall) (t-levs.) [units]
      wphydrometp    ! Covariance < w'h_m' > (momentum levels)     [(m/s)units]

    real( kind = core_rknd ), dimension(d_variables,d_variables), &
      intent(in) :: &
      corr_array_n_cloud, & ! Prescribed norm. space corr. array in cloud    [-]
      corr_array_n_below    ! Prescribed norm. space corr. array below cloud [-]

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params    ! PDF parameters                               [units vary]

    logical, intent(in) :: &
      l_stats_samp    ! Flag to sample statistics

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz,hydromet_dim), intent(inout) :: &
      hydrometp2    ! Variance of a hydrometeor (overall) (m-levs.)   [units^2]

    ! Output Variables
    real( kind = core_rknd ), dimension(d_variables, nz), intent(out) :: &
      mu_x_1_n,    & ! Mean array (normal space): PDF vars. (comp. 1) [un. vary]
      mu_x_2_n,    & ! Mean array (normal space): PDF vars. (comp. 2) [un. vary]
      sigma_x_1_n, & ! Std. dev. array (normal space): PDF vars (comp. 1) [u.v.]
      sigma_x_2_n    ! Std. dev. array (normal space): PDF vars (comp. 2) [u.v.]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
      intent(out) :: &
      corr_array_1_n, & ! Corr. array (normal space):  PDF vars. (comp. 1)   [-]
      corr_array_2_n    ! Corr. array (normal space):  PDF vars. (comp. 2)   [-]

    real( kind = core_rknd ), dimension(d_variables,d_variables,nz), &
      intent(out) :: &
      corr_cholesky_mtx_1, & ! Transposed corr. cholesky matrix, 1st comp. [-]
      corr_cholesky_mtx_2    ! Transposed corr. cholesky matrix, 2nd comp. [-]

    type(hydromet_pdf_parameter), dimension(nz), intent(out) :: &
      hydromet_pdf_params    ! Hydrometeor PDF parameters        [units vary]

    call setup_pdf_parameters( &
      nz, d_variables, dt, &                      ! Intent(in)
      Nc_in_cloud, rcm, cloud_frac, &             ! Intent(in)
      ice_supersat_frac, hydromet, wphydrometp, & ! Intent(in)
      corr_array_n_cloud, corr_array_n_below, &   ! Intent(in)
      pdf_params, l_stats_samp, &                 ! Intent(in)
      hydrometp2, &                               ! Intent(inout)
      mu_x_1_n, mu_x_2_n, &                       ! Intent(out)
      sigma_x_1_n, sigma_x_2_n, &                 ! Intent(out)
      corr_array_1_n, corr_array_2_n, &           ! Intent(out)
      corr_cholesky_mtx_1, corr_cholesky_mtx_2, & ! Intent(out)
      hydromet_pdf_params )                       ! Intent(out)

  end subroutine setup_pdf_parameters_api

  !================================================================================================
  ! stats_init - Initializes the statistics saving functionality of the CLUBB model.
  !================================================================================================

  subroutine stats_init_api( &
    iunit, fname_prefix, fdir, l_stats_in, &
    stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
    nzmax, nlon, nlat, gzt, gzm, nnrad_zt, &
    grad_zt, nnrad_zm, grad_zm, day, month, year, &
    rlon, rlat, time_current, delt, l_silhs_out_in )

    use stats_clubb_utilities, only : stats_init

    implicit none

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
      rlon  ! Longitude(s) [Degrees E]

    real( kind = core_rknd ), dimension(nlat), intent(in) ::  &
      rlat  ! Latitude(s)  [Degrees N]

    real( kind = time_precision ), intent(in) ::  &
      time_current ! Model time                         [s]

    real( kind = core_rknd ), intent(in) ::  &
      delt         ! Timestep (dt_main in CLUBB)         [s]

    logical, intent(in) :: &
      l_silhs_out_in  ! Whether to output SILHS files (stats_lh_zt,stats_lh_sfc) [dimensionless]

    call stats_init( &
      iunit, fname_prefix, fdir, l_stats_in, &
      stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
      nzmax, nlon, nlat, gzt, gzm, nnrad_zt, &
      grad_zt, nnrad_zm, grad_zm, day, month, year, &
      rlon, rlat, time_current, delt, l_silhs_out_in )
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
      itime, stats_nsamp, stats_nout )
  end subroutine stats_begin_timestep_api

  !================================================================================================
  ! stats_end_timestep - Calls statistics to be written to the output format.
  !================================================================================================

  subroutine stats_end_timestep_api

    use stats_clubb_utilities, only : stats_end_timestep

    implicit none

    call stats_end_timestep

  end subroutine stats_end_timestep_api

  !================================================================================================
  ! stats_accumulate_hydromet - Computes stats related the hydrometeors.
  !================================================================================================

  subroutine stats_accumulate_hydromet_api( &
    hydromet, rho_ds_zt )

    use stats_clubb_utilities, only : stats_accumulate_hydromet

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet ! All hydrometeors except for rcm        [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      rho_ds_zt ! Dry, static density (thermo. levs.)      [kg/m^3]

    call stats_accumulate_hydromet( &
      hydromet, rho_ds_zt )
  end subroutine stats_accumulate_hydromet_api

  !================================================================================================
  ! stats_finalize - Close NetCDF files and deallocate scratch space and stats file structures.
  !================================================================================================

  subroutine stats_finalize_api

    use stats_clubb_utilities, only : stats_finalize

    implicit none

    call stats_finalize

  end subroutine stats_finalize_api

  !================================================================================================
  ! stats_init_rad_zm - Initializes array indices for rad_zm variables.
  !================================================================================================

  subroutine stats_init_rad_zm_api( &
    vars_rad_zm, l_error )

    use stats_rad_zm_module, only : stats_init_rad_zm, nvarmax_rad_zm

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zm), intent(in) :: vars_rad_zm

    ! Input/Output Variable
    logical, intent(inout) :: l_error

    call stats_init_rad_zm( &
      vars_rad_zm, l_error )
  end subroutine stats_init_rad_zm_api

  !================================================================================================
  ! stats_init_rad_zt - Initializes array indices for zt.
  !================================================================================================

  subroutine stats_init_rad_zt_api( &
    vars_rad_zt, l_error )

    use stats_rad_zt_module, only : stats_init_rad_zt, nvarmax_rad_zt

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_rad_zt), intent(in) :: vars_rad_zt

    ! Input/Output Variable
    logical, intent(inout) :: l_error

    call stats_init_rad_zt( &
      vars_rad_zt, l_error )
  end subroutine stats_init_rad_zt_api

  !================================================================================================
  ! stats_init_zm - Initializes array indices for zm.
  !================================================================================================

  subroutine stats_init_zm_api( &
    vars_zm, l_error )

    use stats_zm_module, only : stats_init_zm, nvarmax_zm

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_zm), intent(in) :: vars_zm ! zm variable names

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    call stats_init_zm( &
      vars_zm, l_error )

  end subroutine stats_init_zm_api

  !================================================================================================
  ! stats_init_zt - Initializes array indices for zt.
  !================================================================================================

  subroutine stats_init_zt_api( &
    vars_zt, l_error )

    use stats_zt_module, only : stats_init_zt, nvarmax_zt

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_zt), intent(in) :: vars_zt

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    call stats_init_zt( &
      vars_zt, l_error )

  end subroutine stats_init_zt_api

  !================================================================================================
  ! stats_init_sfc - Initializes array indices for sfc.
  !================================================================================================

  subroutine stats_init_sfc_api( &
    vars_sfc, l_error )

    use stats_sfc_module, only : stats_init_sfc, nvarmax_sfc

    implicit none

    ! Input Variable
    character(len= * ), dimension(nvarmax_sfc), intent(in) :: vars_sfc

    ! Input / Output Variable
    logical, intent(inout) :: l_error

    call stats_init_sfc( &
      vars_sfc, l_error )

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
  function zm2zt_scalar_api( azm, k )

    use grid_class, only: zm2zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    integer, intent(in) :: &
      k      ! Vertical level index

    ! Return Variable
    real( kind = core_rknd ) :: &
      zm2zt_scalar_api   ! Variable when interp. to thermo. levels

    zm2zt_scalar_api = zm2zt( azm, k )

  end function zm2zt_scalar_api

  !================================================================================================
  ! zt2zm_scalar - Interpolates a variable from zt to zm grid at one height level
  !================================================================================================
  function zt2zm_scalar_api( azt, k )

    use grid_class, only: zt2zm

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    integer, intent(in) :: &
      k      ! Vertical level index

    ! Return Variable
    real( kind = core_rknd ) :: &
      zt2zm_scalar_api   ! Variable when interp. to momentum levels

    zt2zm_scalar_api = zt2zm( azt, k )

  end function zt2zm_scalar_api

  !================================================================================================
  ! zt2zm_prof - Interpolates a variable (profile) from zt to zm grid
  !================================================================================================
  function zt2zm_prof_api( azt )

    use grid_class, only: zt2zm

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azt    ! Variable on thermodynamic grid levels    [units vary]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      zt2zm_prof_api   ! Variable when interp. to momentum levels

    zt2zm_prof_api = zt2zm( azt )

  end function zt2zm_prof_api

  !================================================================================================
  ! zm2zt_prof - Interpolates a variable (profile) from zm to zt grid
  !================================================================================================
  function zm2zt_prof_api( azm )

    use grid_class, only: zm2zt

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) :: &
      azm    ! Variable on momentum grid levels    [units vary]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      zm2zt_prof_api   ! Variable when interp. to thermo. levels

    zm2zt_prof_api = zm2zt( azm )

  end function zm2zt_prof_api

  !================================================================================================
  ! calculate_thlp2_rad - Computes the contribution of radiative cooling to thlp2
  !================================================================================================
  pure subroutine calculate_thlp2_rad_api &
                  ( nz, rcm_zm, thlprcp, radht_zm, &      ! Intent(in)
                    thlp2_forcing )                       ! Intent(inout)

  use clubb_precision, only: &
    core_rknd                     ! Constant(s)

  use advance_clubb_core_module, only: &
    calculate_thlp2_rad

  implicit none

  ! Input Variables
  integer, intent(in) :: &
    nz                    ! Number of vertical levels                      [-]

  real( kind = core_rknd ), dimension(nz), intent(in) :: &
    rcm_zm, &             ! Cloud water mixing ratio on momentum grid      [kg/kg]
    thlprcp, &            ! thl'rc'                                        [K kg/kg]
    radht_zm              ! SW + LW heating rate (on momentum grid)        [K/s]

  ! Input/Output Variables
  real( kind = core_rknd ), dimension(nz), intent(inout) :: &
    thlp2_forcing         ! <th_l'^2> forcing (momentum levels)            [K^2/s]
  !----------------------------------------------------------------------

    call calculate_thlp2_rad( nz, rcm_zm, thlprcp, radht_zm, &
                    thlp2_forcing )

    return
  end subroutine calculate_thlp2_rad_api

  !================================================================================================
  ! update_xp2_mc - Calculates the effects of rain evaporation on rtp2 and thlp2
  !================================================================================================
  subroutine update_xp2_mc_api( nz, dt, cloud_frac, rcm, rvm, thlm,        &
                            wm, exner, rrm_evap, pdf_params,        &
                            rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc,    &
                            rtpthlp_mc )

    use advance_xp2_xpyp_module, only: &
      update_xp2_mc

    implicit none

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

    type(pdf_parameter), dimension(nz), intent(in) :: &
      pdf_params ! PDF parameters

    !input/output variables
    real( kind = core_rknd ), dimension(nz), intent(inout) :: &
      rtp2_mc, &    !Tendency of <rt'^2> due to evaporation   [(kg/kg)^2/s]
      thlp2_mc, &   !Tendency of <thl'^2> due to evaporation  [K^2/s]
      wprtp_mc, &   !Tendency of <w'rt'> due to evaporation   [m*(kg/kg)/s^2]
      wpthlp_mc, &  !Tendency of <w'thl'> due to evaporation  [m*K/s^2] 
      rtpthlp_mc    !Tendency of <rt'thl'> due to evaporation [K*(kg/kg)/s]

    call update_xp2_mc( nz, dt, cloud_frac, rcm, rvm, thlm,        &
                        wm, exner, rrm_evap, pdf_params,        &
                        rtp2_mc, thlp2_mc, wprtp_mc, wpthlp_mc,    &
                        rtpthlp_mc )
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


end module clubb_api_module
