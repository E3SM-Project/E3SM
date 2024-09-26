!-----------------------------------------------------------------------
! $Id$
!-----------------------------------------------------------------------

module advance_clubb_core_module

! Description:
!   The module containing the `core' of the CLUBB parameterization.
!   It advances CLUBB's equations one model time step.
!
! References:
! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_clubb
!
!  ``A PDF-Based Model for Boundary Layer Clouds. Part I:
!    Method and Model Description'' Golaz, et al. (2002)
!    JAS, Vol. 59, pp. 3540--3551.
!
!                         Copyright Notice:
!
!   This code and the source code it references are (C) 2006-2018.
!
!   The distribution of this code and derived works thereof
!                   should include this notice.
!
!   Portions of this code derived from other sources (Hugh Morrison,
!   ACM TOMS, Numerical Recipes, et cetera) are the intellectual
!   property of their respective authors as noted and are also subject
!   to copyright.
!
!
!
! Cloud Layers Unified By Binormals (CLUBB) user license
! agreement.
!
! Thank you for your interest in CLUBB. We work hard to create a
! code that implements the best software engineering practices,
! is supported to the extent allowed by our limited resources,
! and is available without cost to non-commercial users. You may
! use CLUBB if, in return, you abide by these conditions:
!
! 1. Please cite CLUBB in presentations and publications that
!  contain results obtained using CLUBB.
!
! 2. You may not use any part of CLUBB to create or modify
!  another single-column (1D) model that is not called CLUBB.
!  However, you may modify or augment CLUBB or parts of CLUBB if
!  you include "CLUBB" in the name of the resulting single-column
!  model. For example, a user at MIT might modify CLUBB and call
!  the modified version "CLUBB-MIT." Or, for example, a user of
!  the CLM land-surface model might interface CLM to CLUBB and
!  call it "CLM-CLUBB." This naming convention recognizes the
!  contributions of both sets of developers.
!
! 3. You may implement CLUBB as a parameterization in a large-
!  scale host model that has 2 or 3 spatial dimensions without
!  including "CLUBB" in the combined model name, but please
!  acknowledge in presentations and publications that CLUBB has
!  been included as a parameterization.
!
! 4. You may not provide all or part of CLUBB to anyone without
!  prior permission from Vincent Larson (vlarson@uwm.edu). If
!  you wish to share CLUBB with your collaborators without
!  seeking permission, please ask your collaborators to register
!  as CLUBB users at https://carson.math.uwm.edu/larson-group/clubb_site/ and to
!  download CLUBB from there.
!
! 5. You may not use CLUBB for commercial purposes unless you
!  receive permission from Vincent Larson.
!
! 6. You may not re-license all or any part of CLUBB.
!
! 7. CLUBB is provided "as is" and without warranty.
!
! We hope that CLUBB will develop into a community resource. We
! encourage users to contribute their CLUBB modifications or
! extensions to the CLUBB development group. We will then
! consider them for inclusion in CLUBB. Such contributions will
! benefit all CLUBB users. We would be pleased to acknowledge
! contributors and list their CLUBB-related papers on our "About
! CLUBB" webpage (https://carson.math.uwm.edu/larson-group/clubb_site/about.html) for
! those contributors who so desire.
!
! Thanks so much and best wishes for your research!
!
! The CLUBB Development Group
! (Present and past contributors to the source code include
! Vincent Larson, Chris Golaz, David Schanen, Brian Griffin,
! Joshua Fasching, Adam Smith, and Michael Falk).
!-----------------------------------------------------------------------

  use model_flags, only: ipdf_call_placement 

  implicit none

  public ::  &
    setup_clubb_core, &
    advance_clubb_core, &
    cleanup_clubb_core, &
    set_Lscale_max, &
    calculate_thlp2_rad

  ! Options for the placement of the call to CLUBB's PDF.
  integer, parameter :: &
    ipdf_pre_advance_fields   = 1,   & ! Call before advancing predictive fields
    ipdf_post_advance_fields  = 2,   & ! Call after advancing predictive fields
    ipdf_pre_post_advance_fields = 3   ! Call both before and after advancing
                                       ! predictive fields

  public :: ipdf_post_advance_fields

  private ! Default Scope

  contains

  !-----------------------------------------------------------------------

  !#######################################################################
  !#######################################################################
  ! If you change the argument list of advance_clubb_core you also have to
  ! change the calls to this function in the host models CAM, WRF, SAM
  ! and GFDL.
  !#######################################################################
  !#######################################################################
  subroutine advance_clubb_core &
             ( l_implemented, dt, fcor, sfc_elevation, hydromet_dim, & ! intent(in)
               thlm_forcing, rtm_forcing, um_forcing, vm_forcing, & ! intent(in)
               sclrm_forcing, edsclrm_forcing, wprtp_forcing, &     ! intent(in)
               wpthlp_forcing, rtp2_forcing, thlp2_forcing, &       ! intent(in)
               rtpthlp_forcing, wm_zm, wm_zt, &                     ! intent(in)
               wpthlp_sfc, wprtp_sfc, upwp_sfc, vpwp_sfc, &         ! intent(in)
               wpsclrp_sfc, wpedsclrp_sfc, &                        ! intent(in)
               p_in_Pa, rho_zm, rho, exner, &                       ! intent(in)
               rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &             ! intent(in)
               invrs_rho_ds_zt, thv_ds_zm, thv_ds_zt, hydromet, &   ! intent(in)
               rfrzm, radf, &                                       ! intent(in)
#ifdef CLUBBND_CAM
               varmu, &                                             ! intent(in)
#endif
               wphydrometp, wp2hmp, rtphmp_zt, thlphmp_zt, &        ! intent(in)
               host_dx, host_dy, &                                  ! intent(in)
               um, vm, upwp, vpwp, up2, vp2, &                      ! intent(inout)
               thlm, rtm, wprtp, wpthlp, &                          ! intent(inout)
               wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp, &       ! intent(inout)
               sclrm,   &                                           ! intent(inout)
#ifdef GFDL
               sclrm_trsport_only,  &  ! h1g, 2010-06-16            ! intent(inout)
#endif
               sclrp2, sclrprtp, sclrpthlp, &                       ! intent(inout)
               wpsclrp, edsclrm, &                                  ! intent(inout)
               rcm, cloud_frac, &                                   ! intent(inout)
               wpthvp, wp2thvp, rtpthvp, thlpthvp, &                ! intent(inout)
               sclrpthvp, &                                         ! intent(inout)
               pdf_params, pdf_params_zm, &                         ! intent(inout)
#ifdef GFDL
               RH_crit, & !h1g, 2010-06-16                          ! intent(inout)
               do_liquid_only_in_clubb, &                           ! intent(in)
#endif
#if defined(CLUBB_CAM) || defined(GFDL)
               khzm, khzt, &                                        ! intent(out)
#endif
#ifdef CLUBB_CAM
               qclvar, thlprcp_out, &                               ! intent(out)
#endif
               wprcp, ice_supersat_frac, &                          ! intent(out)
               rcm_in_layer, cloud_cover, &                         ! intent(out)
               upwp_sfc_pert, vpwp_sfc_pert, &                      ! intent(in)
               um_pert, vm_pert, upwp_pert, vpwp_pert )             ! intent(inout)

    ! Description:
    !   Subroutine to advance CLUBB one timestep

    ! References:
    !   https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:overview_clubb
    !
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !     Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !-----------------------------------------------------------------------

    ! Modules to be included

    use constants_clubb, only: &
      em_min, &
      thl_tol, &
      rt_tol, &
      w_tol, &
      w_tol_sqd, &
      ep2, &
      Cp, &
      Lv, &
      ep1, &
      fstderr, &
      zero_threshold, &
      three_halves, &
      one_fourth, &
      one, &
      unused_var, &
      grav, &
      vonk, &
      eps

    use parameters_tunable, only: &
      taumax, & ! Variable(s)
      c_K, &
      mu, &
      Lscale_mu_coef, &
      Lscale_pert_coef, &
      gamma_coef,  &
      gamma_coefb, &
      gamma_coefc, &
      c_K10, &
      c_K10h, &
      C1, C14, &
      C5, C4, &
      C_wp2_splat,&
      C_invrs_tau_bkgnd,&
      C_invrs_tau_sfc,&
      C_invrs_tau_shear,&
      C_invrs_tau_N2

    use parameters_model, only: &
        sclr_dim, & ! Variable(s)
        edsclr_dim, &
        T0, &
        sclr_tol

    use model_flags, only: &
        l_tke_aniso, &  ! Variable(s)
        l_call_pdf_closure_twice, &
        l_host_applies_sfc_fluxes, &
        l_stability_correct_tau_zm, &
        l_do_expldiff_rtm_thlm, &
        l_Lscale_plume_centered, &
        l_use_ice_latent, &
        l_gamma_Skw, &
        l_damp_wp2_using_em, &
        l_advance_xp3, &
        l_predict_upwp_vpwp, &
        l_diag_Lscale_from_tau, &
        l_use_C7_Richardson, &
        l_use_C11_Richardson, &
        l_use_wp3_pr3

    use grid_class, only: &
      gr,  & ! Variable(s)
      zm2zt,  & ! Procedure(s)
      zt2zm, &
      ddzm, &
      ddzt

    use numerical_check, only: &
      parameterization_check, & ! Procedure(s)
      calculate_spurious_source

    use variables_diagnostic_module, only: &
      Skw_zt,  & ! Variable(s)
      Skw_zm, &
      wp4, &
      rtprcp, &
      thlprcp, &
      rcp2, &
      rsat, &
      wprtp2, &
      wp2rtp, &
      wpthlp2, &
      wp2thlp, &
      wprtpthlp, &
      wp2rcp

    use variables_diagnostic_module, only: &
      thvm, &
      em, &
      Lscale, &
      Lscale_up,  &
      Lscale_down, &
      tau_zm, &
      tau_zt, &
      Kh_zm, &
      Kh_zt, &
      vg, &
      ug, &
      um_ref, &
      vm_ref

    use variables_diagnostic_module, only: &
      wp2_zt, &
      thlp2_zt, &
      wpthlp_zt, &
      wprtp_zt, &
      rtp2_zt, &
      rtpthlp_zt, &
      up2_zt, &
      vp2_zt, &
      upwp_zt, &
      vpwp_zt, &
      rtm_ref, &
      thlm_ref

    use variables_diagnostic_module, only: &
      wpedsclrp, &
      sclrprcp,     & ! sclr'rc'
      wp2sclrp,     & ! w'^2 sclr'
      wpsclrp2,     & ! w'sclr'^2
      wpsclrprtp,   & ! w'sclr'rt'
      wpsclrpthlp,  & ! w'sclr'thl'
      wp3_zm,       & ! wp3 interpolated to momentum levels
      Skw_velocity, & ! Skewness velocity       [m/s]
      a3_coef,      & ! The a3 coefficient      [-]
      a3_coef_zt      ! The a3 coefficient interp. to the zt grid [-]

    use variables_diagnostic_module, only: &
      wp3_on_wp2,   & ! Variable(s)
      wp3_on_wp2_zt

    use pdf_parameter_module, only: &
        pdf_parameter, &
        implicit_coefs_terms

#ifdef GFDL
    use advance_sclrm_Nd_module, only: &  ! h1g, 2010-06-16 begin mod
       advance_sclrm_Nd_diffusion_OG, &
       advance_sclrm_Nd_upwind, &
       advance_sclrm_Nd_semi_implicit     ! h1g, 2010-06-16 end mod
#endif

    use advance_xm_wpxp_module, only: &
        advance_xm_wpxp          ! Compute mean/flux terms

    use advance_xp2_xpyp_module, only: &
        advance_xp2_xpyp     ! Computes variance terms

    use surface_varnce_module, only:  &
        calc_surface_varnce ! Procedure

    use mixing_length, only: &
        compute_mixing_length ! Procedure

    use advance_windm_edsclrm_module, only:  &
        advance_windm_edsclrm  ! Procedure(s)

    use saturation, only:  &
        ! Procedure
        sat_mixrat_liq ! Saturation mixing ratio

    use advance_wp2_wp3_module, only:  &
        advance_wp2_wp3 ! Procedure

    use advance_xp3_module, only: &
        advance_xp3    ! Procedure(s)

    use clubb_precision, only:  &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use Skx_module, only: &
        Skx_func,           & ! Procedure(s)
        xp3_LG_2005_ansatz

    use clip_explicit, only: &
      clip_covars_denom ! Procedure(s)

    use T_in_K_module, only: &
        ! Read values from namelist
        thlm2T_in_K ! Procedure

    use sigma_sqd_w_module, only: &
        compute_sigma_sqd_w    ! Procedure(s)

    use stats_clubb_utilities, only: &
        stats_accumulate ! Procedure

    use stats_type_utilities, only:   &
        stat_update_var_pt,   & ! Procedure(s)
        stat_update_var,      &
        stat_begin_update,    &
        stat_begin_update_pt, &
        stat_end_update,      &
        stat_end_update_pt

    use stats_variables, only: &
        irtp2_bt,      & ! Variable(s)
        ithlp2_bt,     &
        irtpthlp_bt,   &
        iwp2_bt,       &
        iwp3_bt,       &
        ivp2_bt,       &
        iup2_bt,       &
        iwprtp_bt,     &
        iwpthlp_bt,    &
        iupwp_bt,      &
        ivpwp_bt,      &
        irtm_bt,       &
        ithlm_bt,      &
        ivm_bt,        &
        ium_bt,        &
        irvm,          &
        irel_humidity, &
        iwpthlp_zt   , &
        itau_zm_simp  

    use stats_variables, only: &
        iwprtp_zt,     &
        iup2_zt,       &
        ivp2_zt,       &
        iupwp_zt,      &
        ivpwp_zt,      &
        ithlp2_sf,     &
        irtp2_sf,      &
        irtpthlp_sf,   &
        iup2_sf,       &
        ivp2_sf,       &
        iwp2_sf,       &
        l_stats_samp,  &
        l_stats,       &
        stats_zt,      &
        stats_zm,      &
        stats_sfc,     &
        irtm_spur_src, &
        ithlm_spur_src

    use stats_variables, only: &
      irfrzm, & ! Variable(s)
      istability_correction

    use stats_variables, only: &
      iLscale_pert_1, & ! Variable(s)
      iLscale_pert_2

    use fill_holes, only: &
      vertical_integral, & ! Procedure(s)
      fill_holes_vertical

    use advance_helper_module, only: &
      calc_stability_correction, & ! Procedure(s)
      compute_Cx_fnc_Richardson, &
      calc_brunt_vaisala_freq_sqd, &
      term_wp2_splat, term_wp3_splat

    use interpolation, only: &
      pvertinterp

    implicit none

    !!! External
    intrinsic :: sqrt, min, max, exp, mod, real

    ! Constant Parameters
    logical, parameter :: &
      l_avg_Lscale = .false.    ! Lscale is calculated in subroutine compute_mixing_length
    ! if l_avg_Lscale is true, compute_mixing_length is called two additional times with
    ! perturbed values of rtm and thlm.  An average value of Lscale
    ! from the three calls to compute_mixing_length is then calculated.
    ! This reduces temporal noise in RICO, BOMEX, LBA, and other cases.

    logical, parameter :: &
      l_iter_xp2_xpyp = .true. ! Set to true when rtp2/thlp2/rtpthlp, et cetera are prognostic

    real( kind = core_rknd ), parameter :: &
      tau_const = 1000._core_rknd

    !!! Input Variables
    logical, intent(in) ::  &
      l_implemented ! True if CLUBB is being run within a large-scale host model,
                    !   rather than a standalone single-column model.

    real( kind = core_rknd ), intent(in) ::  &
      dt  ! Current timestep duration    [s]

    real( kind = core_rknd ), intent(in) ::  &
      fcor,  &          ! Coriolis forcing             [s^-1]
      sfc_elevation     ! Elevation of ground level    [m above MSL]

    integer, intent(in) :: &
      hydromet_dim      ! Total number of hydrometeor species        [#]

    ! Input Variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz) ::  &
      thlm_forcing,    & ! liquid potential temp. forcing (thermodynamic levels)    [K/s]
      rtm_forcing,     & ! total water forcing (thermodynamic levels)        [(kg/kg)/s]
      um_forcing,      & ! eastward wind forcing (thermodynamic levels)     [m/s/s]
      vm_forcing,      & ! northward wind forcing (thermodynamic levels)     [m/s/s]
      wprtp_forcing,   & ! total water turbulent flux forcing (momentum levels)    [m*K/s^2]
      wpthlp_forcing,  & ! liq pot temp turb flux forcing (momentum levels)   [m*(kg/kg)/s^2]
      rtp2_forcing,    & ! total water variance forcing (momentum levels)    [(kg/kg)^2/s]
      thlp2_forcing,   & ! liq pot temp variance forcing (momentum levels)   [K^2/s]
      rtpthlp_forcing, & ! <r_t'th_l'> covariance forcing (momentum levels) [K*(kg/kg)/s]
      wm_zm,           & ! vertical mean wind component on momentum levels  [m/s]
      wm_zt,           & ! vertical mean wind component on thermo. levels   [m/s]
      p_in_Pa,         & ! Air pressure (thermodynamic levels)       [Pa]
      rho_zm,          & ! Air density on momentum levels            [kg/m^3]
      rho,             & ! Air density on thermodynamic levels       [kg/m^3]
      exner,           & ! Exner function (thermodynamic levels)     [-]
      rho_ds_zm,       & ! Dry, static density on momentum levels    [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels     [kg/m^3]
      invrs_rho_ds_zm, & ! Inverse dry, static density on momentum levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inverse dry, static density on thermo levs.  [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on momentum levs. [K]
      thv_ds_zt,       & ! Dry, base-state theta_v on thermo levs.  [K]
      rfrzm              ! Total ice-phase water mixing ratio        [kg/kg]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet           ! Array of hydrometeors                [units vary]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      radf          ! Buoyancy production at cloud top due to longwave radiative cooling [m^2/s^3]

#ifdef CLUBBND_CAM
    real( kind = core_rknd ), intent(in) :: &
      varmu
#endif

    real( kind = core_rknd ), dimension(gr%nz, hydromet_dim), intent(in) :: &
      wphydrometp, & ! Covariance of w and a hydrometeor      [(m/s) <hm units>]
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' > (hm = hydrometeor) [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and hm (on thermo levs.) [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and hm (on thermo levs.)      [K <hm units>]

    real( kind = core_rknd ), intent(in) ::  &
      wpthlp_sfc,   & ! w' theta_l' at surface   [(m K)/s]
      wprtp_sfc,    & ! w' r_t' at surface       [(kg m)/( kg s)]
      upwp_sfc,     & ! u'w' at surface          [m^2/s^2]
      vpwp_sfc        ! v'w' at surface          [m^2/s^2]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz,sclr_dim) :: &
      sclrm_forcing    ! Passive scalar forcing         [{units vary}/s]

    real( kind = core_rknd ), intent(in),  dimension(sclr_dim) ::  &
      wpsclrp_sfc      ! Passive scalar flux at surface         [{units vary} m/s]

    ! Eddy passive scalar variables
    real( kind = core_rknd ), intent(in), dimension(gr%nz,edsclr_dim) :: &
      edsclrm_forcing  ! Eddy-diffusion passive scalar forcing    [{units vary}/s]

    real( kind = core_rknd ), intent(in),  dimension(edsclr_dim) ::  &
      wpedsclrp_sfc    ! Eddy-diffusion passive scalar flux at surface    [{units vary} m/s]

    ! Host model horizontal grid spacing, if part of host model.
    real( kind = core_rknd ), intent(in) :: &
      host_dx,  & ! East-west horizontal grid spacing     [m]
      host_dy     ! North-south horizontal grid spacing   [m]

    !!! Input/Output Variables
    ! These are prognostic or are planned to be in the future
    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  &
      um,      & ! eastward grid-mean wind component (thermodynamic levels)   [m/s]
      upwp,    & ! u'w' (momentum levels)                         [m^2/s^2]
      vm,      & ! northward grid-mean wind component (thermodynamic levels)   [m/s]
      vpwp,    & ! v'w' (momentum levels)                         [m^2/s^2]
      up2,     & ! u'^2 (momentum levels)                         [m^2/s^2]
      vp2,     & ! v'^2 (momentum levels)                         [m^2/s^2]
      rtm,     & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp,   & ! w' r_t' (momentum levels)                      [(kg/kg) m/s]
      thlm,    & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp,  & ! w'th_l' (momentum levels)                      [(m/s) K]
      rtp2,    & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      rtp3,    & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp2,   & ! th_l'^2 (momentum levels)                      [K^2]
      thlp3,   & ! th_l'^3 (thermodynamic levels)                 [K^3]
      rtpthlp, & ! r_t'th_l' (momentum levels)                    [(kg/kg) K]
      wp2,     & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3        ! w'^3 (thermodynamic levels)                    [m^3/s^3]

    ! Passive scalar variables
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    real( kind = core_rknd ), intent(inout), dimension(gr%nz) ::  &
      rcm,        & ! cloud water mixing ratio, r_c (thermo. levels) [kg/kg]
      cloud_frac, & ! cloud fraction (thermodynamic levels)          [-]
      wpthvp,     & ! < w' th_v' > (momentum levels)                 [kg/kg K]
      wp2thvp,    & ! < w'^2 th_v' > (thermodynamic levels)          [m^2/s^2 K]
      rtpthvp,    & ! < r_t' th_v' > (momentum levels)               [kg/kg K]
      thlpthvp      ! < th_l' th_v' > (momentum levels)              [K^2]

    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) :: &
      sclrpthvp    ! < sclr' th_v' > (momentum levels)   [units vary]

    type(pdf_parameter), intent(inout) :: &
      pdf_params,    & ! Fortran structure of PDF parameters on thermodynamic levels    [units vary]
      pdf_params_zm    ! Fortran structure of PDF parameters on momentum levels        [units vary]

#ifdef GFDL
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,sclr_dim) :: &  ! h1g, 2010-06-16
      sclrm_trsport_only  ! Passive scalar concentration due to pure transport [{units vary}/s]
#endif

    ! Eddy passive scalar variable
    real( kind = core_rknd ), intent(inout), dimension(gr%nz,edsclr_dim) :: &
      edsclrm   ! Eddy passive scalar grid-mean (thermo. levels)   [units vary]

    ! Variables that need to be output for use in other parts of the CLUBB
    ! code, such as microphysics (rcm, pdf_params), forcings (rcm), and/or
    ! BUGSrad (cloud_cover).
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  &
      rcm_in_layer, & ! rcm within cloud layer                          [kg/kg]
      cloud_cover     ! cloud cover                                     [-]

    ! Variables that need to be output for use in host models
    real( kind = core_rknd ), intent(out), dimension(gr%nz) ::  &
      wprcp,            & ! w'r_c' (momentum levels)                  [(kg/kg) m/s]
      ice_supersat_frac   ! ice cloud fraction (thermodynamic levels) [-]

    real( kind = core_rknd ), dimension(gr%nz) ::  &
      uprcp,              & ! < u' r_c' >              [(m kg)/(s kg)]
      vprcp                 ! < v' r_c' >              [(m kg)/(s kg)]


#if defined(CLUBB_CAM) || defined(GFDL)
    real( kind = core_rknd ), intent(out), dimension(gr%nz) :: &
      khzt, &       ! eddy diffusivity on thermo levels
      khzm          ! eddy diffusivity on momentum levels
#endif

#ifdef CLUBB_CAM
    real( kind = core_rknd), intent(out), dimension(gr%nz) :: &
      qclvar, &     ! cloud water variance
      thlprcp_out   ! thl'rc'
#endif


#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), intent(inout), dimension(gr%nz, min(1,sclr_dim) , 2) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
! ---> h1g, 2012-06-14
    logical, intent(in)                 ::  do_liquid_only_in_clubb
! <--- h1g, 2012-06-14
#endif

    real( kind = core_rknd ), intent(in), pointer ::  &
      upwp_sfc_pert,     & ! pertubed u'w' at surface          [m^2/s^2]
      vpwp_sfc_pert        ! pertubed v'w' at surface          [m^2/s^2]
    real( kind = core_rknd ), intent(inout), dimension(:), pointer ::  &
      um_pert,      & ! pertubed eastward grid-mean wind component (thermodynamic levels)   [m/s]
      vm_pert,      & ! pertubed northward grid-mean wind component (thermodynamic levels)   [m/s]
      upwp_pert,    & ! pertubed u'w' (momentum levels)                         [m^2/s^2]
      vpwp_pert       ! pertubed v'w' (momentum levels)                         [m^2/s^2]

    ! Local Variables
    integer :: i, k

#ifdef CLUBB_CAM
    integer ::  ixind
#endif

    ! Eric Raut declared this variable solely for output to disk
    real( kind = core_rknd ), dimension(gr%nz) :: &
      rc_coef,    & ! Coefficient of X'r_c' in Eq. (34) (t-levs.)   [K/(kg/kg)]
      rc_coef_zm    ! Coefficient of X'r_c' in Eq. (34) on m-levs.  [K/(kg/kg)]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      Km_zm, & ! Eddy diffusivity for momentum on zm grid levels [m^2/s]
      Kmh_zm   ! Eddy diffusivity for thermodynamic variables [m^2/s]

    logical, parameter ::  &
      l_use_buoy_mod_Km_zm = .false. ! .true. if we use a buoyancy-modified expression for Km_zm

    real( kind = core_rknd ), dimension(gr%nz) :: &
      tau_factor, &         ! factor that includes tau_zm in expression for Km_zm [s]
      Km_zm_denom_term, &   ! term in denominator of Km_zm [-]
      Km_zm_numerator_term  ! term in numerator of Km_zm [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      gamma_Skw_fnc, & ! Gamma as a function of skewness          [-]
      sigma_sqd_w,   & ! PDF width parameter (momentum levels)    [-]
      sigma_sqd_w_zt, & ! PDF width parameter (thermodynamic levels)    [-]
      sqrt_em_zt,    & ! sqrt( em ) on zt levels; where em is TKE [m/s]
      Lscale_pert_1, Lscale_pert_2, & ! For avg. calculation of Lscale  [m]
      thlm_pert_1, thlm_pert_2, &     ! For avg. calculation of Lscale  [K]
      rtm_pert_1, rtm_pert_2,   &     ! For avg. calculation of Lscale  [kg/kg]
      thlm_pert_pos_rt, thlm_pert_neg_rt, &     ! For avg. calculation of Lscale  [K]
      rtm_pert_pos_rt, rtm_pert_neg_rt          ! For avg. calculation of Lscale  [kg/kg]
    !Lscale_weight Uncomment this if you need to use this vairable at some point.

    real( kind = core_rknd ), dimension(gr%nz) :: &
      w_1_zm,        & ! Mean w (1st PDF component)                   [m/s]
      w_2_zm,        & ! Mean w (2nd PDF component)                   [m/s]
      varnce_w_1_zm, & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w_2_zm, & ! Variance of w (2nd PDF component)            [m^2/s^2]
      mixt_frac_zm    ! Weight of 1st PDF component (Sk_w dependent) [-]

    integer :: &
      wprtp_cl_num,   & ! Instance of w'r_t' clipping (1st or 3rd).
      wpthlp_cl_num,  & ! Instance of w'th_l' clipping (1st or 3rd).
      wpsclrp_cl_num, & ! Instance of w'sclr' clipping (1st or 3rd).
      upwp_cl_num,    & ! Instance of u'w' clipping (1st or 2nd).
      vpwp_cl_num       ! Instance of v'w' clipping (1st or 2nd).

    real( kind = core_rknd ), dimension(gr%nz) :: &
      rcp2_zt,              & ! r_c'^2 (on thermo. grid)             [kg^2/kg^2]
      cloud_frac_zm,        & ! Cloud Fraction on momentum grid      [-]
      ice_supersat_frac_zm, & ! Ice Cloud Fraction on momentum grid  [-]
      rtm_zm,               & ! Total water mixing ratio             [kg/kg]
      thlm_zm,              & ! Liquid potential temperature         [kg/kg]
      rcm_zm,               & ! Liquid water mixing ratio on m-levs. [kg/kg]
      sign_rtpthlp,         & ! Sign of the covariance rtpthlp       [-]
      wpsclrp_zt,           & ! Scalar flux on thermo. levels        [un. vary]
      sclrp2_zt               ! Scalar variance on thermo.levels     [un. vary]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      sclrp3    ! <sclr'^3> (thermodynamic levels)    [un. vary]

    real( kind = core_rknd ) :: &
      rtm_integral_before, &
      rtm_integral_after, &
      rtm_integral_forcing, &
      rtm_flux_top, &
      rtm_flux_sfc, &
      rtm_spur_src, &
      thlm_integral_before, &
      thlm_integral_after, &
      thlm_integral_forcing, &
      thlm_flux_top, &
      thlm_flux_sfc, &
      thlm_spur_src, &
      mu_pert_1, mu_pert_2, & ! For l_avg_Lscale
      mu_pert_pos_rt, mu_pert_neg_rt ! For l_Lscale_plume_centered

    !The following variables are defined for use when l_use_ice_latent = .true.
    type(pdf_parameter) :: &
      pdf_params_frz

    type(implicit_coefs_terms), dimension(gr%nz) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    real( kind = core_rknd ), dimension(gr%nz)  :: &
      rtm_frz, &
      thlm_frz

    real( kind = core_rknd ) :: &
      thlm1000, &
      thlm700

    real( kind = core_rknd ), dimension(gr%nz) :: &
      rcm_supersat_adj, & ! Adjustment to rcm due to spurious supersaturation
      rel_humidity        ! Relative humidity after PDF closure [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
       stability_correction,   & ! Stability correction factor
       tau_N2_zm,              & ! Tau with a static stability correction applied to it [s]
       tau_C6_zm,              & ! Tau values used for the C6 (pr1) term in wpxp [s]
       tau_C1_zm,              & ! Tau values used for the C1 (dp1) term in wp2 [s]
       Cx_fnc_Richardson,      & ! Cx_fnc computed from Richardson_num          [-]
       brunt_vaisala_freq_sqd, & ! Buoyancy frequency squared, N^2              [s^-2]
       invrs_tau_zm,           & ! One divided by tau on zm levels              [s^-1]
       invrs_tau_N2_zm,        & ! One divided by tau, including stability effects [s^-1]
       ustar,                  &   ! Friction velocity  [m/s]
       invrs_tau_zm_simp,      & 
       tau_zm_simp,            &
       tau_zt_simp


    real( kind = core_rknd ), parameter :: &
       ufmin = 0.01_core_rknd,       & ! minimum value of friction velocity     [m/s]
       z_displace = 20.0_core_rknd   ! displacement of log law profile above ground   [m]

    real( kind = core_rknd ) :: Lscale_max

    real( kind = core_rknd ) :: newmu

    ! Flag to sample stats in a particular call to subroutine
    ! pdf_closure_driver.
    logical :: l_samp_stats_in_pdf_call

   real( kind = core_rknd ), dimension(gr%nz) :: &
     wp2_splat, &   ! Tendency of <w'2> due to eddies compressing  [m^2/s^3]
     wp3_splat      ! Tendency of <w'3> due to eddies compressing  [m^3/s^4]

    ! Variables associated with upgradient momentum contributions due to cumuli
    !real( kind = core_rknd ), dimension(gr%nz) :: &
    !  Km_Skw_factor ! Factor, with value < 1, that reduces eddy diffusivity,
    !                                          Km_zm, in skewed layers
    !real( kind = core_rknd ),parameter :: &
    !  Km_Skw_thresh = zero_threshold, &  ! Value of Skw at which Skw correction kicks in
    !  Km_Skw_factor_efold = 0.5_core_rknd, & ! E-folding rate of exponential Skw correction
    !  Km_Skw_factor_min   = 0.2_core_rknd    ! Minimum value of Km_Skw_factor

    !----- Begin Code -----

    ! Sanity checks
    if ( clubb_at_least_debug_level( 0 ) ) then

      if ( l_Lscale_plume_centered .and. .not. l_avg_Lscale ) then
        write(fstderr,*) "l_Lscale_plume_centered requires l_avg_Lscale"
        write(fstderr,*) "Fatal error in advance_clubb_core"
        return
      end if

      if ( l_damp_wp2_using_em .and. (C1 /= C14 .or. l_stability_correct_tau_zm) ) then
        write(fstderr,*) "l_damp_wp2_using_em requires C1=C14 and l_stability_correct_tau_zm = F"
        write(fstderr,*) "Fatal error in advance_clubb_core"
        return
      end if

    end if

    ! Determine the maximum allowable value for Lscale (in meters).
    call set_Lscale_max( l_implemented, host_dx, host_dy, & ! intent(in)
                         Lscale_max )                       ! intent(out)

    if ( l_stats .and. l_stats_samp ) then
      ! Spurious source will only be calculated if rtm_ma and thlm_ma are zero.
      ! Therefore, wm must be zero or l_implemented must be true.
       if ( l_implemented .or. ( all( abs(wm_zt) < eps ) .and. &
           all( abs(wm_zm) < eps ) ) ) then
        ! Get the vertical integral of rtm and thlm before this function begins
        ! so that spurious source can be calculated
        rtm_integral_before  &
        = vertical_integral( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                             rtm(2:gr%nz), gr%dzt(2:gr%nz) )

        thlm_integral_before  &
        = vertical_integral( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                             thlm(2:gr%nz), gr%dzt(2:gr%nz) )
      end if
    end if

    !----------------------------------------------------------------
    ! Test input variables
    !----------------------------------------------------------------
    if ( clubb_at_least_debug_level( 2 ) ) then
      call parameterization_check &
           ( thlm_forcing, rtm_forcing, um_forcing,                             & ! intent(in)
             vm_forcing, wm_zm, wm_zt, p_in_Pa,                                 & ! intent(in)
             rho_zm, rho, exner, rho_ds_zm,                                     & ! intent(in)
             rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,                       & ! intent(in)
             thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc,             & ! intent(in)
             vpwp_sfc, um, upwp, vm, vpwp, up2, vp2,                            & ! intent(in)
             rtm, wprtp, thlm, wpthlp, wp2, wp3,                                & ! intent(in)
             rtp2, thlp2, rtpthlp, rcm,                                         & ! intent(in)
             "beginning of ",                                                   & ! intent(in)
             wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp, sclrp2,                & ! intent(in)
             sclrprtp, sclrpthlp, sclrm_forcing, edsclrm, edsclrm_forcing )       ! intent(in)

        if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) return
        end if

    end if
    !-----------------------------------------------------------------------

    if ( l_stats_samp ) then
      call stat_update_var( irfrzm, rfrzm, & ! intent(in)
                            stats_zt ) ! intent(inout)
    end if

    ! Set up budget stats variables.
    if ( l_stats_samp ) then

       call stat_begin_update( iwp2_bt, wp2 / dt, & ! intent(in)
                               stats_zm )           ! intent(inout)
       call stat_begin_update( ivp2_bt, vp2 / dt, & ! intent(in)
                               stats_zm )           ! intent(inout)
       call stat_begin_update( iup2_bt, up2 / dt, & ! intent(in)
                               stats_zm )           ! intent(inout)
       call stat_begin_update( iwprtp_bt, wprtp / dt, & ! intent(in)
                               stats_zm )               ! intent(inout)
       call stat_begin_update( iwpthlp_bt, wpthlp / dt, & ! intent(in)
                               stats_zm )                 ! intent(inout)
       if ( l_predict_upwp_vpwp ) then
          call stat_begin_update( iupwp_bt, upwp / dt, & ! intent(in)
                                  stats_zm )             ! intent(inout)
          call stat_begin_update( ivpwp_bt, vpwp / dt, & ! intent(in)
                                  stats_zm )             ! intent(inout)
       endif ! l_predict_upwp_vpwp
       call stat_begin_update( irtp2_bt, rtp2 / dt, & ! intent(in)
                               stats_zm )             ! intent(inout)
       call stat_begin_update( ithlp2_bt, thlp2 / dt, & ! intent(in)
                               stats_zm )               ! intent(inout)
       call stat_begin_update( irtpthlp_bt, rtpthlp / dt, & ! intent(in)
                               stats_zm )                   ! intent(inout)

       call stat_begin_update( irtm_bt, rtm / dt, & ! intent(in)
                               stats_zt )           ! intent(inout)
       call stat_begin_update( ithlm_bt, thlm / dt, & ! intent(in)
                               stats_zt )             ! intent(inout)
       call stat_begin_update( ium_bt, um / dt, & ! intent(in)
                               stats_zt )         ! intent(inout)
       call stat_begin_update( ivm_bt, vm / dt, & ! intent(in)
                               stats_zt )         ! intent(inout)
       call stat_begin_update( iwp3_bt, wp3 / dt, & ! intent(in)
                               stats_zt )           ! intent(inout)

    endif

    ! SET SURFACE VALUES OF FLUXES (BROUGHT IN)
    ! We only do this for host models that do not apply the flux
    ! elsewhere in the code (e.g. WRF).  In other cases the _sfc variables will
    ! only be used to compute the variance at the surface. -dschanen 8 Sept 2009
    if ( .not. l_host_applies_sfc_fluxes ) then

      wpthlp(1) = wpthlp_sfc
      wprtp(1)  = wprtp_sfc
      upwp(1)   = upwp_sfc
      vpwp(1)   = vpwp_sfc
      if (associated(upwp_sfc_pert) .and. associated(vpwp_sfc_pert) .and. &
           associated(upwp_pert) .and. associated(vpwp_pert)) then
         upwp_pert(1) = upwp_sfc_pert
         vpwp_pert(1) = vpwp_sfc_pert
      end if

      ! Set fluxes for passive scalars (if enabled)
      if ( sclr_dim > 0 ) then
        wpsclrp(1,1:sclr_dim)   = wpsclrp_sfc(1:sclr_dim)
      end if

      if ( edsclr_dim > 0 ) then
        wpedsclrp(1,1:edsclr_dim) = wpedsclrp_sfc(1:edsclr_dim)
      end if

    else

      wpthlp(1) = 0.0_core_rknd
      wprtp(1)  = 0.0_core_rknd
      upwp(1)   = 0.0_core_rknd
      vpwp(1)   = 0.0_core_rknd

      ! Set fluxes for passive scalars (if enabled)
      if ( sclr_dim > 0 ) then
        wpsclrp(1,1:sclr_dim) = 0.0_core_rknd
      end if

      if ( edsclr_dim > 0 ) then
        wpedsclrp(1,1:edsclr_dim) = 0.0_core_rknd
      end if

    end if ! ~l_host_applies_sfc_fluxes

#ifdef CLUBBND_CAM
    newmu = varmu
#else
    newmu = mu
#endif

    if ( ipdf_call_placement == ipdf_pre_advance_fields &
         .or. ipdf_call_placement == ipdf_pre_post_advance_fields ) then

       ! Sample stats in this call to subroutine pdf_closure_driver for
       ! both of these options (ipdf_pre_advance_fields and
       ! ipdf_pre_post_advance_fields).
       if ( ipdf_call_placement == ipdf_pre_advance_fields ) then
          l_samp_stats_in_pdf_call = .true.
       elseif ( ipdf_call_placement == ipdf_pre_post_advance_fields ) then
          l_samp_stats_in_pdf_call = .true.
       endif

       !########################################################################
       !#######                     CALL CLUBB's PDF                     #######
       !#######   AND OUTPUT PDF PARAMETERS AND INTEGRATED QUANTITITES   #######
       !########################################################################
       !call pdf_closure_driver( dt, hydromet_dim, rtm, wprtp,  & ! Intent(in)
       call pdf_closure_driver( dt, hydromet_dim, wprtp,       & ! Intent(in)
                                thlm, wpthlp, rtp2, rtp3,      & ! Intent(in)
                                thlp2, thlp3, rtpthlp, wp2,    & ! Intent(in)
                                wp3, wm_zm, wm_zt,             & ! Intent(in)
                                um, up2, upwp,                 & ! Intent(in)
                                vm, vp2, vpwp,                 & ! Intent(in)
                                p_in_Pa, exner,                & ! Intent(in)
                                thv_ds_zm, thv_ds_zt,          & ! Intent(in)
                                rfrzm, hydromet, wphydrometp,  & ! Intent(in)
                                wp2hmp, rtphmp_zt, thlphmp_zt, & ! Intent(in)
                                sclrm, wpsclrp, sclrp2,        & ! Intent(in)
                                sclrprtp, sclrpthlp,           & ! Intent(in)
                                l_samp_stats_in_pdf_call,      & ! Intent(in)
                                rtm,                           & ! Intent(i/o)
#ifdef GFDL
                                RH_crit(k, : , :),             & ! Intent(i/o)
                                do_liquid_only_in_clubb,       & ! Intent(in)
#endif
                                rcm, cloud_frac,               & ! Intent(out)
                                ice_supersat_frac, wprcp,      & ! Intent(out)
                                sigma_sqd_w, wpthvp, wp2thvp,  & ! Intent(out)
                                rtpthvp, thlpthvp, rc_coef,    & ! Intent(out)
                                rcm_in_layer, cloud_cover,     & ! Intent(out)
                                rcp2_zt, thlprcp, rc_coef_zm,  & ! Intent(out)
                                rtm_frz, thlm_frz, sclrpthvp,  & ! Intent(out)
                                wp4, wp2rtp, wprtp2, wp2thlp,  & ! Intent(out)
                                wpthlp2, wprtpthlp, wp2rcp,    & ! Intent(out)
                                rtprcp, rcp2,                  & ! Intent(out)
                                uprcp, vprcp,                  & ! Intent(out)
                                Skw_velocity,                  & ! Intent(out)
                                cloud_frac_zm,                 & ! Intent(out)
                                ice_supersat_frac_zm,          & ! Intent(out)
                                rtm_zm, thlm_zm, rcm_zm,       & ! Intent(out)
                                rcm_supersat_adj,              & ! Intent(out)
                                wp2sclrp, wpsclrp2, sclrprcp,  & ! Intent(out)
                                wpsclrprtp, wpsclrpthlp,       & ! Intent(out)
                                pdf_params, pdf_params_frz,    & ! Intent(out)
                                pdf_params_zm,                 & ! Intent(out)
                                pdf_implicit_coefs_terms )       ! Intent(out)

    endif ! ipdf_call_placement == ipdf_pre_advance_fields
          ! or ipdf_call_placement == ipdf_pre_post_advance_fields

    ! Interpolate wp3 to momentum levels, and wp2 to thermodynamic levels
    ! and then compute Skw for m & t grid.
    wp2_zt = max( zm2zt( wp2 ), w_tol_sqd ) ! Positive definite quantity
    wp3_zm = zt2zm( wp3 )

    Skw_zt(1:gr%nz) = Skx_func( wp2_zt(1:gr%nz), wp3(1:gr%nz), w_tol )
    Skw_zm(1:gr%nz) = Skx_func( wp2(1:gr%nz), wp3_zm(1:gr%nz), w_tol )

    if ( ipdf_call_placement == ipdf_post_advance_fields ) then

       ! Calculate sigma_sqd_w here in order to avoid having to pass it in
       ! and out of subroutine advance_clubb_core.
       if ( l_gamma_Skw .and. &
            abs(gamma_coef-gamma_coefb) > abs(gamma_coef+gamma_coefb)*eps/2) then

          gamma_Skw_fnc = gamma_coefb + (gamma_coef-gamma_coefb) &
                *exp( -(1.0_core_rknd/2.0_core_rknd) * (Skw_zm/gamma_coefc)**2 )
       else
          gamma_Skw_fnc = gamma_coef
       endif

       ! Compute sigma_sqd_w (dimensionless PDF width parameter)
       sigma_sqd_w &
       = compute_sigma_sqd_w( gamma_Skw_fnc, wp2, thlp2, rtp2, &
                              up2, vp2, wpthlp, wprtp, upwp, vpwp )

       ! Smooth in the vertical using interpolation
       sigma_sqd_w = zt2zm( zm2zt( sigma_sqd_w ) )
       sigma_sqd_w = max( zero_threshold, sigma_sqd_w ) ! Pos. def. quantity

    endif ! ipdf_call_placement == ipdf_post_advance_fields

    ! Compute the a3 coefficient (formula 25 in `Equations for CLUBB')
    ! Note:  a3 has been modified because the wp3 turbulent advection term is
    !        now discretized on its own.  This removes the "- 3" from the end.
!   a3_coef = 3.0_core_rknd * sigma_sqd_w*sigma_sqd_w  &
!      + 6.0_core_rknd*(1.0_core_rknd-sigma_sqd_w)*sigma_sqd_w  &
!      + (1.0_core_rknd-sigma_sqd_w)*(1.0_core_rknd-sigma_sqd_w)

    ! This is a simplified version of the formula above.
    ! Note:  a3 has been modified because the wp3 turbulent advection term is
    !        now discretized on its own.
    a3_coef = -2._core_rknd * ( 1._core_rknd - sigma_sqd_w )**2 + 3.0_core_rknd

    ! We found we obtain fewer spikes in wp3 when we clip a3 to be no greater
    ! than -1.4 -dschanen 4 Jan 2011
    !a3_coef = max( a3_coef, -1.4_core_rknd ) ! Known magic number
    a3_coef = max( a3_coef, 1.6_core_rknd ) ! Known magic number

    a3_coef_zt = zm2zt( a3_coef )

    ! Interpolate thlp2, rtp2, and rtpthlp to thermodynamic levels.
    thlp2_zt   = max( zm2zt( thlp2 ), thl_tol**2 ) ! Positive def. quantity
    rtp2_zt    = max( zm2zt( rtp2 ), rt_tol**2 )   ! Positive def. quantity
    rtpthlp_zt = zm2zt( rtpthlp )

    ! Compute wp3 / wp2 on zt levels.  Always use the interpolated value in the
    ! denominator since it's less likely to create spikes
    wp3_on_wp2_zt = ( wp3(1:gr%nz) / max( wp2_zt(1:gr%nz), w_tol_sqd ) )

    ! Clip wp3_on_wp2_zt if it's too large
    do k=1, gr%nz
      if( wp3_on_wp2_zt(k) < 0._core_rknd ) then
        wp3_on_wp2_zt = max( -1000._core_rknd, wp3_on_wp2_zt )
      else
        wp3_on_wp2_zt = min( 1000._core_rknd, wp3_on_wp2_zt )
      end if
    end do

    ! Compute wp3_on_wp2 by interpolating wp3_on_wp2_zt
    wp3_on_wp2 = zt2zm( wp3_on_wp2_zt )

    ! Smooth again as above
    wp3_on_wp2_zt = zm2zt( wp3_on_wp2 )

      !----------------------------------------------------------------
      ! Compute thvm
      !----------------------------------------------------------------

      thvm = thlm + ep1 * thv_ds_zt * rtm &
                  + ( Lv/(Cp*exner) - ep2 * thv_ds_zt ) * rcm

      !----------------------------------------------------------------
      ! Compute tke (turbulent kinetic energy)
      !----------------------------------------------------------------

      if ( .not. l_tke_aniso ) then
        ! tke is assumed to be 3/2 of wp2
        em = three_halves * wp2
      else
        em = 0.5_core_rknd * ( wp2 + vp2 + up2 )
      end if

      sqrt_em_zt = SQRT( MAX( em_min, zm2zt( em ) ) )

      !----------------------------------------------------------------
      ! Compute mixing length
      !----------------------------------------------------------------

      if ( .not. l_diag_Lscale_from_tau ) then ! compute Lscale 1st, using buoyant parcel calc

      if ( l_avg_Lscale .and. .not. l_Lscale_plume_centered ) then

         ! Call compute length two additional times with perturbed values
         ! of rtm and thlm so that an average value of Lscale may be calculated.

         do k = 1, gr%nz, 1
            sign_rtpthlp(k) = sign( one, rtpthlp(k) )
         enddo

         if ( l_use_ice_latent ) then

            ! Include the effects of ice in the length scale calculation

            rtm_pert_1  = rtm_frz &
                          + Lscale_pert_coef * sqrt( max( rtp2, rt_tol**2 ) )
            thlm_pert_1 = thlm_frz &
                          + sign_rtpthlp * Lscale_pert_coef &
                            * sqrt( max( thlp2, thl_tol**2 ) )
            mu_pert_1   = newmu / Lscale_mu_coef

            rtm_pert_2  = rtm_frz &
                          - Lscale_pert_coef * sqrt( max( rtp2, rt_tol**2 ) )
            thlm_pert_2 = thlm_frz &
                          - sign_rtpthlp * Lscale_pert_coef &
                            * sqrt( max( thlp2, thl_tol**2 ) )
            mu_pert_2   = newmu * Lscale_mu_coef

         else

            rtm_pert_1  = rtm &
                          + Lscale_pert_coef * sqrt( max( rtp2, rt_tol**2 ) )
            thlm_pert_1 = thlm &
                          + sign_rtpthlp * Lscale_pert_coef &
                            * sqrt( max( thlp2, thl_tol**2 ) )
            mu_pert_1   = newmu / Lscale_mu_coef

            rtm_pert_2  = rtm &
                          - Lscale_pert_coef * sqrt( max( rtp2, rt_tol**2 ) )
            thlm_pert_2 = thlm &
                          - sign_rtpthlp * Lscale_pert_coef &
                            * sqrt( max( thlp2, thl_tol**2 ) )
            mu_pert_2   = newmu * Lscale_mu_coef

         endif

         call compute_mixing_length( thvm, thlm_pert_1,                            & !intent(in)
                              rtm_pert_1, em, Lscale_max, p_in_Pa,                 & !intent(in)
                              exner, thv_ds_zt, mu_pert_1, l_implemented,          & !intent(in)
                              Lscale_pert_1, Lscale_up, Lscale_down )                !intent(out)

         call compute_mixing_length( thvm, thlm_pert_2,                            & !intent(in)
                              rtm_pert_2, em, Lscale_max, p_in_Pa,                 & !intent(in)
                              exner, thv_ds_zt, mu_pert_2, l_implemented,          & !intent(in)
                              Lscale_pert_2, Lscale_up, Lscale_down )                !intent(out)

      else if ( l_avg_Lscale .and. l_Lscale_plume_centered ) then
        ! Take the values of thl and rt based one 1st or 2nd plume

        do k = 1, gr%nz, 1
          sign_rtpthlp(k) = sign(1.0_core_rknd, rtpthlp(k))
        end do

        if ( l_use_ice_latent ) then
          where ( pdf_params_frz%rt_1 > pdf_params_frz%rt_2 )
            rtm_pert_pos_rt = pdf_params_frz%rt_1 &
                       + Lscale_pert_coef * sqrt( max( pdf_params_frz%varnce_rt_1, rt_tol**2 ) )
            thlm_pert_pos_rt = pdf_params_frz%thl_1 + ( sign_rtpthlp * Lscale_pert_coef &
                       * sqrt( max( pdf_params_frz%varnce_thl_1, thl_tol**2 ) ) )
            thlm_pert_neg_rt = pdf_params_frz%thl_2 - ( sign_rtpthlp * Lscale_pert_coef &
                       * sqrt( max( pdf_params_frz%varnce_thl_2, thl_tol**2 ) ) )
            rtm_pert_neg_rt = pdf_params_frz%rt_2 &
                       - Lscale_pert_coef * sqrt( max( pdf_params_frz%varnce_rt_2, rt_tol**2 ) )
            !Lscale_weight = pdf_params%mixt_frac
          elsewhere
            rtm_pert_pos_rt = pdf_params_frz%rt_2 &
                       + Lscale_pert_coef * sqrt( max( pdf_params_frz%varnce_rt_2, rt_tol**2 ) )
            thlm_pert_pos_rt = pdf_params_frz%thl_2 + ( sign_rtpthlp * Lscale_pert_coef &
                       * sqrt( max( pdf_params_frz%varnce_thl_2, thl_tol**2 ) ) )
            thlm_pert_neg_rt = pdf_params_frz%thl_1 - ( sign_rtpthlp * Lscale_pert_coef &
                       * sqrt( max( pdf_params_frz%varnce_thl_1, thl_tol**2 ) ) )
            rtm_pert_neg_rt = pdf_params_frz%rt_1 &
                       - Lscale_pert_coef * sqrt( max( pdf_params_frz%varnce_rt_1, rt_tol**2 ) )
            !Lscale_weight = 1.0_core_rknd - pdf_params%mixt_frac
          end where
        else
          where ( pdf_params%rt_1 > pdf_params%rt_2 )
            rtm_pert_pos_rt = pdf_params%rt_1 &
                       + Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_1, rt_tol**2 ) )
            thlm_pert_pos_rt = pdf_params%thl_1 + ( sign_rtpthlp * Lscale_pert_coef &
                       * sqrt( max( pdf_params%varnce_thl_1, thl_tol**2 ) ) )
            thlm_pert_neg_rt = pdf_params%thl_2 - ( sign_rtpthlp * Lscale_pert_coef &
                       * sqrt( max( pdf_params%varnce_thl_2, thl_tol**2 ) ) )
            rtm_pert_neg_rt = pdf_params%rt_2 &
                       - Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_2, rt_tol**2 ) )
            !Lscale_weight = pdf_params%mixt_frac
          elsewhere
            rtm_pert_pos_rt = pdf_params%rt_2 &
                       + Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_2, rt_tol**2 ) )
            thlm_pert_pos_rt = pdf_params%thl_2 + ( sign_rtpthlp * Lscale_pert_coef &
                       * sqrt( max( pdf_params%varnce_thl_2, thl_tol**2 ) ) )
            thlm_pert_neg_rt = pdf_params%thl_1 - ( sign_rtpthlp * Lscale_pert_coef &
                       * sqrt( max( pdf_params%varnce_thl_1, thl_tol**2 ) ) )
            rtm_pert_neg_rt = pdf_params%rt_1 &
                       - Lscale_pert_coef * sqrt( max( pdf_params%varnce_rt_1, rt_tol**2 ) )
            !Lscale_weight = 1.0_core_rknd - pdf_params%mixt_frac
          end where
        end if
        mu_pert_pos_rt  = newmu / Lscale_mu_coef
        mu_pert_neg_rt  = newmu * Lscale_mu_coef

        ! Call length with perturbed values of thl and rt
        call compute_mixing_length( thvm, thlm_pert_pos_rt,                    & !intent(in)
                           rtm_pert_pos_rt, em, Lscale_max, p_in_Pa,           & !intent(in)
                           exner, thv_ds_zt, mu_pert_pos_rt, l_implemented,    & !intent(in)
                           Lscale_pert_1, Lscale_up, Lscale_down )               !intent(out)

        call compute_mixing_length( thvm, thlm_pert_neg_rt,                    & !intent(in)
                           rtm_pert_neg_rt, em, Lscale_max, p_in_Pa,           & !intent(in)
                           exner, thv_ds_zt, mu_pert_neg_rt, l_implemented,    & !intent(in)
                           Lscale_pert_2, Lscale_up, Lscale_down )               !intent(out)
      else
        Lscale_pert_1 = unused_var ! Undefined
        Lscale_pert_2 = unused_var ! Undefined

      end if ! l_avg_Lscale

      if ( l_stats_samp ) then
        call stat_update_var( iLscale_pert_1, Lscale_pert_1, & ! intent(in)
                              stats_zt )                             ! intent(inout)
        call stat_update_var( iLscale_pert_2, Lscale_pert_2, & ! intent(in)
                              stats_zt )                             ! intent(inout)
      end if ! l_stats_samp

      ! ********** NOTE: **********
      ! This call to compute_mixing_length must be last.  Otherwise, the values of
      ! Lscale_up and Lscale_down in stats will be based on perturbation length scales
      ! rather than the mean length scale.

      ! Diagnose CLUBB's turbulent mixing length scale.
      call compute_mixing_length( thvm, thlm,                           & !intent(in)
                           rtm, em, Lscale_max, p_in_Pa,                & !intent(in)
                           exner, thv_ds_zt, newmu, l_implemented,      & !intent(in)
                           Lscale, Lscale_up, Lscale_down )               !intent(out)

      if ( l_avg_Lscale ) then
        if ( l_Lscale_plume_centered ) then
          ! Weighted average of mean, pert_1, & pert_2
!       Lscale = 0.5_core_rknd * ( Lscale + Lscale_weight*Lscale_pert_1 &
!                                  + (1.0_core_rknd-Lscale_weight)*Lscale_pert_2 )

          ! Weighted average of just the perturbed values
!       Lscale = Lscale_weight*Lscale_pert_1 + (1.0_core_rknd-Lscale_weight)*Lscale_pert_2

          ! Un-weighted average of just the perturbed values
          Lscale = 0.5_core_rknd*( Lscale_pert_1 + Lscale_pert_2 )
        else
          Lscale = (1.0_core_rknd/3.0_core_rknd) * ( Lscale + Lscale_pert_1 + Lscale_pert_2 )
        end if
      end if

      !----------------------------------------------------------------
      ! Dissipation time
      !----------------------------------------------------------------

      ! Calculate CLUBB's turbulent eddy-turnover time scale as
      !   CLUBB's length scale divided by a velocity scale.
      tau_zt = MIN( Lscale / sqrt_em_zt, taumax )
      tau_zm = MIN( ( MAX( zt2zm( Lscale ), zero_threshold )  &
                     / SQRT( MAX( em_min, em ) ) ), taumax )
! End Vince Larson's replacement.


      else ! l_diag_Lscale_from_tau = .true., diagnose simple tau and Lscale.

        call calc_brunt_vaisala_freq_sqd( thlm, exner, rtm, rcm, p_in_Pa, thvm, & ! intent(in)
                                          brunt_vaisala_freq_sqd )   ! intent(out)

        ustar = max( ( upwp_sfc**2 + vpwp_sfc**2 )**(one_fourth), ufmin )

        invrs_tau_zm = C_invrs_tau_bkgnd / tau_const &
         + C_invrs_tau_sfc * ( ustar / vonk ) / ( gr%zm - sfc_elevation + z_displace ) &
         + C_invrs_tau_shear * zt2zm( zm2zt( sqrt( (ddzt( um ))**2 + (ddzt( vm ))**2 ) ) )  &
         + C_invrs_tau_N2 * sqrt( max( zero_threshold, &
              zt2zm( zm2zt( brunt_vaisala_freq_sqd ) ) - 1e-4_core_rknd) )

        invrs_tau_zm_simp = C_invrs_tau_bkgnd  / tau_const & 
         + C_invrs_tau_sfc * ( ustar / vonk ) / ( gr%zm - sfc_elevation + z_displace ) & 
         + C_invrs_tau_N2 * zt2zm( zm2zt( sqrt( (ddzt( um ))**2 + (ddzt( vm ))**2 ) ) )

        if ( gr%zm(1) - sfc_elevation + z_displace < eps ) then
             stop  "Lowest zm grid level is below ground in CLUBB."
        end if

        tau_zm = one / invrs_tau_zm
        tau_zm_simp = one / invrs_tau_zm_simp  

        tau_zt = zm2zt( tau_zm )

        invrs_tau_N2_zm = invrs_tau_zm  &
                          + C_invrs_tau_N2 * sqrt( max( zero_threshold, brunt_vaisala_freq_sqd ) )

        tau_N2_zm = tau_zm

!        tau_zt = tau_const / &
!                     ( one + 0.1_core_rknd * tau_const * &
!                             sqrt( max( zero_threshold, brunt_vaisala_freq_sqd ) ) )
!        tau_zm = max( zero_threshold, zt2zm( tau_zt ) )

        Lscale = tau_zt * sqrt_em_zt

      end if ! l_diag_Lscale_from_tau

      ! Modification to damp noise in stable region
! Vince Larson commented out because it may prevent turbulence from
!    initiating in unstable regions.  7 Jul 2007
!       do k = 1, gr%nz
!         if ( wp2(k) <= 0.005_core_rknd ) then
!           tau_zt(k) = taumin
!           tau_zm(k) = taumin
!         end if
!       end do
! End Vince Larson's commenting.

      !----------------------------------------------------------------
      ! Eddy diffusivity coefficient
      !----------------------------------------------------------------
      ! c_K is 0.548 usually (Duynkerke and Driedonks 1987)
      ! CLUBB uses a smaller value to better fit empirical data.

      ! Calculate CLUBB's eddy diffusivity as
      !   CLUBB's length scale times a velocity scale.
      Kh_zt = c_K * Lscale * sqrt_em_zt
      Kh_zm = c_K * max( zt2zm( Lscale ), zero_threshold )  &
                  * sqrt( max( em, em_min ) )

#if defined(CLUBB_CAM) || defined(GFDL)
      khzt(:) = Kh_zt(:)
      khzm(:) = Kh_zm(:)
#endif

      ! Vertical compression of eddies causes gustiness (increase in up2 and vp2)
      call term_wp2_splat( C_wp2_splat, gr%nz, dt, wp2, wp2_zt, tau_zm, & ! Intent(in)
                           wp2_splat )                                ! Intent(out)
      ! Vertical compression of eddies also diminishes w'3
      call term_wp3_splat( C_wp2_splat, gr%nz, dt, wp2, wp3, tau_zt, & ! Intent(in)
                           wp3_splat )                             ! Intent(out)

      !----------------------------------------------------------------
      ! Set Surface variances
      !----------------------------------------------------------------

      ! Surface variances should be set here, before the call to either
      ! advance_xp2_xpyp or advance_wp2_wp3.
      ! Surface effects should not be included with any case where the lowest
      ! level is not the ground level.  Brian Griffin.  December 22, 2005.
      if ( abs(gr%zm(1)-sfc_elevation) <= abs(gr%zm(1)+sfc_elevation)*eps/2) then

        ! Reflect surface varnce changes in budget
        if ( l_stats_samp ) then
          call stat_begin_update_pt( ithlp2_sf, 1,      &      ! intent(in)
           thlp2(1) / dt,    &                                 ! intent(in)
                                     stats_zm )                      ! intent(inout)
          call stat_begin_update_pt( irtp2_sf, 1,       &      ! intent(in)
            rtp2(1) / dt,    &                                 ! intent(in)
                                     stats_zm )                      ! intent(inout)
          call stat_begin_update_pt( irtpthlp_sf, 1,    &      ! intent(in)
            rtpthlp(1) / dt, &                                 ! intent(in)
                                     stats_zm )                      ! intent(inout)
          call stat_begin_update_pt( iup2_sf, 1,        &      ! intent(in)
            up2(1) / dt,     &                                 ! intent(in)
                                     stats_zm )                      ! intent(inout)
          call stat_begin_update_pt( ivp2_sf, 1,        &      ! intent(in)
            vp2(1) / dt,     &                                 ! intent(in)
                                     stats_zm )                      ! intent(inout)
          call stat_begin_update_pt( iwp2_sf, 1,        &      ! intent(in)
            wp2(1) / dt,     &                                 ! intent(in)
                                     stats_zm )                      ! intent(inout)
        end if

        ! Diagnose surface variances based on surface fluxes.
        call calc_surface_varnce( upwp_sfc, vpwp_sfc, wpthlp_sfc, wprtp_sfc, &      ! intent(in)
                             um(2), vm(2), Lscale_up(2), wpsclrp_sfc,        &      ! intent(in)
                             wp2_splat(1), tau_zm(1),                        &      ! intent(in)
                             wp2(1), up2(1), vp2(1),                         &      ! intent(out)
                             thlp2(1), rtp2(1), rtpthlp(1),                  &      ! intent(out)
                             sclrp2(1,1:sclr_dim),                           &      ! intent(out)
                             sclrprtp(1,1:sclr_dim),                         &      ! intent(out)
                             sclrpthlp(1,1:sclr_dim) )                              ! intent(out)

        if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Error calling calc_surface_varnce"
            return
          end if
        end if

        ! Update surface stats
        if ( l_stats_samp ) then
          call stat_end_update_pt( ithlp2_sf, 1, &                ! intent(in)
            thlp2(1) / dt, &                                      ! intent(in)
                                   stats_zm )                           ! intent(inout)
          call stat_end_update_pt( irtp2_sf, 1, &                 ! intent(in)
            rtp2(1) / dt, &                                       ! intent(in)
                                   stats_zm )                           ! intent(inout)
          call stat_end_update_pt( irtpthlp_sf, 1, &              ! intent(in)
            rtpthlp(1) / dt, &                                    ! intent(in)
                                   stats_zm )                           ! intent(inout)
          call stat_end_update_pt( iup2_sf, 1, &                  ! intent(in)
            up2(1) / dt, &                                        ! intent(in)
                                   stats_zm )                           ! intent(inout)
          call stat_end_update_pt( ivp2_sf, 1, &                  ! intent(in)
            vp2(1) / dt, &                                        ! intent(in)
                                   stats_zm )                           ! intent(inout)
          call stat_end_update_pt( iwp2_sf, 1, &                  ! intent(in)
            wp2(1) / dt, &                                        ! intent(in)
                                   stats_zm )                           ! intent(inout)
        end if

      else

        ! Variances for cases where the lowest level is not at the surface.
        ! Eliminate surface effects on lowest level variances.
        wp2(1)     = w_tol_sqd
        up2(1)     = w_tol_sqd
        vp2(1)     = w_tol_sqd
        thlp2(1)   = thl_tol**2
        rtp2(1)    = rt_tol**2
        rtpthlp(1) = 0.0_core_rknd

        do i = 1, sclr_dim, 1
          sclrp2(1,i)    = 0.0_core_rknd
          sclrprtp(1,i)  = 0.0_core_rknd
          sclrpthlp(1,i) = 0.0_core_rknd
        end do

      end if ! gr%zm(1) == sfc_elevation


      !#######################################################################
      !############## ADVANCE PROGNOSTIC VARIABLES ONE TIMESTEP ##############
      !#######################################################################

      if ( l_stats_samp ) then
        call stat_update_var( irvm, rtm - rcm, & !intent(in)
                              stats_zt )               !intent(inout)

        ! Output relative humidity (q/q∗ where q∗ is the saturation mixing ratio over liquid)
        ! Added an extra check for irel_humidity > 0; otherwise, if both irsat = 0 and
        ! irel_humidity = 0, rsat is not computed, leading to a floating-point exception
        ! when stat_update_var is called for rel_humidity.  ldgrant
        if ( irel_humidity > 0 ) then
          ! Recompute rsat and rel_humidity. They might have changed.
          rsat = sat_mixrat_liq( p_in_Pa, thlm2T_in_K( thlm, exner, rcm ) )
          rel_humidity = (rtm - rcm) / rsat

          call stat_update_var( irel_humidity, rel_humidity, &             ! intent(in)
                                stats_zt)                                  ! intent(inout)
        end if ! irel_humidity > 0
      end if ! l_stats_samp

      !----------------------------------------------------------------
      ! Advance rtm/wprtp and thlm/wpthlp one time step
      !----------------------------------------------------------------
      if ( l_call_pdf_closure_twice ) then
        w_1_zm        = pdf_params_zm%w_1
        w_2_zm        = pdf_params_zm%w_2
        varnce_w_1_zm = pdf_params_zm%varnce_w_1
        varnce_w_2_zm = pdf_params_zm%varnce_w_2
        mixt_frac_zm = pdf_params_zm%mixt_frac
      else
        w_1_zm        = zt2zm( pdf_params%w_1 )
        w_2_zm        = zt2zm( pdf_params%w_2 )
        varnce_w_1_zm = zt2zm( pdf_params%varnce_w_1 )
        varnce_w_2_zm = zt2zm( pdf_params%varnce_w_2 )
        mixt_frac_zm  = zt2zm( pdf_params%mixt_frac )
      end if

      ! Determine stability correction factor
      stability_correction = calc_stability_correction( thlm, Lscale, em, exner, rtm, rcm, & ! In
                                                        p_in_Pa, thvm ) ! In
      if ( l_stats_samp ) then
        call stat_update_var( istability_correction, stability_correction, & ! In
                              stats_zm ) ! In/Out
      end if

      ! Here we determine if we're using tau_zm or tau_N2_zm, which is tau
      ! that has been stability corrected for stably stratified regions.
      ! -dschanen 7 Nov 2014
      if ( l_stability_correct_tau_zm ) then
        ! Determine the static stability corrected version of tau_zm
        ! Create a damping time scale that is more strongly damped at the
        ! altitudes where the Brunt-Vaisala frequency (N^2) is large.
        tau_N2_zm = tau_zm / stability_correction
        tau_C6_zm = tau_N2_zm
        tau_C1_zm = tau_N2_zm

      else
        tau_N2_zm = unused_var
        tau_C6_zm = tau_zm
        tau_C1_zm = tau_zm

      end if ! l_stability_correction

      ! Cx_fnc_Richardson is only used if one of these flags is true,
      ! otherwise its value is irrelevant, set it to 0 to avoid NaN problems
      if ( l_use_C7_Richardson .or. l_use_C11_Richardson .or. l_use_wp3_pr3 ) then
          call compute_Cx_Fnc_Richardson( thlm, um, vm, em, Lscale, exner, rtm, &
                                          rcm, p_in_Pa, thvm, rho_ds_zm, &
                                          Cx_fnc_Richardson )
      else
          Cx_fnc_Richardson = 0.0
      end if

      ! Advance the prognostic equations for
      !   the scalar grid means (rtm, thlm, sclrm) and
      !   scalar turbulent fluxes (wprtp, wpthlp, and wpsclrp)
      !   by one time step.
      ! advance_xm_wpxp_bad_wp2 ! Test error comment, DO NOT modify or move
      call advance_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2,              & ! intent(in)
                            Lscale, wp3_on_wp2, wp3_on_wp2_zt, Kh_zt, Kh_zm, & ! intent(in)
                            tau_C6_zm, Skw_zm, wp2rtp, rtpthvp, rtm_forcing, & ! intent(in)
                            wprtp_forcing, rtm_ref, wp2thlp, thlpthvp,       & ! intent(in)
                            thlm_forcing, wpthlp_forcing, thlm_ref,          & ! intent(in)
                            rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm,           & ! intent(in)
                            invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2,         & ! intent(in)
                            w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm,    & ! intent(in)
                            mixt_frac_zm, l_implemented, em, wp2sclrp,       & ! intent(in)
                            sclrpthvp, sclrm_forcing, sclrp2, exner, rcm,    & ! intent(in)
                            p_in_Pa, thvm, Cx_fnc_Richardson,                & ! intent(in)
                            pdf_implicit_coefs_terms,                        & ! intent(in)
                            um_forcing, vm_forcing, ug, vg, wpthvp,          & ! intent(in)
                            fcor, um_ref, vm_ref, up2, vp2,                  & ! intent(in)
                            uprcp, vprcp, rc_coef,                           & ! intent(in)
                            rtm, wprtp, thlm, wpthlp,                        & ! intent(inout)
                            sclrm, wpsclrp, um, upwp, vm, vpwp,              & ! intent(inout)
                            um_pert, vm_pert, upwp_pert, vpwp_pert)            ! intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Error calling advance_xm_wpxp"
            return
          end if
      end if

      ! Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
      ! This code won't work unless rtm >= 0 !!!
      ! We do not clip rcm_in_layer because rcm_in_layer only influences
      ! radiation, and we do not want to bother recomputing it.  6 Aug 2009
      call clip_rcm( rtm, 'rtm < rcm in advance_xm_wpxp',             & ! intent(in)
                     rcm )                                              ! intent(inout)

#ifdef GFDL
      call advance_sclrm_Nd_diffusion_OG( dt, &  ! h1g, 2012-06-16     ! intent(in)
                                          sclrm, sclrm_trsport_only, & ! intent(inout)
                                          Kh_zm,  cloud_frac )         ! intent(in)
#endif

      !----------------------------------------------------------------
      ! Compute some of the variances and covariances.  These include the variance of
      ! total water (rtp2), liquid water potential temperature (thlp2), their
      ! covariance (rtpthlp), and the variance of horizontal wind (up2 and vp2).
      ! The variance of vertical velocity is computed later.
      !----------------------------------------------------------------

      ! We found that certain cases require a time tendency to run
      ! at shorter timesteps so these are prognosed now.

      ! We found that if we call advance_xp2_xpyp first, we can use a longer timestep.

      ! Advance the prognostic equations
      !   for scalar variances and covariances,
      !   plus the horizontal wind variances by one time step, by one time step.
      call advance_xp2_xpyp( tau_zm, wm_zm, rtm, wprtp, thlm,        & ! intent(in)
                             wpthlp, wpthvp, um, vm, wp2, wp2_zt,    & ! intent(in)
                             wp3, upwp, vpwp, sigma_sqd_w, Skw_zm,   & ! intent(in)
                             wprtp2, wpthlp2, wprtpthlp,             & ! intent(in)
                             Kh_zt, rtp2_forcing, thlp2_forcing,     & ! intent(in)
                             rtpthlp_forcing, rho_ds_zm, rho_ds_zt,  & ! intent(in)
                             invrs_rho_ds_zm, thv_ds_zm, cloud_frac, & ! intent(in)
                             Lscale, wp3_on_wp2, wp3_on_wp2_zt,      & ! intent(in)
                             pdf_implicit_coefs_terms,               & ! intent(in)
                             l_iter_xp2_xpyp, dt,                    & ! intent(in)
                             sclrm, wpsclrp,                         & ! intent(in)
                             wpsclrp2, wpsclrprtp, wpsclrpthlp,      & ! intent(in)
                             wp2_splat,                              & ! intent(in)
                             rtp2, thlp2, rtpthlp, up2, vp2,         & ! intent(inout)
                             sclrp2, sclrprtp, sclrpthlp)              ! intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Error calling advance_xp2_xpyp"
            return
          end if
      end if

      !----------------------------------------------------------------
      ! Covariance clipping for wprtp, wpthlp, wpsclrp, upwp, and vpwp
      ! after subroutine advance_xp2_xpyp updated xp2.
      !----------------------------------------------------------------

      wprtp_cl_num   = 2 ! Second instance of w'r_t' clipping.
      wpthlp_cl_num  = 2 ! Second instance of w'th_l' clipping.
      wpsclrp_cl_num = 2 ! Second instance of w'sclr' clipping.
      if ( l_predict_upwp_vpwp ) then
         upwp_cl_num = 2 ! Second instance of u'w' clipping.
         vpwp_cl_num = 2 ! Second instance of v'w' clipping.
      else
         upwp_cl_num = 1 ! First instance of u'w' clipping.
         vpwp_cl_num = 1 ! First instance of v'w' clipping.
      endif ! l_predict_upwp_vpwp

      call clip_covars_denom( dt, rtp2, thlp2, up2, vp2, wp2,           & ! intent(in)
                              sclrp2, wprtp_cl_num, wpthlp_cl_num,      & ! intent(in)
                              wpsclrp_cl_num, upwp_cl_num, vpwp_cl_num, & ! intent(in)
                              wprtp, wpthlp, upwp, vpwp, wpsclrp,       & ! intent(inout)
                              upwp_pert, vpwp_pert )                      ! intent(inout)


      !----------------------------------------------------------------
      ! Advance the 2nd- and 3rd-order moments
      !   of vertical velocity (wp2, wp3) by one timestep.
      !----------------------------------------------------------------

      ! advance_wp2_wp3_bad_wp2 ! Test error comment, DO NOT modify or move
      call advance_wp2_wp3 &
           ( dt, sfc_elevation, sigma_sqd_w, wm_zm,              & ! intent(in)
             wm_zt, a3_coef, a3_coef_zt, wp3_on_wp2, wp4,        & ! intent(in)
             wpthvp, wp2thvp, um, vm, upwp, vpwp,                & ! intent(in)
             up2, vp2, Kh_zm, Kh_zt, tau_zm, tau_zt,             & ! intent(in)
             tau_C1_zm, Skw_zm, Skw_zt, rho_ds_zm,               & ! intent(in)
             rho_ds_zt, invrs_rho_ds_zm,                         & ! intent(in)
             invrs_rho_ds_zt, radf, thv_ds_zm,                   & ! intent(in)
             thv_ds_zt, pdf_params%mixt_frac, Cx_fnc_Richardson, & ! intent(in)
             wp2_splat, wp3_splat,                               & ! intent(in)
             pdf_implicit_coefs_terms,                           & ! intent(in)
             wprtp, wpthlp, rtp2, thlp2,                         & ! intent(in)
             wp2, wp3, wp3_zm, wp2_zt )                            ! intent(inout)

      if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "Error calling advance_wp2_wp3"
            return
          end if
      end if

      !----------------------------------------------------------------
      ! Covariance clipping for wprtp, wpthlp, wpsclrp, upwp, and vpwp
      ! after subroutine advance_wp2_wp3 updated wp2.
      !----------------------------------------------------------------

      wprtp_cl_num   = 3 ! Third instance of w'r_t' clipping.
      wpthlp_cl_num  = 3 ! Third instance of w'th_l' clipping.
      wpsclrp_cl_num = 3 ! Third instance of w'sclr' clipping.
      if ( l_predict_upwp_vpwp ) then
         upwp_cl_num = 3 ! Third instance of u'w' clipping.
         vpwp_cl_num = 3 ! Third instance of v'w' clipping.
      else
         upwp_cl_num = 2 ! Second instance of u'w' clipping.
         vpwp_cl_num = 2 ! Second instance of v'w' clipping.
      endif ! l_predict_upwp_vpwp

      call clip_covars_denom( dt, rtp2, thlp2, up2, vp2, wp2,           & ! intent(in)
                              sclrp2, wprtp_cl_num, wpthlp_cl_num,      & ! intent(in)
                              wpsclrp_cl_num, upwp_cl_num, vpwp_cl_num, & ! intent(in)
                              wprtp, wpthlp, upwp, vpwp, wpsclrp,       & ! intent(inout)
                              upwp_pert, vpwp_pert )                      ! intent(inout)

      !----------------------------------------------------------------
      ! Advance or otherwise calculate <thl'^3>, <rt'^3>, and
      ! <sclr'^3>.
      !----------------------------------------------------------------
      if ( l_advance_xp3 ) then

         ! Advance <rt'^3>, <thl'^3>, and <sclr'^3> one model timestep using a
         ! simplified form of the <x'^3> predictive equation.  The simplified
         ! <x'^3> equation can either be advanced from its previous value or
         ! calculated using a steady-state approximation.
         call advance_xp3( dt, rtm, thlm, rtp2, thlp2, wprtp,  & ! Intent(in)
                           wpthlp, wprtp2, wpthlp2, rho_ds_zm, & ! Intent(in)
                           invrs_rho_ds_zt, tau_zt,            & ! Intent(in)
                           sclrm, sclrp2, wpsclrp, wpsclrp2,   & ! Intent(in)
                           rtp3, thlp3, sclrp3                 ) ! Intent(inout)

      else

         ! Use the Larson and Golaz (2005) ansatz for the ADG1 PDF to
         ! calculate <rt'^3>, <thl'^3>, and <sclr'^3>.
         Skw_zt(1:gr%nz) = Skx_func( wp2_zt(1:gr%nz), wp3(1:gr%nz), w_tol )

         wpthlp_zt = zm2zt( wpthlp )
         wprtp_zt  = zm2zt( wprtp )
         thlp2_zt  = max( zm2zt( thlp2 ), thl_tol**2 ) ! Positive def. quantity
         rtp2_zt   = max( zm2zt( rtp2 ), rt_tol**2 )   ! Positive def. quantity

         sigma_sqd_w_zt = max( zm2zt( sigma_sqd_w ), zero_threshold )

         thlp3 = xp3_LG_2005_ansatz( Skw_zt, wpthlp_zt, wp2_zt, &
                                     thlp2_zt, sigma_sqd_w_zt, thl_tol )

         rtp3 = xp3_LG_2005_ansatz( Skw_zt, wprtp_zt, wp2_zt, &
                                    rtp2_zt, sigma_sqd_w_zt, rt_tol )

         do i = 1, sclr_dim, 1

            wpsclrp_zt = zm2zt( wpsclrp(:,i) )
            sclrp2_zt  = max( zm2zt( sclrp2(:,i) ), sclr_tol(i)**2 )

            sclrp3(:,i) = xp3_LG_2005_ansatz( Skw_zt, wpsclrp_zt, wp2_zt, &
                                              sclrp2_zt, sigma_sqd_w_zt, &
                                              sclr_tol(i) )

         enddo ! i = 1, sclr_dim

      endif ! l_advance_xp3

      !----------------------------------------------------------------
      ! Advance the horizontal mean winds (um, vm),
      !   the mean of the eddy-diffusivity scalars (i.e. edsclrm),
      !   and their fluxes (upwp, vpwp, wpedsclrp) by one time step.
      !----------------------------------------------------------------

      if ( l_use_buoy_mod_Km_zm ) then

         tau_factor = ( ( one - C5 ) / C4 ) * tau_zm
         Km_zm_denom_term = tau_factor * ( grav / T0 ) * &
                              wpthvp / max( 10._core_rknd*w_tol_sqd, wp2 )
         Km_zm_numerator_term = 0.02_core_rknd * 0.5_core_rknd * ( grav / T0 ) &
                                * tau_zm**2 * ddzt( thlm )
         Km_zm = c_K10 * tau_factor * wp2 * &
                           ( one - min( 0.9_core_rknd, Km_zm_numerator_term ) ) / &
                           ( one - min( 0.9_core_rknd, Km_zm_denom_term ) )
         ! Old method to account for upgradient contribution due to cumuli
         !Km_Skw_factor = exp( - (Skw_zm - Km_Skw_thresh) / Km_Skw_factor_efold )
         !Km_Skw_factor = max( Km_Skw_factor_min, Km_Skw_factor )
         !Km_Skw_factor = min( one, Km_Skw_factor )
         !Km_zm = Km_zm * Km_Skw_factor
      else

        Km_zm = Kh_zm * c_K10   ! Coefficient for momentum

      end if

      Kmh_zm = Kh_zm * c_K10h ! Coefficient for thermo

      if ( l_do_expldiff_rtm_thlm ) then
        edsclrm(:,edsclr_dim-1)=thlm(:)
        edsclrm(:,edsclr_dim)=rtm(:)
      endif


      call advance_windm_edsclrm( dt, wm_zt, Km_zm, Kmh_zm, ug, vg, um_ref, vm_ref, & ! intent(in)
                                  wp2, up2, vp2, um_forcing, vm_forcing,        & ! intent(in)
                                  edsclrm_forcing,                              & ! intent(in)
                                  rho_ds_zm, invrs_rho_ds_zt,                   & ! intent(in)
                                  fcor, l_implemented,                          & ! intent(in)
                                  um, vm, edsclrm,                              & ! intent(inout)
                                  upwp, vpwp, wpedsclrp, &                        ! intent(inout)
                                  um_pert, vm_pert, upwp_pert, vpwp_pert)         ! intent(inout)

      if ( l_do_expldiff_rtm_thlm ) then
        call pvertinterp(gr%nz, p_in_Pa, 70000.0_core_rknd, thlm, thlm700)
        call pvertinterp(gr%nz, p_in_Pa, 100000.0_core_rknd, thlm, thlm1000)
        if ( thlm700 - thlm1000 < 20.0_core_rknd ) then
          thlm(:) = edsclrm(:,edsclr_dim-1)
          rtm(:) = edsclrm(:,edsclr_dim)
        end if
      end if

      ! Eric Raut: this seems dangerous to call without any attached flag.
      ! Hence the preprocessor.
#ifdef CLUBB_CAM
      do ixind=1,edsclr_dim
        call fill_holes_vertical(2,0.0_core_rknd,"zt",rho_ds_zt,rho_ds_zm,edsclrm(:,ixind))
      enddo
#endif

    if ( ipdf_call_placement == ipdf_post_advance_fields &
         .or. ipdf_call_placement == ipdf_pre_post_advance_fields ) then

       ! Sample stats in this call to subroutine pdf_closure_driver for
       ! ipdf_post_advance_fields, but not for ipdf_pre_post_advance_fields
       ! because stats were sampled during the first call to subroutine
       ! pdf_closure_driver.
       if ( ipdf_call_placement == ipdf_post_advance_fields ) then
          l_samp_stats_in_pdf_call = .true.
       elseif ( ipdf_call_placement == ipdf_pre_post_advance_fields ) then
          l_samp_stats_in_pdf_call = .false.
       endif

       !########################################################################
       !#######                     CALL CLUBB's PDF                     #######
       !#######   AND OUTPUT PDF PARAMETERS AND INTEGRATED QUANTITITES   #######
       !########################################################################
       ! Given CLUBB's prognosed moments, diagnose CLUBB's PDF parameters
       !   and quantities integrated over that PDF, including
       !   quantities related to clouds, buoyancy, and turbulent advection.
       call pdf_closure_driver( dt, hydromet_dim, wprtp,       & ! Intent(in)
                                thlm, wpthlp, rtp2, rtp3,      & ! Intent(in)
                                thlp2, thlp3, rtpthlp, wp2,    & ! Intent(in)
                                wp3, wm_zm, wm_zt,             & ! Intent(in)
                                um, up2, upwp,                 & ! Intent(in)
                                vm, vp2, vpwp,                 & ! Intent(in)
                                p_in_Pa, exner,                & ! Intent(in)
                                thv_ds_zm, thv_ds_zt,          & ! Intent(in)
                                rfrzm, hydromet, wphydrometp,  & ! Intent(in)
                                wp2hmp, rtphmp_zt, thlphmp_zt, & ! Intent(in)
                                sclrm, wpsclrp, sclrp2,        & ! Intent(in)
                                sclrprtp, sclrpthlp,           & ! Intent(in)
                                l_samp_stats_in_pdf_call,      & ! Intent(in)
                                rtm,                           & ! Intent(i/o)
#ifdef GFDL
                                RH_crit(k, : , :),             & ! Intent(i/o)
                                do_liquid_only_in_clubb,       & ! Intent(in)
#endif
                                rcm, cloud_frac,               & ! Intent(out)
                                ice_supersat_frac, wprcp,      & ! Intent(out)
                                sigma_sqd_w, wpthvp, wp2thvp,  & ! Intent(out)
                                rtpthvp, thlpthvp, rc_coef,    & ! Intent(out)
                                rcm_in_layer, cloud_cover,     & ! Intent(out)
                                rcp2_zt, thlprcp, rc_coef_zm,  & ! Intent(out)
                                rtm_frz, thlm_frz, sclrpthvp,  & ! Intent(out)
                                wp4, wp2rtp, wprtp2, wp2thlp,  & ! Intent(out)
                                wpthlp2, wprtpthlp, wp2rcp,    & ! Intent(out)
                                rtprcp, rcp2,                  & ! Intent(out)
                                uprcp, vprcp,                  & ! Intent(out)
                                Skw_velocity,                  & ! Intent(out)
                                cloud_frac_zm,                 & ! Intent(out)
                                ice_supersat_frac_zm,          & ! Intent(out)
                                rtm_zm, thlm_zm, rcm_zm,       & ! Intent(out)
                                rcm_supersat_adj,              & ! Intent(out)
                                wp2sclrp, wpsclrp2, sclrprcp,  & ! Intent(out)
                                wpsclrprtp, wpsclrpthlp,       & ! Intent(out)
                                pdf_params, pdf_params_frz,    & ! Intent(out)
                                pdf_params_zm,                 & ! Intent(out)
                                pdf_implicit_coefs_terms )       ! Intent(out)

    endif ! ipdf_call_placement == ipdf_post_advance_fields
          ! or ipdf_call_placement == ipdf_pre_post_advance_fields

#ifdef CLUBB_CAM
      qclvar(:) = rcp2_zt(:)
      thlprcp_out(:) = thlprcp(:)
#endif


      !#######################################################################
      !#############            ACCUMULATE STATISTICS            #############
      !#######################################################################

      if ( l_stats_samp ) then

         call stat_end_update( iwp2_bt, wp2 / dt, & ! intent(in)
                               stats_zm )           ! intent(inout)
         call stat_end_update( ivp2_bt, vp2 / dt, & ! intent(in)
                               stats_zm )           ! intent(inout)
         call stat_end_update( iup2_bt, up2 / dt, & ! intent(in)
                               stats_zm )           ! intent(inout)
         call stat_end_update( iwprtp_bt, wprtp / dt, & ! intent(in)
                               stats_zm )               ! intent(inout)
         call stat_end_update( iwpthlp_bt, wpthlp / dt, & ! intent(in)
                               stats_zm )                 ! intent(inout)
         if ( l_predict_upwp_vpwp ) then
            call stat_end_update( iupwp_bt, upwp / dt, & ! intent(in)
                                  stats_zm )             ! intent(inout)
            call stat_end_update( ivpwp_bt, vpwp / dt, & ! intent(in)
                                  stats_zm )             ! intent(inout)
         endif ! l_predict_upwp_vpwp
         call stat_end_update( irtp2_bt, rtp2 / dt, & ! intent(in)
                               stats_zm )             ! intent(inout)
         call stat_end_update( ithlp2_bt, thlp2 / dt, & ! intent(in)
                               stats_zm )               ! intent(inout)
         call stat_end_update( irtpthlp_bt, rtpthlp / dt, & ! intent(in)
                               stats_zm )                   ! intent(inout)

         call stat_end_update( irtm_bt, rtm / dt, & ! intent(in)
                               stats_zt )           ! intent(inout)
         call stat_end_update( ithlm_bt, thlm / dt, & ! intent(in)
                               stats_zt )             ! intent(inout)
         call stat_end_update( ium_bt, um / dt, & ! intent(in)
                               stats_zt )         ! intent(inout)
         call stat_end_update( ivm_bt, vm / dt, & ! intent(in)
                               stats_zt )         ! intent(inout)
         call stat_end_update( iwp3_bt, wp3 / dt, & ! intent(in)
                               stats_zt )           ! intent(inout)

      endif ! l_stats_samp


      if ( iwpthlp_zt > 0 ) then
        wpthlp_zt  = zm2zt( wpthlp )
      end if

      if ( iwprtp_zt > 0 ) then
        wprtp_zt   = zm2zt( wprtp )
      end if

      if ( iup2_zt > 0 ) then
        up2_zt = max( zm2zt( up2 ), w_tol_sqd )
      end if

      if (ivp2_zt > 0 ) then
        vp2_zt = max( zm2zt( vp2 ), w_tol_sqd )
      end if

      if ( iupwp_zt > 0 ) then
        upwp_zt = zm2zt( upwp )
      end if

      if ( ivpwp_zt > 0 ) then
        vpwp_zt = zm2zt( vpwp )
      end if

      call stats_accumulate &
           ( um, vm, upwp, vpwp, up2, vp2,                          & ! intent(in)
             thlm, rtm, wprtp, wpthlp,                              & ! intent(in)
             wp2, wp3, rtp2, rtp3, thlp2, thlp3, rtpthlp,           & ! intent(in)
             wpthvp, wp2thvp, rtpthvp, thlpthvp,                    & ! intent(in)
             p_in_Pa, exner, rho, rho_zm,                           & ! intent(in)
             rho_ds_zm, rho_ds_zt, thv_ds_zm, thv_ds_zt,            & ! intent(in)
             wm_zt, wm_zm, rcm, wprcp, rc_coef, rc_coef_zm,         & ! intent(in)
             rcm_zm, rtm_zm, thlm_zm, cloud_frac, ice_supersat_frac,& ! intent(in)
             cloud_frac_zm, ice_supersat_frac_zm, rcm_in_layer,     & ! intent(in)
             cloud_cover, rcm_supersat_adj, sigma_sqd_w,            & ! intent(in)
             pdf_params, pdf_params_zm, sclrm, sclrp2,              & ! intent(in)
             sclrprtp, sclrpthlp, sclrm_forcing, sclrpthvp,         & ! intent(in)
             wpsclrp, edsclrm, edsclrm_forcing                      ) ! intent(in)


      if ( clubb_at_least_debug_level( 2 ) ) then
        call parameterization_check &
           ( thlm_forcing, rtm_forcing, um_forcing,                             & ! intent(in)
             vm_forcing, wm_zm, wm_zt, p_in_Pa,                                 & ! intent(in)
             rho_zm, rho, exner, rho_ds_zm,                                     & ! intent(in)
             rho_ds_zt, invrs_rho_ds_zm, invrs_rho_ds_zt,                       & ! intent(in)
             thv_ds_zm, thv_ds_zt, wpthlp_sfc, wprtp_sfc, upwp_sfc,             & ! intent(in)
             vpwp_sfc, um, upwp, vm, vpwp, up2, vp2,                            & ! intent(in)
             rtm, wprtp, thlm, wpthlp, wp2, wp3,                                & ! intent(in)
             rtp2, thlp2, rtpthlp, rcm,                                         & ! intent(in)
            "end of ",                                                          & ! intent(in)
             wpsclrp_sfc, wpedsclrp_sfc, sclrm, wpsclrp, sclrp2,                & ! intent(in)
             sclrprtp, sclrpthlp, sclrm_forcing, edsclrm, edsclrm_forcing       ) ! intent(in)

        if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) return
        end if

      end if

      if ( l_stats .and. l_stats_samp ) then
        ! Spurious source will only be calculated if rtm_ma and thlm_ma are zero.
        ! Therefore, wm must be zero or l_implemented must be true.
        if ( l_implemented .or. &
             (all( abs(wm_zt) < eps ) .and. all( abs(wm_zm) < eps ))) then
          ! Calculate the spurious source for rtm
          rtm_flux_top = rho_ds_zm(gr%nz) * wprtp(gr%nz)

          if ( .not. l_host_applies_sfc_fluxes ) then
            rtm_flux_sfc = rho_ds_zm(1) * wprtp_sfc
          else
            rtm_flux_sfc = 0.0_core_rknd
          end if

          rtm_integral_after  &
          = vertical_integral( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                               rtm(2:gr%nz), gr%dzt(2:gr%nz) )

          rtm_integral_forcing  &
          = vertical_integral( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                               rtm_forcing(2:gr%nz), gr%dzt(2:gr%nz) )

          rtm_spur_src  &
          = calculate_spurious_source( rtm_integral_after, &
                                       rtm_integral_before, &
                                       rtm_flux_top, rtm_flux_sfc, &
                                       rtm_integral_forcing, &
                                       dt )

          ! Calculate the spurious source for thlm
          thlm_flux_top = rho_ds_zm(gr%nz) * wpthlp(gr%nz)

          if ( .not. l_host_applies_sfc_fluxes ) then
            thlm_flux_sfc = rho_ds_zm(1) * wpthlp_sfc
          else
            thlm_flux_sfc = 0.0_core_rknd
          end if

          thlm_integral_after  &
          = vertical_integral( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                               thlm(2:gr%nz), gr%dzt(2:gr%nz) )

          thlm_integral_forcing  &
          = vertical_integral( (gr%nz - 2 + 1), rho_ds_zt(2:gr%nz), &
                               thlm_forcing(2:gr%nz), gr%dzt(2:gr%nz) )

          thlm_spur_src  &
          = calculate_spurious_source( thlm_integral_after, &
                                       thlm_integral_before, &
                                       thlm_flux_top, thlm_flux_sfc, &
                                       thlm_integral_forcing, &
                                       dt )
        else ! If l_implemented is false, we don't want spurious source output
          rtm_spur_src = -9999.0_core_rknd
          thlm_spur_src = -9999.0_core_rknd
        end if

        ! Write the var to stats
        call stat_update_var_pt( irtm_spur_src, 1, rtm_spur_src,   & ! intent(in)
                                 stats_sfc )                               ! intent(inout)
        call stat_update_var_pt( ithlm_spur_src, 1, thlm_spur_src, & ! intent(in)
                                 stats_sfc )                               ! intent(inout)
      end if

      return
    end subroutine advance_clubb_core

  !=============================================================================
  subroutine pdf_closure_driver( dt, hydromet_dim, wprtp,       & ! Intent(in)
                                 thlm, wpthlp, rtp2, rtp3,      & ! Intent(in)
                                 thlp2, thlp3, rtpthlp, wp2,    & ! Intent(in)
                                 wp3, wm_zm, wm_zt,             & ! Intent(in)
                                 um, up2, upwp,                 & ! Intent(in)
                                 vm, vp2, vpwp,                 & ! Intent(in)
                                 p_in_Pa, exner,                & ! Intent(in)
                                 thv_ds_zm, thv_ds_zt,          & ! Intent(in)
                                 rfrzm, hydromet, wphydrometp,  & ! Intent(in)
                                 wp2hmp, rtphmp_zt, thlphmp_zt, & ! Intent(in)
                                 sclrm, wpsclrp, sclrp2,        & ! Intent(in)
                                 sclrprtp, sclrpthlp,           & ! Intent(in)
                                 l_samp_stats_in_pdf_call,      & ! Intent(in)
                                 rtm,                           & ! Intent(i/o)
#ifdef GFDL
                                 RH_crit(k, : , :),             & ! Intent(i/o)
                                 do_liquid_only_in_clubb,       & ! Intent(in)
#endif
                                 rcm, cloud_frac,               & ! Intent(out)
                                 ice_supersat_frac, wprcp,      & ! Intent(out)
                                 sigma_sqd_w, wpthvp, wp2thvp,  & ! Intent(out)
                                 rtpthvp, thlpthvp, rc_coef,    & ! Intent(out)
                                 rcm_in_layer, cloud_cover,     & ! Intent(out)
                                 rcp2_zt, thlprcp, rc_coef_zm,  & ! Intent(out)
                                 rtm_frz, thlm_frz, sclrpthvp,  & ! Intent(out)
                                 wp4, wp2rtp, wprtp2, wp2thlp,  & ! Intent(out)
                                 wpthlp2, wprtpthlp, wp2rcp,    & ! Intent(out)
                                 rtprcp, rcp2,                  & ! Intent(out)
                                 uprcp, vprcp,                  & ! Intent(out)
                                 Skw_velocity,                  & ! Intent(out)
                                 cloud_frac_zm,                 & ! Intent(out)
                                 ice_supersat_frac_zm,          & ! Intent(out)
                                 rtm_zm, thlm_zm, rcm_zm,       & ! Intent(out)
                                 rcm_supersat_adj,              & ! Intent(out)
                                 wp2sclrp, wpsclrp2, sclrprcp,  & ! Intent(out)
                                 wpsclrprtp, wpsclrpthlp,       & ! Intent(out)
                                 pdf_params, pdf_params_frz,    & ! Intent(out)
                                 pdf_params_zm,                 & ! Intent(out)
                                 pdf_implicit_coefs_terms )       ! Intent(out)

    use grid_class, only: &
        gr,    & ! Variable(s)
        zt2zm, & ! Procedure(s)
        zm2zt

    use constants_clubb, only: &
        w_tol,          & ! Variable(s)
        w_tol_sqd,      &
        rt_tol,         &
        thl_tol,        &
        Cp,             &
        Lv,             &
        Ls,             &
        p0,             &
        kappa,          &
        fstderr,        &
        zero,           &
        zero_threshold, &
        eps

    use pdf_parameter_module, only: &
        pdf_parameter,        & ! Variable Type
        implicit_coefs_terms, & ! Variable Type
        init_pdf_params         ! Procedure(s)

    use parameters_model, only: &
        sclr_dim,               & ! Variable(s)
        sclr_tol,               &
        ts_nudge,               &
        rtm_min,                &
        rtm_nudge_max_altitude

    use parameters_tunable, only: &
        gamma_coef,  & ! Variable(s)
        gamma_coefb, &
        gamma_coefc

    use pdf_closure_module, only: &
        pdf_closure,                & ! Procedure(s)
        calc_vert_avg_cf_component

    use Skx_module, only: &
        Skx_func    ! Procedure(s)

    use sigma_sqd_w_module, only: &
        compute_sigma_sqd_w    ! Procedure(s)

    use pdf_utilities, only: &
        compute_mean_binormal    ! Procedure(s)

    use T_in_K_module, only: &
        thlm2T_in_K    ! Procedure(s)

    use saturation, only:  &
        sat_mixrat_liq    ! Procedure(s)

    use array_index, only: &
        iirr    ! Variable(s)

    use model_flags, only: &
        l_gamma_Skw,                  & ! Variable(s)
        l_call_pdf_closure_twice,     &
        l_trapezoidal_rule_zm,        &
        l_trapezoidal_rule_zt,        &
        l_use_cloud_cover,            &
        l_use_ice_latent,             &
        l_rcm_supersat_adj,           &
        l_rtm_nudge,                  &
        l_explicit_turbulent_adv_wp3

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use variables_diagnostic_module, only: &
        sigma_sqd_w_zt, & ! Variable(s)
        rtm_ref

    use stats_type_utilities, only: &
        stat_update_var,    & ! Procedure(s)
        stat_update_var_pt

    use stats_variables, only: &
        l_stats_samp,        & ! Variable(s)
        stats_zt,            &
        stats_zm,            &
        iSkw_zm,             &
        iSkw_zt,             &
        iSkrt_zm,            &
        iSkrt_zt,            &
        iSkthl_zm,           &
        iSkthl_zt,           &
        iSkw_velocity,       &
        igamma_Skw_fnc,      &
        iF_w,                &
        iF_rt,               &
        iF_thl,              &
        imin_F_w,            &
        imax_F_w,            &
        imin_F_rt,           &
        imax_F_rt,           &
        imin_F_thl,          &
        imax_F_thl,          &
        ircp2,               &
        iwp4,                &
        ircm_refined,        &
        icloud_frac_refined

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    !!! External
    intrinsic :: sqrt, min, max, exp, mod, real

    logical, parameter :: &
      l_refine_grid_in_cloud = .false., & ! Compute cloud_frac and rcm on a refined grid

      l_interactive_refined  = .false.    ! Should the refined grid code feed into the model?
                                          ! Only has meaning if l_refined_grid_in_cloud is .true.

    real( kind = core_rknd ), parameter :: &
      chi_at_liq_sat = 0._core_rknd  ! Value of chi(s) at saturation with
                                     ! respect to ice (zero for liquid)

    !!! Input Variables
    real( kind = core_rknd ), intent(in) ::  &
      dt  ! Current timestep duration    [s]

    integer, intent(in) :: &
      hydromet_dim      ! Total number of hydrometeors          [#]

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
!      rtm,       & ! total water mixing ratio, r_t (thermo. levels) [kg/kg]
      wprtp,     & ! w' r_t' (momentum levels)                      [(kg/kg)m/s]
      thlm,      & ! liq. water pot. temp., th_l (thermo. levels)   [K]
      wpthlp,    & ! w' th_l' (momentum levels)                     [(m/s) K]
      rtp2,      & ! r_t'^2 (momentum levels)                       [(kg/kg)^2]
      rtp3,      & ! r_t'^3 (thermodynamic levels)                  [(kg/kg)^3]
      thlp2,     & ! th_l'^2 (momentum levels)                      [K^2]
      thlp3,     & ! th_l'^3 (thermodynamic levels)                 [K^3]
      rtpthlp,   & ! r_t' th_l' (momentum levels)                   [(kg/kg) K]
      wp2,       & ! w'^2 (momentum levels)                         [m^2/s^2]
      wp3,       & ! w'^3 (thermodynamic levels)                    [m^3/s^3]
      wm_zm,     & ! w mean wind component on momentum levels       [m/s]
      wm_zt,     & ! w mean wind component on thermo. levels        [m/s]
      p_in_Pa,   & ! Air pressure (thermodynamic levels)            [Pa]
      exner,     & ! Exner function (thermodynamic levels)          [-]
      thv_ds_zm, & ! Dry, base-state theta_v on momentum levs.      [K]
      thv_ds_zt, & ! Dry, base-state theta_v on thermo. levs.       [K]
      rfrzm        ! Total ice-phase water mixing ratio             [kg/kg]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      um,          & ! Grid-mean eastward wind     [m/s]
      up2,         & ! u'^2                        [(m/s)^2]
      upwp,        & ! u'w'                        [(m/s)^2]
      vm,          & ! Grid-mean northward wind    [m/s]
      vp2,         & ! v'^2                        [(m/s)^2]
      vpwp           ! v'w'                        [(m/s)^2]

    ! Hydrometeor variables
    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim), intent(in) :: &
      hydromet       ! Mean of hydrometeor fields               [units vary]

    real( kind = core_rknd ), dimension(gr%nz, hydromet_dim), intent(in) :: &
      wphydrometp, & ! Covariance of w and a hydrometeor      [(m/s) <hm units>]
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' >    [(m/s)^2 <hm units>]
      rtphmp_zt,   & ! Covariance of rt and hm (on t-levs.) [(kg/kg) <hm units>]
      thlphmp_zt     ! Covariance of thl and hm (on t-levs.)      [K <hm units>]

    ! Passive scalar variables
    real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(in) :: &
      sclrm,     & ! Passive scalar mean (thermo. levels) [units vary]
      wpsclrp,   & ! w'sclr' (momentum levels)            [{units vary} m/s]
      sclrp2,    & ! sclr'^2 (momentum levels)            [{units vary}^2]
      sclrprtp,  & ! sclr'rt' (momentum levels)           [{units vary} (kg/kg)]
      sclrpthlp    ! sclr'thl' (momentum levels)          [{units vary} K]

    logical, intent(in) :: &
      l_samp_stats_in_pdf_call    ! Sample stats in this call to this subroutine

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) ::  &
      rtm    ! total water mixing ratio, r_t (thermo. levels) [kg/kg]

#ifdef GFDL
    ! hlg, 2010-06-16
    real( kind = core_rknd ), dimension(gr%nz, min(1,sclr_dim) , 2), intent(inout) :: &
      RH_crit  ! critical relative humidity for droplet and ice nucleation
! ---> h1g, 2012-06-14
    logical, intent(in)                 ::  do_liquid_only_in_clubb
! <--- h1g, 2012-06-14
#endif

    !!! Output Variables
    ! Variables being passed back to and out of advance_clubb_core.
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      rcm,               & ! mean r_c (thermodynamic levels)        [kg/kg]
      cloud_frac,        & ! cloud fraction (thermodynamic levels)  [-]
      ice_supersat_frac, & ! ice supersat. frac. (thermo. levels)   [-]
      wprcp,             & ! < w'r_c' > (momentum levels)           [m/s kg/kg]
      sigma_sqd_w,       & ! PDF width parameter (momentum levels)  [-]
      wpthvp,            & ! < w' th_v' > (momentum levels)         [kg/kg K]
      wp2thvp,           & ! < w'^2 th_v' > (thermodynamic levels)  [m^2/s^2 K]
      rtpthvp,           & ! < r_t' th_v' > (momentum levels)       [kg/kg K]
      thlpthvp,          & ! < th_l' th_v' > (momentum levels)      [K^2]
      rc_coef,           & ! Coefficient of X'r_c' (thermo. levs.)  [K/(kg/kg)]
      rcm_in_layer,      & ! rcm in cloud layer                     [kg/kg]
      cloud_cover,       & ! cloud cover                            [-]
      rcp2_zt,           & ! r_c'^2 (on thermo. grid)               [kg^2/kg^2]
      thlprcp,           & ! < th_l' r_c' > (momentum levels)       [K kg/kg]
      rc_coef_zm,        & ! Coefficient of X'r_c' on m-levs.       [K/(kg/kg)]
      rtm_frz,           & ! rtm adjusted to include hydrometeors   [kg/kg]
      thlm_frz             ! thlm adjusted to include hydrometeors  [K]

    ! Variable being passed back to and out of advance_clubb_core.
    real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(out) :: &
      sclrpthvp    ! < sclr' th_v' > (momentum levels)   [units vary]

    ! Variables being passed back to only advance_clubb_core (for statistics).
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      wp4,       & ! < w'^4 > (momentum levels)               [m^4/s^4]
      wp2rtp,    & ! < w'^2 r_t' > (thermodynamic levels)     [m^2/s^2 kg/kg]
      wprtp2,    & ! < w' r_t'^2 > (thermodynamic levels)     [m/s kg^2/kg^2]
      wp2thlp,   & ! < w'^2 th_l' > (thermodynamic levels)    [m^2/s^2 K]
      wpthlp2,   & ! < w' th_l'^2 > (thermodynamic levels)    [m/s K^2]
      wprtpthlp, & ! < w' r_t' th_l' > (thermodynamic levels) [m/s kg/kg K]
      wp2rcp,    & ! < w'^2 r_c' > (thermodynamic levels)     [m^2/s^2 kg/kg]
      rtprcp,    & ! < r_t' r_c' > (momentum levels)          [kg^2/kg^2]
      rcp2         ! Variance of r_c (momentum levels)        [kg^2/kg^2]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      uprcp,              & ! < u' r_c' >              [(m kg)/(s kg)]
      vprcp                 ! < v' r_c' >              [(m kg)/(s kg)]

    ! Variables being passed back to only advance_clubb_core (for statistics).
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      Skw_velocity,         & ! Skewness velocity                        [m/s]
      cloud_frac_zm,        & ! Cloud Fraction on momentum levels        [-]
      ice_supersat_frac_zm, & ! Ice supersat. frac. on momentum levels   [-]
      rtm_zm,               & ! Total water mixing ratio at mom. levs.   [kg/kg]
      thlm_zm,              & ! Liquid water pot. temp. at mom. levs.    [K]
      rcm_zm,               & ! rcm at momentum levels                   [kg/kg]
      rcm_supersat_adj        ! Adjust. to rcm due to spurious supersat. [kg/kg]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(out) :: &
      wp2sclrp,    & ! < w'^2 sclr' > (thermodynamic levels)      [units vary]
      wpsclrp2,    & ! < w' sclr'^2 > (thermodynamic levels)      [units vary]
      sclrprcp,    & ! < sclr' r_c' > (momentum levels)           [units vary]
      wpsclrprtp,  & ! < w' sclr' r_t' > (thermodynamic levels)   [units vary]
      wpsclrpthlp    ! < w' sclr' th_l' > (thermodynamic levels)  [units vary]

    ! Variable being passed back to and out of advance_clubb_core.
    type(pdf_parameter), intent(inout) :: &
      pdf_params,     & ! PDF parameters                          [units vary]
      pdf_params_frz    ! Output for use in pert. Lscale calc.

    ! Variable being passed back to only advance_clubb_core.
    type(pdf_parameter), intent(inout) :: &
      pdf_params_zm    ! PDF parameters   [units vary]

    type(implicit_coefs_terms), dimension(gr%nz), intent(out) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    !!! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      wp2_zt,        & ! wp2 interpolated to thermodynamic levels    [m^2/s^2]
      wp3_zm,        & ! wp3 interpolated to momentum levels         [m^3/s^3]
      rtp2_zt,       & ! rtp2 interpolated to thermodynamic levels   [kg^2/kg^2]
      rtp3_zm,       & ! rtp3 interpolated to momentum levels        [kg^3/kg^3]
      thlp2_zt,      & ! thlp2 interpolated to thermodynamic levels  [K^2]
      thlp3_zm,      & ! thlp3 interpolated to momentum levels       [K^3]
      rtpthlp_zt,    & ! rtpthlp interp. to thermodynamic levels     [kg/kg K]
      gamma_Skw_fnc, & ! Gamma as a function of skewness             [-]
      Skw_zt,        & ! Skewness of w on thermodynamic levels       [-]
      Skw_zm,        & ! Skewness of w on momentum levels            [-]
      Skrt_zt,       & ! Skewness of rt on thermodynamic levels      [-]
      Skrt_zm,       & ! Skewness of rt on momentum levels           [-]
      Skthl_zt,      & ! Skewness of thl on thermodynamic levels     [-]
      Skthl_zm         ! Skewness of thl on momentum levels          [-]

    ! Interpolated values for optional second call to PDF closure.
    real( kind = core_rknd ), dimension(gr%nz) :: &
      p_in_Pa_zm, & ! Pressure interpolated to momentum levels  [Pa]
      exner_zm      ! Exner interpolated to momentum levels     [-]

    real( kind = core_rknd ), dimension(gr%nz,hydromet_dim) :: &
      wphydrometp_zt, & ! Covariance of w and hm (on t-levs.) [(m/s) <hm units>]
      wp2hmp_zm,      & ! Moment <w'^2 hm'> (on m-levs.)    [(m/s)^2 <hm units>]
      rtphmp,         & ! Covariance of rt and hm           [(kg/kg) <hm units>]
      thlphmp           ! Covariance of thl and hm                [K <hm units>]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      wpsclrp_zt,   & ! w' sclr' on thermo. levels
      sclrp2_zt,    & ! sclr'^2 on thermo. levels
      sclrprtp_zt,  & ! sclr' r_t' on thermo. levels
      sclrpthlp_zt    ! sclr' th_l' on thermo. levels

    ! These local variables are declared because they originally belong on the
    ! momentum grid levels, but pdf_closure outputs them on the thermodynamic
    ! grid levels.
    real( kind = core_rknd ), dimension(gr%nz) :: &
      wp4_zt,      & ! w'^4 (on thermo. grid)           [m^4/s^4]
      wpthvp_zt,   & ! Buoyancy flux (on thermo. grid)  [(K m)/s]
      rtpthvp_zt,  & ! r_t' th_v' (on thermo. grid)     [(kg K)/kg]
      thlpthvp_zt, & ! th_l' th_v' (on thermo. grid)    [K^2]
      wprcp_zt,    & ! w' r_c' (on thermo. grid)        [(m kg)/(s kg)]
      rtprcp_zt,   & ! r_t' r_c' (on thermo. grid)      [(kg^2)/(kg^2)]
      thlprcp_zt,  & ! th_l' r_c' (on thermo. grid)     [(K kg)/kg]
      uprcp_zt,    & ! u' r_c' (on thermo. grid)        [(m kg)/(s kg)]
      vprcp_zt       ! v' r_c' (on thermo. grid)        [(m kg)/(s kg)]

    real( kind = core_rknd ), dimension(gr%nz, sclr_dim) :: &
      sclrpthvp_zt, & ! sclr'th_v' (on thermo. grid)
      sclrprcp_zt     ! sclr'rc' (on thermo. grid)

    real( kind = core_rknd ), dimension(gr%nz) :: &
      wprtp2_zm,    & ! < w' r_t'^2 > on momentum levels      [m/s kg^2/kg^2]
      wp2rtp_zm,    & ! < w'^2 r_t' > on momentum levels      [m^2/s^2 kg/kg]
      wpthlp2_zm,   & ! < w' th_l'^2 > on momentum levels     [m/s K^2]
      wp2thlp_zm,   & ! < w'^2 th_l' > on momentum levels     [m^2/s^2 K]
      wprtpthlp_zm, & ! < w' r_t' th_l' > on momentum levels  [m/s kg/kg K]
      wp2thvp_zm,   & ! < w'^2 th_v' > on momentum levels     [m^2/s^2 K]
      wp2rcp_zm       ! < w'^2 r_c' > on momentum levles      [m^2/s^2 kg/kg]

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      wpsclrprtp_zm,  & ! w'sclr'rt' on momentum grid
      wpsclrp2_zm,    & ! w'sclr'^2 on momentum grid
      wpsclrpthlp_zm, & ! w'sclr'thl' on momentum grid
      wp2sclrp_zm,    & ! w'^2 sclr' on momentum grid
      sclrm_zm          ! Passive scalar mean on momentum grid

    ! Output from new PDF for recording statistics.
    real( kind = core_rknd ), dimension(gr%nz) :: &
      F_w,   & ! Parameter for the spread of the PDF component means of w    [-]
      F_rt,  & ! Parameter for the spread of the PDF component means of rt   [-]
      F_thl    ! Parameter for the spread of the PDF component means of thl  [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      min_F_w,   & ! Minimum allowable value of parameter F_w      [-]
      max_F_w,   & ! Maximum allowable value of parameter F_w      [-]
      min_F_rt,  & ! Minimum allowable value of parameter F_rt     [-]
      max_F_rt,  & ! Maximum allowable value of parameter F_rt     [-]
      min_F_thl, & ! Minimum allowable value of parameter F_thl    [-]
      max_F_thl    ! Maximum allowable value of parameter F_thl    [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      F_w_zm,       &
      F_rt_zm,      &
      F_thl_zm,     &
      min_F_w_zm,   &
      max_F_w_zm,   &
      min_F_rt_zm,  &
      max_F_rt_zm,  &
      min_F_thl_zm, &
      max_F_thl_zm

    type(implicit_coefs_terms), dimension(gr%nz) :: &
      pdf_implicit_coefs_terms_zm,     &
      pdf_implicit_coefs_terms_frz,    &
      pdf_implicit_coefs_terms_zm_frz

    ! The following variables are defined for use when l_use_ice_latent = .true.
    type(pdf_parameter ):: &
      pdf_params_zm_frz

    real( kind = core_rknd ), dimension(gr%nz)  :: &
      wp4_zt_frz, &
      wprtp2_frz, &
      wp2rtp_frz, &
      wpthlp2_frz, &
      wp2thlp_frz, &
      wprtpthlp_frz, &
      cloud_frac_frz, &
      ice_supersat_frac_frz, &
      rcm_frz, &
      wpthvp_frz, &
      wpthvp_zt_frz, &
      wp2thvp_frz, &
      wp2thvp_zm_frz, &
      rtpthvp_frz, &
      rtpthvp_zt_frz, &
      thlpthvp_frz, &
      thlpthvp_zt_frz, &
      wprcp_zt_frz, &
      wp2rcp_frz, &
      uprcp_zt_frz, &
      vprcp_zt_frz

    real( kind = core_rknd ), dimension(gr%nz)  :: &
      rtprcp_zt_frz, &
      thlprcp_zt_frz, &
      rcp2_zt_frz, &
      rc_coef_frz, &
      wp4_frz, &
      wprtp2_zm_frz, &
      wp2rtp_zm_frz, &
      wpthlp2_zm_frz, &
      wp2thlp_zm_frz, &
      wprtpthlp_zm_frz, &
      cloud_frac_zm_frz, &
      ice_supersat_frac_zm_frz, &
      rcm_zm_frz, &
      wprcp_frz, &
      wp2rcp_zm_frz, &
      rtprcp_frz, &
      thlprcp_frz, &
      rcp2_frz, &
      rtm_zm_frz, &
      thlm_zm_frz, &
      rc_coef_zm_frz, &
      uprcp_frz, &
      vprcp_frz

    real( kind = core_rknd ), dimension(gr%nz) :: &
      F_w_frz,          &
      F_rt_frz,         &
      F_thl_frz,        &
      min_F_w_frz,      &
      max_F_w_frz,      &
      min_F_rt_frz,     &
      max_F_rt_frz,     &
      min_F_thl_frz,    &
      max_F_thl_frz,    &
      F_w_zm_frz,       &
      F_rt_zm_frz,      &
      F_thl_zm_frz,     &
      min_F_w_zm_frz,   &
      max_F_w_zm_frz,   &
      min_F_rt_zm_frz,  &
      max_F_rt_zm_frz,  &
      min_F_thl_zm_frz, &
      max_F_thl_zm_frz

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) :: &
      wpsclrprtp_frz, &
      wpsclrp2_frz, &
      sclrpthvp_zt_frz, &
      wpsclrpthlp_frz, &
      sclrprcp_zt_frz, &
      wp2sclrp_frz, &
      wpsclrprtp_zm_frz, &
      wpsclrp2_zm_frz, &
      sclrpthvp_frz, &
      wpsclrpthlp_zm_frz, &
      sclrprcp_frz, &
      wp2sclrp_zm_frz

    real( kind = core_rknd ) :: &
      cloud_frac_1_refined, & ! cloud_frac_1 computed on refined grid
      cloud_frac_2_refined, & ! cloud_frac_2 computed on refined grid
      rc_1_refined, &         ! rc_1 computed on refined grid
      rc_2_refined, &         ! rc_2 computed on refined grid
      cloud_frac_refined, &   ! cloud_frac gridbox mean on refined grid
      rcm_refined             ! rcm gridbox mean on refined grid

    real( kind = core_rknd ), dimension(gr%nz) :: &
      rrm,              & ! Rain water mixing ratio
      rsat,             & ! Saturation mixing ratio from mean rt and thl.
      rel_humidity        ! Relative humidity after PDF closure [-]

    logical :: l_spur_supersat   ! Spurious supersaturation?

    integer :: i, k

    !---------------------------------------------------------------------------
    ! Interpolate wp3, rtp3, and thlp3 to momentum levels, and wp2, rtp2, and
    ! thlp2 to thermodynamic levels, and then compute Skw, Skrt, and Skthl for
    ! both the momentum and thermodynamic grid levels.
    !---------------------------------------------------------------------------

    wp2_zt   = max( zm2zt( wp2 ), w_tol_sqd ) ! Positive definite quantity
    wp3_zm   = zt2zm( wp3 )
    thlp2_zt = max( zm2zt( thlp2 ), thl_tol**2 ) ! Positive definite quantity
    thlp3_zm = zt2zm( thlp3 )
    rtp2_zt  = max( zm2zt( rtp2 ), rt_tol**2 ) ! Positive definite quantity
    rtp3_zm  = zt2zm( rtp3 )

    Skw_zt(1:gr%nz) = Skx_func( wp2_zt(1:gr%nz), wp3(1:gr%nz), w_tol )
    Skw_zm(1:gr%nz) = Skx_func( wp2(1:gr%nz), wp3_zm(1:gr%nz), w_tol )

    Skthl_zt(1:gr%nz) = Skx_func( thlp2_zt(1:gr%nz), thlp3(1:gr%nz), thl_tol )
    Skthl_zm(1:gr%nz) = Skx_func( thlp2(1:gr%nz), thlp3_zm(1:gr%nz), thl_tol )

    Skrt_zt(1:gr%nz) = Skx_func( rtp2_zt(1:gr%nz), rtp3(1:gr%nz), rt_tol )
    Skrt_zm(1:gr%nz) = Skx_func( rtp2(1:gr%nz), rtp3_zm(1:gr%nz), rt_tol )

    if ( l_stats_samp .and. l_samp_stats_in_pdf_call ) then
       call stat_update_var( iSkw_zt, Skw_zt, & ! In
                             stats_zt ) ! In/Out
       call stat_update_var( iSkw_zm, Skw_zm, &
                             stats_zm ) ! In/Out
       call stat_update_var( iSkthl_zt, Skthl_zt, &
                             stats_zt ) ! In/Out
       call stat_update_var( iSkthl_zm, Skthl_zm, &
                             stats_zm ) ! In/Out
       call stat_update_var( iSkrt_zt, Skrt_zt, &
                             stats_zt ) ! In/Out
       call stat_update_var( iSkrt_zm, Skrt_zm, &
                             stats_zm ) ! In/Out
    endif

    ! The right hand side of this conjunction is only for reducing cpu time,
    ! since the more complicated formula is mathematically equivalent
    if ( l_gamma_Skw &
         .and. abs( gamma_coef - gamma_coefb ) &
               > abs( gamma_coef + gamma_coefb ) * eps/2 ) then
      !----------------------------------------------------------------
      ! Compute gamma as a function of Skw  - 14 April 06 dschanen
      !----------------------------------------------------------------

      gamma_Skw_fnc = gamma_coefb + (gamma_coef-gamma_coefb) &
            *exp( -(1.0_core_rknd/2.0_core_rknd) * (Skw_zm/gamma_coefc)**2 )

    else

      gamma_Skw_fnc = gamma_coef

    end if

    ! Compute sigma_sqd_w (dimensionless PDF width parameter)
    sigma_sqd_w = compute_sigma_sqd_w( gamma_Skw_fnc, wp2, thlp2, rtp2, &
                                       up2, vp2, wpthlp, wprtp, upwp, vpwp )

    if ( l_stats_samp .and. l_samp_stats_in_pdf_call ) then
      call stat_update_var( igamma_Skw_fnc, gamma_Skw_fnc, & ! intent(in)
                            stats_zm )                       ! intent(inout)
    endif

    ! Smooth in the vertical using interpolation
    sigma_sqd_w = zt2zm( zm2zt( sigma_sqd_w ) )
    sigma_sqd_w = max( zero_threshold, sigma_sqd_w ) ! Pos. def. quantity

    ! Interpolate the the stats_zt grid
    sigma_sqd_w_zt = max( zm2zt( sigma_sqd_w ), zero_threshold )  ! Pos. def. quantity

    !---------------------------------------------------------------------------
    ! Interpolate thlp2, rtp2, and rtpthlp to thermodynamic levels,
    !---------------------------------------------------------------------------

    ! Interpolate variances to the stats_zt grid (statistics and closure)
    thlp2_zt   = max( zm2zt( thlp2 ), thl_tol**2 ) ! Positive def. quantity
    rtp2_zt    = max( zm2zt( rtp2 ), rt_tol**2 )   ! Positive def. quantity
    rtpthlp_zt = zm2zt( rtpthlp )

    ! Compute skewness velocity for stats output purposes
    if ( iSkw_velocity > 0 ) then
      Skw_velocity = ( 1.0_core_rknd / ( 1.0_core_rknd - sigma_sqd_w(1:gr%nz) ) ) &
                   * ( wp3_zm(1:gr%nz) / max( wp2(1:gr%nz), w_tol_sqd ) )
    end if

    !----------------------------------------------------------------
    ! Call closure scheme
    !----------------------------------------------------------------

    ! Put passive scalar input on the t grid for the PDF
    do i = 1, sclr_dim, 1
      wpsclrp_zt(:,i)   = zm2zt( wpsclrp(:,i) )
      sclrp2_zt(:,i)    = max( zm2zt( sclrp2(:,i) ), zero_threshold ) ! Pos. def. quantity
      sclrprtp_zt(:,i)  = zm2zt( sclrprtp(:,i) )
      sclrpthlp_zt(:,i) = zm2zt( sclrpthlp(:,i) )
    end do ! i = 1, sclr_dim, 1

    ! Interpolate hydrometeor mixed moments to momentum levels.
    do i = 1, hydromet_dim, 1
       wphydrometp_zt(:,i) = zm2zt( wphydrometp(:,i) )
    enddo ! i = 1, hydromet_dim, 1


    call pdf_closure &
         ( hydromet_dim, p_in_Pa, exner, thv_ds_zt,        & ! intent(in)
           wm_zt, wp2_zt, wp3, sigma_sqd_w_zt,             & ! intent(in)
           Skw_zt, Skthl_zt, Skrt_zt,                      & ! intent(in)
           rtm, rtp2_zt, zm2zt( wprtp ),                   & ! intent(in)
           thlm, thlp2_zt, zm2zt( wpthlp ),                & ! intent(in)
           um, zm2zt( up2 ), zm2zt( upwp ),                & ! intent(in)
           vm, zm2zt( vp2 ), zm2zt( vpwp ),                & ! intent(in)
           rtpthlp_zt,                                     & ! intent(in)
           sclrm, wpsclrp_zt, sclrp2_zt,                   & ! intent(in)
           sclrprtp_zt, sclrpthlp_zt,                      & ! intent(in)
#ifdef GFDL
           RH_crit, do_liquid_only_in_clubb,               & ! intent(in)
#endif
           wphydrometp_zt, wp2hmp,                         & ! intent(in)
           rtphmp_zt, thlphmp_zt,                          & ! intent(in)
           wp4_zt, wprtp2, wp2rtp,                         & ! intent(out)
           wpthlp2, wp2thlp, wprtpthlp,                    & ! intent(out)
           cloud_frac, ice_supersat_frac,                  & ! intent(out)
           rcm, wpthvp_zt, wp2thvp, rtpthvp_zt,            & ! intent(out)
           thlpthvp_zt, wprcp_zt, wp2rcp, rtprcp_zt,       & ! intent(out)
           thlprcp_zt, rcp2_zt,                            & ! intent(out)
           uprcp_zt, vprcp_zt,                             & ! intent(out)
           pdf_params, pdf_implicit_coefs_terms,           & ! intent(out)
           F_w, F_rt, F_thl,                               & ! intent(out)
           min_F_w, max_F_w,                               & ! intent(out)
           min_F_rt, max_F_rt,                             & ! intent(out)
           min_F_thl, max_F_thl,                           & ! intent(out)
           wpsclrprtp, wpsclrp2, sclrpthvp_zt,             & ! intent(out)
           wpsclrpthlp, sclrprcp_zt, wp2sclrp,             & ! intent(out)
           rc_coef                                         ) ! intent(out)

    ! Subroutine may produce NaN values, and if so, return
    if ( clubb_at_least_debug_level( 0 ) ) then
       if ( err_code == clubb_fatal_error ) then
          write(fstderr,*) "After pdf_closure"
          return
       endif
    endif

    ! Stats output
    if ( l_stats_samp .and. l_samp_stats_in_pdf_call ) then
       call stat_update_var( iF_w, F_w, stats_zt )
       call stat_update_var( iF_rt, F_rt, stats_zt )
       call stat_update_var( iF_thl, F_thl, stats_zt )
       call stat_update_var( imin_F_w, min_F_w, stats_zt )
       call stat_update_var( imax_F_w, max_F_w, stats_zt )
       call stat_update_var( imin_F_rt, min_F_rt, stats_zt )
       call stat_update_var( imax_F_rt, max_F_rt, stats_zt )
       call stat_update_var( imin_F_thl, min_F_thl, stats_zt )
       call stat_update_var( imax_F_thl, max_F_thl, stats_zt )
    endif

    if ( l_refine_grid_in_cloud ) then

      ! Compute cloud_frac and rcm on a refined grid to improve parameterization
      ! of subgrid clouds
      do k=1, gr%nz

        if ( pdf_params%chi_1(k)/pdf_params%stdev_chi_1(k) > -1._core_rknd ) then

          ! Recalculate cloud_frac and r_c for each PDF component

          call calc_vert_avg_cf_component &
               ( gr%nz, k, gr%zt, pdf_params%chi_1, &                    ! Intent(in)
                 pdf_params%stdev_chi_1, (/(chi_at_liq_sat,i=1,gr%nz)/), & ! Intent(in)
                 cloud_frac_1_refined, rc_1_refined )                   ! Intent(out)

          call calc_vert_avg_cf_component &
               ( gr%nz, k, gr%zt, pdf_params%chi_2, &                     ! Intent(in)
                 pdf_params%stdev_chi_2, (/(chi_at_liq_sat,i=1,gr%nz)/), &  ! Intent(in)
                 cloud_frac_2_refined, rc_2_refined )                    ! Intent(out)

          cloud_frac_refined = compute_mean_binormal &
                               ( cloud_frac_1_refined, cloud_frac_2_refined, &
                                 pdf_params%mixt_frac(k) )

          rcm_refined = compute_mean_binormal &
                        ( rc_1_refined, rc_2_refined, pdf_params%mixt_frac(k) )

          if ( l_interactive_refined ) then
            ! I commented out the lines that modify the values in pdf_params, as it seems that
            ! these values need to remain consistent with the rest of the PDF.
            ! Eric Raut Jun 2014
            ! Replace pdf_closure estimates with refined estimates
            ! pdf_params%rc_1(k) = rc_1_refined
            ! pdf_params%rc_2(k) = rc_2_refined
            rcm(k) = rcm_refined

            ! pdf_params%cloud_frac_1(k) = cloud_frac_1_refined
            ! pdf_params%cloud_frac_2(k) = cloud_frac_2_refined
            cloud_frac(k) = cloud_frac_refined
          end if

        else
          ! Set these equal to the non-refined values so we have something to
          ! output to stats!
          cloud_frac_refined = cloud_frac(k)
          rcm_refined = rcm(k)
        end if ! pdf_params%chi_1(k)/pdf_params%stdev_chi_1(k) > -1._core_rknd

        ! Stats output
        if ( l_stats_samp .and. l_samp_stats_in_pdf_call ) then
          call stat_update_var_pt( icloud_frac_refined, k, cloud_frac_refined, stats_zt )
          call stat_update_var_pt( ircm_refined, k, rcm_refined, stats_zt )
        end if

      end do ! k=1, gr%nz

    end if ! l_refine_grid_in_cloud

    if( l_rtm_nudge ) then
      ! Nudge rtm to prevent excessive drying
      where( rtm < rtm_min .and. gr%zt < rtm_nudge_max_altitude )
        rtm = rtm + (rtm_ref - rtm) * ( dt / ts_nudge )
      end where
    end if


    if ( l_call_pdf_closure_twice ) then

      ! Call pdf_closure a second time on momentum levels, to
      ! output (rather than interpolate) the variables which
      ! belong on the momentum levels.

      ! Interpolate sclrm to the momentum level for use in
      ! the second call to pdf_closure
      do i = 1, sclr_dim
        sclrm_zm(:,i) = zt2zm( sclrm(:,i) )
        ! Clip if extrap. causes sclrm_zm to be less than sclr_tol
        sclrm_zm(gr%nz,i) = max( sclrm_zm(gr%nz,i), sclr_tol(i) )
      end do ! i = 1, sclr_dim

      ! Interpolate pressure, p_in_Pa, to momentum levels.
      ! The pressure at thermodynamic level k = 1 has been set to be the surface
      ! (or model lower boundary) pressure.  Since the surface (or model lower
      ! boundary) is located at momentum level k = 1, the pressure there is
      ! p_sfc, which is p_in_Pa(1).  Thus, p_in_Pa_zm(1) = p_in_Pa(1).
      p_in_Pa_zm(:) = zt2zm( p_in_Pa )
      p_in_Pa_zm(1) = p_in_Pa(1)

      ! Clip pressure if the extrapolation leads to a negative value of pressure
      p_in_Pa_zm(gr%nz) = max( p_in_Pa_zm(gr%nz), 0.5_core_rknd*p_in_Pa(gr%nz) )
      ! Set exner at momentum levels, exner_zm, based on p_in_Pa_zm.
      exner_zm(:) = (p_in_Pa_zm(:)/p0)**kappa

      rtm_zm = zt2zm( rtm )
      ! Clip if extrapolation at the top level causes rtm_zm to be < rt_tol
      rtm_zm(gr%nz) = max( rtm_zm(gr%nz), rt_tol )
      thlm_zm = zt2zm( thlm )
      ! Clip if extrapolation at the top level causes thlm_zm to be < thl_tol
      thlm_zm(gr%nz) = max( thlm_zm(gr%nz), thl_tol )

      ! Interpolate hydrometeor mixed moments to momentum levels.
      do i = 1, hydromet_dim, 1
         rtphmp(:,i)    = zt2zm( rtphmp_zt(:,i) )
         thlphmp(:,i)   = zt2zm( thlphmp_zt(:,i) )
         wp2hmp_zm(:,i) = zt2zm( wp2hmp(:,i) )
      enddo ! i = 1, hydromet_dim, 1

      ! Call pdf_closure to output the variables which belong on the momentum grid.

      call pdf_closure &
           ( hydromet_dim, p_in_Pa_zm, exner_zm, thv_ds_zm,        & ! intent(in)
             wm_zm, wp2, wp3_zm, sigma_sqd_w,                      & ! intent(in)
             Skw_zm, Skthl_zm, Skrt_zm,                            & ! intent(in)
             rtm_zm, rtp2, wprtp,                                  & ! intent(in)
             thlm_zm, thlp2, wpthlp,                               & ! intent(in)
             zt2zm( um ), up2, upwp,                               & ! intent(in)
             zt2zm( vm ), vp2, vpwp,                               & ! intent(in)
             rtpthlp,                                              & ! intent(in)
             sclrm_zm, wpsclrp, sclrp2,                            & ! intent(in)
             sclrprtp, sclrpthlp,                                  & ! intent(in)
#ifdef GFDL
             RH_crit,  do_liquid_only_in_clubb,                    & ! intent(in)
#endif
             wphydrometp, wp2hmp_zm,                               & ! intent(in)
             rtphmp, thlphmp,                                      & ! intent(in)
             wp4, wprtp2_zm, wp2rtp_zm,                            & ! intent(out)
             wpthlp2_zm, wp2thlp_zm, wprtpthlp_zm,                 & ! intent(out)
             cloud_frac_zm, ice_supersat_frac_zm,                  & ! intent(out)
             rcm_zm, wpthvp, wp2thvp_zm, rtpthvp,                  & ! intent(out)
             thlpthvp, wprcp, wp2rcp_zm, rtprcp,                   & ! intent(out)
             thlprcp, rcp2,                                        & ! intent(out)
             uprcp, vprcp,                                         & ! intent(out)
             pdf_params_zm, pdf_implicit_coefs_terms_zm,           & ! intent(out)
             F_w_zm, F_rt_zm, F_thl_zm,                            & ! intent(out)
             min_F_w_zm, max_F_w_zm,                               & ! intent(out)
             min_F_rt_zm, max_F_rt_zm,                             & ! intent(out)
             min_F_thl_zm, max_F_thl_zm,                           & ! intent(out)
             wpsclrprtp_zm, wpsclrp2_zm, sclrpthvp,                & ! intent(out)
             wpsclrpthlp_zm, sclrprcp, wp2sclrp_zm,                & ! intent(out)
             rc_coef_zm                                            ) ! intent(out)

      ! Subroutine may produce NaN values, and if so, return
      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "After pdf_closure"
            return
         endif
      endif

    else ! l_call_pdf_closure_twice is false

      ! Interpolate momentum variables output from the first call to
      ! pdf_closure back to momentum grid.
      ! Since top momentum level is higher than top thermo level,
      ! Set variables at top momentum level to 0.
      if ( l_explicit_turbulent_adv_wp3 .or. iwp4 > 0 ) then
         wp4 = max( zt2zm( wp4_zt ), zero_threshold )  ! Pos. def. quantity
         wp4(gr%nz) = zero
         ! Set wp4 to 0 at the lowest momentum level (momentum level 1).
         ! When the l_explicit_turbulent_adv_wp3 flag is enabled, wp4 (including
         ! its value at momentum level 1) is used interactively in the code.
         ! The value of wp4 at momentum level 1 is found by interpolation of
         ! the values produced by the PDF for wp4_zt at thermodynamic levels
         ! 1 and 2.  This value is unreliable at thermodynamic level 1.
         wp4(1) = zero
      endif

#ifndef CLUBB_CAM
      ! CAM-CLUBB needs cloud water variance thus always compute this
      if ( ircp2 > 0 ) then
#endif
         rcp2 = max( zt2zm( rcp2_zt ), zero_threshold )  ! Pos. def. quantity
#ifndef CLUBB_CAM
         rcp2(gr%nz) = zero
      endif
#endif

      wpthvp            = zt2zm( wpthvp_zt )
      wpthvp(gr%nz)     = 0.0_core_rknd
      thlpthvp          = zt2zm( thlpthvp_zt )
      thlpthvp(gr%nz)   = 0.0_core_rknd
      rtpthvp           = zt2zm( rtpthvp_zt )
      rtpthvp(gr%nz)    = 0.0_core_rknd
      wprcp             = zt2zm( wprcp_zt )
      wprcp(gr%nz)      = 0.0_core_rknd
      rc_coef_zm        = zt2zm( rc_coef )
      rc_coef_zm(gr%nz) = 0.0_core_rknd
      rtprcp            = zt2zm( rtprcp_zt )
      rtprcp(gr%nz)     = 0.0_core_rknd
      thlprcp           = zt2zm( thlprcp_zt )
      thlprcp(gr%nz)    = 0.0_core_rknd
      uprcp             = zt2zm( uprcp_zt )
      uprcp(gr%nz)      = 0.0_core_rknd
      vprcp             = zt2zm( vprcp_zt )
      vprcp(gr%nz)      = 0.0_core_rknd

      ! Initialize variables to avoid uninitialized variables.
      cloud_frac_zm   = 0.0_core_rknd
      ice_supersat_frac_zm = 0.0_core_rknd
      rcm_zm = 0.0_core_rknd
      rtm_zm = 0.0_core_rknd
      thlm_zm = 0.0_core_rknd

      ! Interpolate passive scalars back onto the m grid
      do i = 1, sclr_dim
        sclrpthvp(:,i)       = zt2zm( sclrpthvp_zt(:,i) )
        sclrpthvp(gr%nz,i) = 0.0_core_rknd
        sclrprcp(:,i)        = zt2zm( sclrprcp_zt(:,i) )
        sclrprcp(gr%nz,i)  = 0.0_core_rknd
      end do ! i=1, sclr_dim

    end if ! l_call_pdf_closure_twice

    ! If l_trapezoidal_rule_zt is true, call trapezoidal_rule_zt for
    ! thermodynamic-level variables output from pdf_closure.
    ! ldgrant June 2009
    if ( l_trapezoidal_rule_zt ) then
      call trapezoidal_rule_zt &
           ( l_call_pdf_closure_twice,                    & ! intent(in)
             wprtp2, wpthlp2,                             & ! intent(inout)
             wprtpthlp, cloud_frac, ice_supersat_frac,    & ! intent(inout)
             rcm, wp2thvp, wpsclrprtp, wpsclrp2,          & ! intent(inout)
             wpsclrpthlp, pdf_params,                     & ! intent(inout)
             wprtp2_zm, wpthlp2_zm,                       & ! intent(inout)
             wprtpthlp_zm, cloud_frac_zm,                 & ! intent(inout)
             ice_supersat_frac_zm, rcm_zm, wp2thvp_zm,    & ! intent(inout)
             wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm,  & ! intent(inout)
             pdf_params_zm )                                ! intent(inout)
    end if ! l_trapezoidal_rule_zt

    ! If l_trapezoidal_rule_zm is true, call trapezoidal_rule_zm for
    ! the important momentum-level variabes output from pdf_closure.
    ! ldgrant Feb. 2010
    if ( l_trapezoidal_rule_zm ) then
      call trapezoidal_rule_zm &
         ( wpthvp_zt, thlpthvp_zt, rtpthvp_zt, & ! intent(in)
           wpthvp, thlpthvp, rtpthvp )           ! intent(inout)
    end if ! l_trapezoidal_rule_zm

    ! Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
    ! This code won't work unless rtm >= 0 !!!
    ! We do not clip rcm_in_layer because rcm_in_layer only influences
    ! radiation, and we do not want to bother recomputing it.
    ! Code is duplicated from below to ensure that relative humidity
    ! is calculated properly.  3 Sep 2009
    call clip_rcm( rtm, 'rtm < rcm after pdf_closure', & ! intent (in)
                   rcm )                                 ! intent (inout)

    ! Compute variables cloud_cover and rcm_in_layer.
    ! Added July 2009
    call compute_cloud_cover &
       ( pdf_params, cloud_frac, rcm, & ! intent(in)
         cloud_cover, rcm_in_layer )    ! intent(out)

    ! Use cloud_cover and rcm_in_layer to help boost cloud_frac and rcm to help
    ! increase cloudiness at coarser grid resolutions.
    if ( l_use_cloud_cover ) then
      cloud_frac = cloud_cover
      rcm = rcm_in_layer
    end if

    ! Clip cloud fraction here if it still exceeds 1.0 due to round off
    cloud_frac = min( 1.0_core_rknd, cloud_frac )
    ! Ditto with ice cloud fraction
    ice_supersat_frac = min( 1.0_core_rknd, ice_supersat_frac )

    if (l_use_ice_latent) then
      !A third call to pdf_closure, with terms modified to include the effects
      !of latent heating due to ice.  Thlm and rtm add the effects of ice, and
      !the terms are all renamed with "_frz" appended. The modified terms will
      !be fed into the calculations of the turbulence terms. storer-3/14/13

      !Also added rain for completeness. storer-3/4/14

      if ( iirr > 0 ) then
        rrm = hydromet(:,iirr)
      else
        rrm = zero
      end if

      thlm_frz = thlm - (Lv / (Cp*exner) ) * rrm - (Ls / (Cp*exner) ) * rfrzm
      rtm_frz = rtm + rrm + rfrzm


      call pdf_closure &
           ( hydromet_dim, p_in_Pa, exner, thv_ds_zt,                  & ! intent(in)
             wm_zt, wp2_zt, wp3, sigma_sqd_w_zt,                       & ! intent(in)
             Skw_zt, Skthl_zt, Skrt_zt,                                & ! intent(in)
             rtm_frz, rtp2_zt, zm2zt( wprtp ),                         & ! intent(in)
             thlm_frz, thlp2_zt, zm2zt( wpthlp ),                      & ! intent(in)
             um, zm2zt( up2 ), zm2zt( upwp ),                          & ! intent(in)
             vm, zm2zt( vp2 ), zm2zt( vpwp ),                          & ! intent(in)
             rtpthlp_zt,                                               & ! intent(in)
             sclrm, wpsclrp_zt, sclrp2_zt,                             & ! intent(in)
             sclrprtp_zt, sclrpthlp_zt,                                & ! intent(in)
#ifdef GFDL
             RH_crit, do_liquid_only_in_clubb,                         & ! intent(in)
#endif
             wphydrometp_zt, wp2hmp,                                   & ! intent(in)
             rtphmp_zt, thlphmp_zt,                                    & ! intent(in)
             wp4_zt_frz, wprtp2_frz, wp2rtp_frz,                       & ! intent(out)
             wpthlp2_frz, wp2thlp_frz, wprtpthlp_frz,                  & ! intent(out)
             cloud_frac_frz, ice_supersat_frac_frz,                    & ! intent(out)
             rcm_frz, wpthvp_zt_frz, wp2thvp_frz, rtpthvp_zt_frz,      & ! intent(out)
             thlpthvp_zt_frz, wprcp_zt_frz, wp2rcp_frz, rtprcp_zt_frz, & ! intent(out)
             thlprcp_zt_frz, rcp2_zt_frz,                              & ! intent(out)
             uprcp_zt_frz, vprcp_zt_frz,                               & ! intent(out)
             pdf_params_frz, pdf_implicit_coefs_terms_frz,             & ! intent(out)
             F_w_frz, F_rt_frz, F_thl_frz,                             & ! intent(out)
             min_F_w_frz, max_F_w_frz,                                 & ! intent(out)
             min_F_rt_frz, max_F_rt_frz,                               & ! intent(out)
             min_F_thl_frz, max_F_thl_frz,                             & ! intent(out)
             wpsclrprtp_frz, wpsclrp2_frz, sclrpthvp_zt_frz,           & ! intent(out)
             wpsclrpthlp_frz, sclrprcp_zt_frz, wp2sclrp_frz,           & ! intent(out)
             rc_coef_frz                                               ) ! intent(out)

      ! Subroutine may produce NaN values, and if so, return
      if ( clubb_at_least_debug_level( 0 ) ) then
         if ( err_code == clubb_fatal_error ) then
            write(fstderr,*) "After pdf_closure"
            return
         endif
      endif


      if( l_rtm_nudge ) then
        ! Nudge rtm to prevent excessive drying
        where( rtm < rtm_min .and. gr%zt < rtm_nudge_max_altitude )
          rtm = rtm + (rtm_ref - rtm) * ( dt / ts_nudge )
        end where
      end if

      rtm_zm_frz = zt2zm( rtm_frz )
      ! Clip if extrapolation at the top level causes rtm_zm to be < rt_tol
      rtm_zm_frz(gr%nz) = max( rtm_zm_frz(gr%nz), rt_tol )
      thlm_zm_frz = zt2zm( thlm_frz )
      ! Clip if extrapolation at the top level causes thlm_zm to be < thl_tol
      thlm_zm_frz(gr%nz) = max( thlm_zm_frz(gr%nz), thl_tol )

      if ( l_call_pdf_closure_twice ) then

        call init_pdf_params( gr%nz, pdf_params_zm_frz )

        ! Call pdf_closure again to output the variables which belong on the momentum grid.
        call pdf_closure &
             ( hydromet_dim, p_in_Pa_zm, exner_zm, thv_ds_zm,         & ! intent(in)
               wm_zm, wp2, wp3_zm, sigma_sqd_w,                       & ! intent(in)
               Skw_zm, Skthl_zm, Skrt_zm,                             & ! intent(in)
               rtm_zm_frz, rtp2, wprtp,                               & ! intent(in)
               thlm_zm_frz, thlp2, wpthlp,                            & ! intent(in)
               zt2zm( um ), up2, upwp,                                & ! intent(in)
               zt2zm( vm ), vp2, vpwp,                                & ! intent(in)
               rtpthlp,                                               & ! intent(in)
               sclrm_zm, wpsclrp, sclrp2,                             & ! intent(in)
               sclrprtp, sclrpthlp,                                   & ! intent(in)
#ifdef GFDL
               RH_crit, do_liquid_only_in_clubb,                      & ! intent(in)
#endif
               wphydrometp, wp2hmp_zm,                                & ! intent(in)
               rtphmp, thlphmp,                                       & ! intent(in)
               wp4_frz, wprtp2_zm_frz, wp2rtp_zm_frz,                 & ! intent(out)
               wpthlp2_zm_frz, wp2thlp_zm_frz, wprtpthlp_zm_frz,      & ! intent(out)
               cloud_frac_zm_frz, ice_supersat_frac_zm_frz,           & ! intent(out)
               rcm_zm_frz, wpthvp_frz, wp2thvp_zm_frz, rtpthvp_frz,   & ! intent(out)
               thlpthvp_frz, wprcp_frz, wp2rcp_zm_frz, rtprcp_frz,    & ! intent(out)
               thlprcp_frz, rcp2_frz,                                 & ! intent(out)
               uprcp_frz, vprcp_frz,                                  & ! intent(out)
               pdf_params_zm_frz, pdf_implicit_coefs_terms_zm_frz,    & ! intent(out)
               F_w_zm_frz, F_rt_zm_frz, F_thl_zm_frz,                 & ! intent(out)
               min_F_w_zm_frz, max_F_w_zm_frz,                        & ! intent(out)
               min_F_rt_zm_frz, max_F_rt_zm_frz,                      & ! intent(out)
               min_F_thl_zm_frz, max_F_thl_zm_frz,                    & ! intent(out)
               wpsclrprtp_zm_frz, wpsclrp2_zm_frz, sclrpthvp_frz,     & ! intent(out)
               wpsclrpthlp_zm_frz, sclrprcp_frz, wp2sclrp_zm_frz,     & ! intent(out)
               rc_coef_zm_frz                                         ) ! intent(out)

        ! Subroutine may produce NaN values, and if so, return
        if ( clubb_at_least_debug_level( 0 ) ) then
           if ( err_code == clubb_fatal_error ) then
              write(fstderr,*) "After pdf_closure"
              return
           endif
        endif

      else ! l_call_pdf_closure_twice is false

        wpthvp_frz            = zt2zm( wpthvp_zt_frz )
        wpthvp_frz(gr%nz)   = 0.0_core_rknd
        thlpthvp_frz          = zt2zm( thlpthvp_zt_frz )
        thlpthvp_frz(gr%nz) = 0.0_core_rknd
        rtpthvp_frz           = zt2zm( rtpthvp_zt_frz )
        rtpthvp_frz(gr%nz)  = 0.0_core_rknd

      end if ! l_call_pdf_closure_twice

      if ( l_trapezoidal_rule_zt ) then
        call trapezoidal_rule_zt &
           ( l_call_pdf_closure_twice,                                & ! intent(in)
             wprtp2_frz, wpthlp2_frz,                                 & ! intent(inout)
             wprtpthlp_frz, cloud_frac_frz, ice_supersat_frac_frz,    & ! intent(inout)
             rcm_frz, wp2thvp_frz, wpsclrprtp_frz, wpsclrp2_frz,      & ! intent(inout)
             wpsclrpthlp_frz, pdf_params_frz,                         & ! intent(inout)
             wprtp2_zm_frz, wpthlp2_zm_frz,                           & ! intent(inout)
             wprtpthlp_zm_frz, cloud_frac_zm_frz,                     & ! intent(inout)
             ice_supersat_frac_zm_frz, rcm_zm_frz, wp2thvp_zm_frz,    & ! intent(inout)
             wpsclrprtp_zm_frz, wpsclrp2_zm_frz, wpsclrpthlp_zm_frz,  & ! intent(inout)
             pdf_params_zm_frz                                        ) ! intent(inout)
      end if ! l_trapezoidal_rule_zt

        ! If l_trapezoidal_rule_zm is true, call trapezoidal_rule_zm for
        ! the important momentum-level variabes output from pdf_closure.
        ! ldgrant Feb. 2010
        if ( l_trapezoidal_rule_zm ) then
          call trapezoidal_rule_zm &
             ( wpthvp_zt_frz, thlpthvp_zt_frz, rtpthvp_zt_frz, & ! intent(in)
               wpthvp_frz, thlpthvp_frz, rtpthvp_frz )           ! intent(inout)
        end if ! l_trapezoidal_rule_zm

        wpthvp = wpthvp_frz
        wp2thvp = wp2thvp_frz
        thlpthvp = thlpthvp_frz
        rtpthvp = rtpthvp_frz

      end if ! l_use_ice_latent = .true.

      rsat = sat_mixrat_liq( p_in_Pa, thlm2T_in_K( thlm, exner, rcm ) )
      rel_humidity = (rtm - rcm) / rsat

      rcm_supersat_adj = zero
      if ( l_rcm_supersat_adj ) then
        ! +PAB mods, take remaining supersaturation that may exist
        !   after CLUBB PDF call and add it to rcm.  Supersaturation
        !   may exist after PDF call due to issues with calling PDF on the
        !   thermo grid and momentum grid and the interpolation between the two
        l_spur_supersat = .false.
        do k = 2, gr%nz
          if (rel_humidity(k) > 1.0_core_rknd) then
            rcm_supersat_adj(k) = (rtm(k) - rcm(k)) - rsat(k)
            rcm(k) = rcm(k) + rcm_supersat_adj(k)
            l_spur_supersat = .true.
          end if
        enddo

        if ( clubb_at_least_debug_level( 1 ) .and. l_spur_supersat ) then
          write(fstderr,*) 'Warning: spurious supersaturation was removed after pdf_closure!'
        end if

      end if ! l_rcm_supersat_adj


    return

  end subroutine pdf_closure_driver

  !=============================================================================
    subroutine setup_clubb_core &
               ( nzmax, T0_in, ts_nudge_in,               & ! intent(in)
                 hydromet_dim_in, sclr_dim_in,            & ! intent(in)
                 sclr_tol_in, edsclr_dim_in, params,      & ! intent(in)
                 l_host_applies_sfc_fluxes,               & ! intent(in)
                 l_uv_nudge, saturation_formula,          & ! intent(in)
                 l_input_fields,                          & ! intent(in)
#ifdef GFDL
                 I_sat_sphum,                             & ! intent(in)  h1g, 2010-06-16
#endif
                 l_implemented, grid_type, deltaz,        & ! intent(in)
                 zm_init, zm_top,                         & ! intent(in)
                 momentum_heights, thermodynamic_heights, & ! intent(in)
                 sfc_elevation                            & ! intent(in)
#ifdef GFDL
                 , cloud_frac_min                         & ! intent(in)  h1g, 2010-06-16
#endif
                )

      ! Description:
      !   Subroutine to set up the model for execution.
      !
      ! References:
      !   None
      !---------------------------------------------------------------------

      use grid_class, only: &
          setup_grid, & ! Procedure
          gr ! Variable(s)

      use parameter_indices, only:  &
          nparams ! Variable(s)

      use parameters_tunable, only: &
          setup_parameters ! Procedure

      use parameters_model, only: &
          setup_parameters_model ! Procedure

      use variables_diagnostic_module, only: &
          setup_diagnostic_variables ! Procedure

      use variables_prognostic_module, only: &
          setup_prognostic_variables ! Procedure

      use constants_clubb, only:  &
          fstderr  ! Variable(s)

      use error_code, only: &
          clubb_at_least_debug_level,  & ! Procedures
          initialize_error_headers,    &
          err_code,                    & ! Error Indicator
          clubb_fatal_error              ! Constant

      use model_flags, only: &
          setup_model_flags, & ! Subroutine
          l_use_ice_latent, & ! Variable(s)
          l_predict_upwp_vpwp, &
          l_explicit_turbulent_adv_wpxp

      use pdf_closure_module, only: &
          iiPDF_ADG1,     & ! Variable(s)
          iiPDF_ADG2,     &
          iiPDF_3D_Luhar, &
          iiPDF_new,      &
          iiPDF_TSDADG,   &
          iiPDF_LY93,     &
          iiPDF_type

#ifdef MKL
      use csr_matrix_module, only: &
          initialize_csr_matrix, & ! Subroutine
          intlc_5d_5d_ja_size      ! Variable

      use gmres_wrap, only: &
          gmres_init              ! Subroutine

      use gmres_cache, only: &
          gmres_cache_temp_init, & ! Subroutine
          gmres_idx_wp2wp3         ! Variable
#endif /* MKL */

      use clubb_precision, only: &
          core_rknd ! Variable(s)

      implicit none

      ! Input Variables

      ! Grid definition
      integer, intent(in) :: nzmax  ! Vertical grid levels            [#]
      !                      Only true when used in a host model
      !                      CLUBB determines what nzmax should be
      !                      given zm_init and zm_top when
      !                      running in standalone mode.

      real( kind = core_rknd ), intent(in) ::  &
        sfc_elevation  ! Elevation of ground level    [m AMSL]

      ! Flag to see if CLUBB is running on it's own,
      ! or if it's implemented as part of a host model.
      logical, intent(in) :: l_implemented   ! (T/F)

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

      logical, intent(in) ::  &
        l_input_fields    ! Flag for whether LES input fields are being used

#ifdef GFDL
      logical, intent(in) :: &  ! h1g, 2010-06-16 begin mod
         I_sat_sphum

      real( kind = core_rknd ), intent(in) :: &
         cloud_frac_min         ! h1g, 2010-06-16 end mod
#endif

      ! Local variables
      integer :: begin_height, end_height

      !----- Begin Code -----

      call initialize_error_headers

      ! Sanity check for the saturation formula
      select case ( trim( saturation_formula ) )
      case ( "bolton", "Bolton" )
        ! Using the Bolton 1980 approximations for SVP over vapor/ice

      case ( "flatau", "Flatau" )
        ! Using the Flatau, et al. polynomial approximation for SVP over vapor/ice

      case ( "gfdl", "GFDL" )   ! h1g, 2010-06-16
        ! Using the GFDL SVP formula (Goff-Gratch)

        ! Add new saturation formulas after this

      case default
        write(fstderr,*) "Error in setup_clubb_core."
        write(fstderr,*) "Unknown approx. of saturation vapor pressure: "// &
          trim( saturation_formula )
        err_code = clubb_fatal_error
        return
      end select

      ! Check for the type of two component normal (double Gaussian) PDF being
      ! used for w, rt, and theta-l (or w, chi, and eta).
      if ( iiPDF_type < iiPDF_ADG1 .or. iiPDF_type > iiPDF_LY93 ) then
         write(fstderr,*) "Error in setup_clubb_core."
         write(fstderr,*) "Unknown type of double Gaussian PDF selected."
         write(fstderr,*) "iiPDF_type = ", iiPDF_type
         err_code = clubb_fatal_error
         return
      endif ! iiPDF_type < iiPDF_ADG1 or iiPDF_type > iiPDF_lY93

      ! The ADG2 and 3D Luhar PDFs can only be used as part of input fields.
      if ( iiPDF_type == iiPDF_ADG2 ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "Error in setup_clubb_core."
            write(fstderr,*) "The ADG2 PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_ADG2

      if ( iiPDF_type == iiPDF_3D_Luhar ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "Error in setup_clubb_core."
            write(fstderr,*) "The 3D Luhar PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_3D_Luhar

      ! This also currently applies to the new PDF until it has been fully
      ! implemented.
      if ( iiPDF_type == iiPDF_new ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "Error in setup_clubb_core."
            write(fstderr,*) "The new PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_new

      ! This also currently applies to the TSDADG PDF until it has been fully
      ! implemented.
      if ( iiPDF_type == iiPDF_TSDADG ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "Error in setup_clubb_core."
            write(fstderr,*) "The new TSDADG PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_TSDADG

      ! This also applies to Lewellen and Yoh (1993).
      if ( iiPDF_type == iiPDF_LY93 ) then
         if ( .not. l_input_fields ) then
            write(fstderr,*) "Error in setup_clubb_core."
            write(fstderr,*) "The Lewellen and Yoh PDF can only be used with" &
                             // " input fields (l_input_fields = .true.)."
            write(fstderr,*) "iiPDF_type = ", iiPDF_type
            write(fstderr,*) "l_input_fields = ", l_input_fields
            err_code = clubb_fatal_error
            return
         endif ! .not. l_input_fields
      endif ! iiPDF_type == iiPDF_LY93

      ! Check the option for the placement of the call to CLUBB's PDF.
      if ( ipdf_call_placement < ipdf_pre_advance_fields &
           .or. ipdf_call_placement > ipdf_pre_post_advance_fields ) then
         write(fstderr,*) "Invalid option selected for ipdf_call_placement"
         err_code = clubb_fatal_error
         return
      endif

      ! When ipdf_call_placement = ipdf_post_advance_fields, additional
      ! variables need to be passed out of advance_clubb_core, saved, and passed
      ! in again during the next model timestep.  In order to keep from needing
      ! to pass out and save variables that are not normally used in CLUBB, some
      ! options (that are typically turned off anyway) will be disallowed when
      ! ipdf_call_placement = ipdf_post_advance_fields.
      if ( ipdf_call_placement == ipdf_post_advance_fields ) then
         if ( l_use_ice_latent ) then
            write(fstderr,*) "The l_use_ice_latent option is incompatible" &
                             // " with the ipdf_post_advance_fields option."
            err_code = clubb_fatal_error
            return
         endif ! l_use_ice_latent
      endif ! ipdf_call_placement == ipdf_post_advance_fields

      ! The l_predict_upwp_vpwp flag requires that the ADG1 PDF is used
      ! implicitly in subroutine advance_xm_wpxp.
      if ( l_predict_upwp_vpwp ) then

         ! When l_predict_upwp_vpwp is enabled, the
         ! l_explicit_turbulent_adv_wpxp flag must be turned off.
         ! Otherwise, explicit turbulent advection would require PDF parameters
         ! for u and v to be calculated in PDF closure.  These would be needed
         ! to calculate integrated fields such as wp2up, etc.
         if ( l_explicit_turbulent_adv_wpxp ) then
            write(fstderr,*) "The l_explicit_turbulent_adv_wpxp option" &
                             // " is not currently set up for use with the" &
                             // " l_predict_upwp_vpwp code."
            err_code = clubb_fatal_error
            return
         endif ! l_explicit_turbulent_adv_wpxp

         ! When l_predict_upwp_vpwp is enabled, the PDF type must be set to
         ! the ADG1 PDF.  The other PDFs are not currently set up to calculate
         ! variables needed for implicit or semi-implicit turbulent advection,
         ! such as coef_wp2up_implicit, etc.
         if ( iiPDF_type /= iiPDF_ADG1 ) then
            write(fstderr,*) "Currently, only the ADG1 PDF is set up for use" &
                             // " with the l_predict_upwp_vpwp code."
            err_code = clubb_fatal_error
            return
         endif ! iiPDF_type /= iiPDF_ADG1

         if (  ipdf_call_placement == ipdf_post_advance_fields ) then
            write(fstderr,*) "Currently,  ipdf_call_placement == ipdf_post_advance_fields" &
                             // " is incompatible with the l_predict_upwp_vpwp code" &
                             // " because uprcp has not been fed through advance_clubb_core."
            err_code = clubb_fatal_error
            return
         endif

      endif ! l_predict_upwp_vpwp

      ! Setup grid
      call setup_grid( nzmax, sfc_elevation, l_implemented,     & ! intent(in)
                       grid_type, deltaz, zm_init, zm_top,      & ! intent(in)
                       momentum_heights, thermodynamic_heights, & ! intent(in)
                       begin_height, end_height                 ) ! intent(out)

      if ( clubb_at_least_debug_level( 0 ) ) then
        if ( err_code == clubb_fatal_error ) then

          write(fstderr,*) "Error in setup_clubb_core"

          write(fstderr,*) "Intent(in)"

          write(fstderr,*) "deltaz = ", deltaz
          write(fstderr,*) "zm_init = ", zm_init
          write(fstderr,*) "zm_top = ", zm_top
          write(fstderr,*) "momentum_heights = ", momentum_heights
          write(fstderr,*) "thermodynamic_heights = ",  &
              thermodynamic_heights
          write(fstderr,*) "T0_in = ", T0_in
          write(fstderr,*) "ts_nudge_in = ", ts_nudge_in
          write(fstderr,*) "params = ", params
          return

        end if
    end if

      ! Setup flags
#ifdef GFDL
      call setup_model_flags &
           ( l_host_applies_sfc_fluxes,      & ! intent(in)
             l_uv_nudge, saturation_formula, & ! intent(in)
             I_sat_sphum )                     ! intent(in)  h1g, 2010-06-16

#else
      call setup_model_flags &
           ( l_host_applies_sfc_fluxes,      & ! intent(in)
             l_uv_nudge, saturation_formula )  ! intent(in)
#endif


      ! Define model constant parameters
#ifdef GFDL
      call setup_parameters_model( T0_in, ts_nudge_in,                         & ! intent(in)
                                   hydromet_dim_in,                            & ! intent(in)
                                   sclr_dim_in, sclr_tol_in, edsclr_dim_in,    & ! intent(in)
                                   cloud_frac_min )                 ! intent(in)  h1g, 2010-06-16
#else
      call setup_parameters_model( T0_in, ts_nudge_in,                       & ! intent(in)
                                   hydromet_dim_in,                          & ! intent(in)
                                   sclr_dim_in, sclr_tol_in, edsclr_dim_in )   ! intent(in)
#endif

      ! Define tunable constant parameters
      call setup_parameters &
           ( deltaz, params, gr%nz,                                & ! intent(in)
             grid_type, momentum_heights(begin_height:end_height), & ! intent(in)
             thermodynamic_heights(begin_height:end_height) )        ! intent(in)

      if ( clubb_at_least_debug_level( 0 ) ) then
          if ( err_code == clubb_fatal_error ) then

            write(fstderr,*) "Error in setup_clubb_core"

            write(fstderr,*) "Intent(in)"

            write(fstderr,*) "deltaz = ", deltaz
            write(fstderr,*) "zm_init = ", zm_init
            write(fstderr,*) "zm_top = ", zm_top
            write(fstderr,*) "momentum_heights = ", momentum_heights
            write(fstderr,*) "thermodynamic_heights = ",  &
              thermodynamic_heights
            write(fstderr,*) "T0_in = ", T0_in
            write(fstderr,*) "ts_nudge_in = ", ts_nudge_in
            write(fstderr,*) "params = ", params

            return

          end if
      end if

#ifdef GFDL
! setup  prognostic_variables
      call setup_prognostic_variables( gr%nz ) ! intent(in)  h1g, 2010-06-16
#else
      if ( .not. l_implemented ) then
        call setup_prognostic_variables( gr%nz ) ! intent(in)
      end if
#endif

      ! The diagnostic variables need to be
      ! declared, allocated, initialized, and deallocated whether CLUBB
      ! is part of a larger model or not.
      call setup_diagnostic_variables( gr%nz )  ! intent(in)

#ifdef MKL
      ! Initialize the CSR matrix class.
      if ( l_gmres ) then
        call initialize_csr_matrix
      end if

      if ( l_gmres ) then
        call gmres_cache_temp_init( gr%nz ) ! intent(in)
        call gmres_init( (2 * gr%nz), intlc_5d_5d_ja_size ) ! intent(in)
      end if
#endif /* MKL */

      return
    end subroutine setup_clubb_core

    !----------------------------------------------------------------------------
    subroutine cleanup_clubb_core( l_implemented )
      !
      ! Description:
      !   Frees memory used by the model itself.
      !
      ! References:
      !   None
      !---------------------------------------------------------------------------
      use parameters_model, only: sclr_tol ! Variable

      use variables_diagnostic_module, only: &
        cleanup_diagnostic_variables ! Procedure

      use variables_prognostic_module, only: &
        cleanup_prognostic_variables ! Procedure

      use grid_class, only: &
        cleanup_grid ! Procedure

      use parameters_tunable, only: &
        cleanup_nu ! Procedure

      implicit none

      ! Flag to see if CLUBB is running on it's own,
      ! or if it's implemented as part of a host model.
      logical, intent(in) :: l_implemented   ! (T/F)

      !----- Begin Code -----
#ifdef GFDL
      ! cleanup  prognostic_variables
      call  cleanup_prognostic_variables( )  ! h1g, 2010-06-16
#else
      if ( .not. l_implemented ) then
        call cleanup_prognostic_variables( )
      end if
#endif

      ! The diagnostic variables need to be
      ! declared, allocated, initialized, and deallocated whether CLUBB
      ! is part of a larger model or not.
      call cleanup_diagnostic_variables( )

      ! De-allocate the array for the passive scalar tolerances
      deallocate( sclr_tol )

      ! De-allocate the arrays for the grid
      call cleanup_grid( )

      ! De-allocate the arrays for nu
      call cleanup_nu( )

      return
    end subroutine cleanup_clubb_core

    !-----------------------------------------------------------------------
    subroutine trapezoidal_rule_zt &
               ( l_call_pdf_closure_twice,                    & ! intent(in)
                 wprtp2, wpthlp2,                             & ! intent(inout)
                 wprtpthlp, cloud_frac, ice_supersat_frac,    & ! intent(inout)
                 rcm, wp2thvp, wpsclrprtp, wpsclrp2,          & ! intent(inout)
                 wpsclrpthlp, pdf_params,                     & ! intent(inout)
                 wprtp2_zm, wpthlp2_zm,                       & ! intent(inout)
                 wprtpthlp_zm, cloud_frac_zm,                 & ! intent(inout)
                 ice_supersat_frac_zm, rcm_zm, wp2thvp_zm,    & ! intent(inout)
                 wpsclrprtp_zm, wpsclrp2_zm, wpsclrpthlp_zm,  & ! intent(inout)
                 pdf_params_zm )                                ! intent(inout)
      !
      ! Description:
      !   This subroutine takes the output variables on the thermo.
      !   grid and either: interpolates them to the momentum grid, or uses the
      !   values output from the second call to pdf_closure on momentum levels if
      !   l_call_pdf_closure_twice is true.  It then calls the function
      !   trapezoid_zt to recompute the variables on the thermo. grid.
      !
      !   ldgrant June 2009
      !
      ! Note:
      !   The argument variables in the last 5 lines of the subroutine
      !   (wprtp2_zm through pdf_params_zm) are declared intent(inout) because
      !   if l_call_pdf_closure_twice is true, these variables will already have
      !   values from pdf_closure on momentum levels and will not be altered in
      !   this subroutine.  However, if l_call_pdf_closure_twice is false, these
      !   variables will not have values yet and will be interpolated to
      !   momentum levels in this subroutine.
      ! References:
      !   None
      !-----------------------------------------------------------------------

      use constants_clubb, only: &
          fstderr  ! Constant(s)

      use stats_variables, only: &
          iwprtp2, & ! Varibles
          iwprtpthlp, &
          iwpthlp2, &
          iwprtp2, &
          iwpsclrp2, &
          iwpsclrprtp, &
          iwpsclrpthlp, &
          l_stats

      use grid_class, only: &
          gr, & ! Variable
          zt2zm ! Procedure

      use parameters_model, only: &
          sclr_dim ! Number of passive scalar variables

      use pdf_parameter_module, only: &
          pdf_parameter ! Derived data type

      use clubb_precision, only: &
          core_rknd ! Variable(s)

      implicit none

      ! Constant parameters
      logical, parameter :: &
        l_apply_rule_to_pdf_params = .false. ! Apply the trapezoidal rule to pdf_params

      ! Input variables
      logical, intent(in) :: l_call_pdf_closure_twice

      ! Input/Output variables
      ! Thermodynamic level variables output from the first call to pdf_closure
      real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
        wprtp2,             & ! w'rt'^2                   [m kg^2/kg^2]
        wpthlp2,            & ! w'thl'^2                  [m K^2/s]
        wprtpthlp,          & ! w'rt'thl'                 [m kg K/kg s]
        cloud_frac,         & ! Cloud Fraction            [-]
        ice_supersat_frac,  & ! Ice Cloud Fraction        [-]
        rcm,                & ! Liquid water mixing ratio [kg/kg]
        wp2thvp               ! w'^2 th_v'                [m^2 K/s^2]

      real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(inout) :: &
        wpsclrprtp,  & ! w'sclr'rt'
        wpsclrp2,    & ! w'sclr'^2
        wpsclrpthlp    ! w'sclr'thl'

      type (pdf_parameter), intent(inout) :: &
        pdf_params ! PDF parameters [units vary]

      ! Thermo. level variables brought to momentum levels either by
      ! interpolation (in subroutine trapezoidal_rule_zt) or by
      ! the second call to pdf_closure (in subroutine advance_clubb_core)
      real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
        wprtp2_zm,            & ! w'rt'^2 on momentum grid                   [m kg^2/kg^2]
        wpthlp2_zm,           & ! w'thl'^2 on momentum grid                  [m K^2/s]
        wprtpthlp_zm,         & ! w'rt'thl' on momentum grid                 [m kg K/kg s]
        cloud_frac_zm,        & ! Cloud Fraction on momentum grid            [-]
        ice_supersat_frac_zm, & ! Ice Cloud Fraction on momentum grid        [-]
        rcm_zm,               & ! Liquid water mixing ratio on momentum grid [kg/kg]
        wp2thvp_zm              ! w'^2 th_v' on momentum grid                [m^2 K/s^2]

      real( kind = core_rknd ), dimension(gr%nz,sclr_dim), intent(inout) :: &
        wpsclrprtp_zm,  & ! w'sclr'rt' on momentum grid
        wpsclrp2_zm,    & ! w'sclr'^2 on momentum grid
        wpsclrpthlp_zm    ! w'sclr'thl' on momentum grid

      type (pdf_parameter), intent(inout) :: &
        pdf_params_zm ! PDF parameters on momentum grid [units vary]

      ! Local variables

      ! Components of PDF_parameters on the momentum grid (_zm) and on the thermo. grid (_zt)
      real( kind = core_rknd ), dimension(gr%nz) :: &
        w_1_zt,          & ! Mean of w for 1st normal distribution                 [m/s]
        w_1_zm,          & ! Mean of w for 1st normal distribution                 [m/s]
        w_2_zm,          & ! Mean of w for 2nd normal distribution                 [m/s]
        w_2_zt,          & ! Mean of w for 2nd normal distribution                 [m/s]
        varnce_w_1_zm,   & ! Variance of w for 1st normal distribution         [m^2/s^2]
        varnce_w_1_zt,   & ! Variance of w for 1st normal distribution         [m^2/s^2]
        varnce_w_2_zm,   & ! Variance of w for 2nd normal distribution         [m^2/s^2]
        varnce_w_2_zt,   & ! Variance of w for 2nd normal distribution         [m^2/s^2]
        rt_1_zm,         & ! Mean of r_t for 1st normal distribution             [kg/kg]
        rt_1_zt,         & ! Mean of r_t for 1st normal distribution             [kg/kg]
        rt_2_zm,         & ! Mean of r_t for 2nd normal distribution             [kg/kg]
        rt_2_zt,         & ! Mean of r_t for 2nd normal distribution             [kg/kg]
        varnce_rt_1_zm,  & ! Variance of r_t for 1st normal distribution     [kg^2/kg^2]
        varnce_rt_1_zt,  & ! Variance of r_t for 1st normal distribution     [kg^2/kg^2]
        varnce_rt_2_zm,  & ! Variance of r_t for 2nd normal distribution     [kg^2/kg^2]
        varnce_rt_2_zt,  & ! Variance of r_t for 2nd normal distribution     [kg^2/kg^2]
        crt_1_zm,        & ! Coefficient for s'                                      [-]
        crt_1_zt,        & ! Coefficient for s'                                      [-]
        crt_2_zm           ! Coefficient for s'                                      [-]

      real( kind = core_rknd ), dimension(gr%nz) :: &
        crt_2_zt,        & ! Coefficient for s'                                      [-]
        cthl_1_zm,       & ! Coefficient for s'                                    [1/K]
        cthl_1_zt,       & ! Coefficient for s'                                    [1/K]
        cthl_2_zm,       & ! Coefficient for s'                                    [1/K]
        cthl_2_zt,       & ! Coefficient for s'                                    [1/K]
        thl_1_zm,        & ! Mean of th_l for 1st normal distribution                [K]
        thl_1_zt,        & ! Mean of th_l for 1st normal distribution                [K]
        thl_2_zm,        & ! Mean of th_l for 2nd normal distribution                [K]
        thl_2_zt,        & ! Mean of th_l for 2nd normal distribution
        varnce_thl_1_zm, & ! Variance of th_l for 1st normal distribution          [K^2]
        varnce_thl_1_zt, & ! Variance of th_l for 1st normal distribution          [K^2]
        varnce_thl_2_zm, & ! Variance of th_l for 2nd normal distribution          [K^2]
        varnce_thl_2_zt    ! Variance of th_l for 2nd normal distribution          [K^2]

      real( kind = core_rknd ), dimension(gr%nz) :: &
        mixt_frac_zm,   & ! Weight of 1st normal distribution (Sk_w dependent)      [-]
        mixt_frac_zt,   & ! Weight of 1st normal distribution (Sk_w dependent)      [-]
        rc_1_zm,         & ! Mean of r_c for 1st normal distribution             [kg/kg]
        rc_1_zt,         & ! Mean of r_c for 1st normal distribution             [kg/kg]
        rc_2_zm,         & ! Mean of r_c for 2nd normal distribution             [kg/kg]
        rc_2_zt,         & ! Mean of r_c for 2nd normal distribution             [kg/kg]
        rsatl_1_zm,        & ! Mean of r_sl for 1st normal distribution            [kg/kg]
        rsatl_1_zt,        & ! Mean of r_sl for 1st normal distribution            [kg/kg]
        rsatl_2_zm,        & ! Mean of r_sl for 2nd normal distribution            [kg/kg]
        rsatl_2_zt,        & ! Mean of r_sl for 2nd normal distribution            [kg/kg]
        cloud_frac_1_zm, & ! Cloud fraction for 1st normal distribution              [-]
        cloud_frac_1_zt, & ! Cloud fraction for 1st normal distribution              [-]
        cloud_frac_2_zm, & ! Cloud fraction for 2nd normal distribution              [-]
        cloud_frac_2_zt, & ! Cloud fraction for 2nd normal distribution              [-]
        chi_1_zm,          & ! Mean of chi(s) for 1st normal distribution               [kg/kg]
        chi_1_zt,          & ! Mean of chi(s) for 1st normal distribution               [kg/kg]
        chi_2_zm,          & ! Mean of chi(s) for 2nd normal distribution               [kg/kg]
        chi_2_zt,          & ! Mean of chi(s) for 2nd normal distribution               [kg/kg]
        stdev_chi_1_zm       ! Standard deviation of chi(s) for 1st normal distribution [kg/kg]

      real( kind = core_rknd ), dimension(gr%nz) :: &
        stdev_chi_1_zt,    & ! Standard deviation of chi(s) for 1st normal distribution [kg/kg]
        stdev_chi_2_zm,    & ! Standard deviation of chi(s) for 2nd normal distribution [kg/kg]
        stdev_chi_2_zt,    & ! Standard deviation of chi(s) for 2nd normal distribution [kg/kg]
        stdev_eta_1_zm,    & ! Standard deviation of eta(t) for 1st normal distribution [kg/kg]
        stdev_eta_1_zt,    & ! Standard deviation of eta(t) for 1st normal distribution [kg/kg]
        stdev_eta_2_zm,    & ! Standard deviation of eta(t) for 2nd normal distribution [kg/kg]
        stdev_eta_2_zt,    & ! Standard deviation of eta(t) for 2nd normal distribution [kg/kg]
        corr_w_rt_1_zm,    & ! PDF comp. correlation of w and r_t for 1st normal        [-]
        corr_w_rt_1_zt,    & ! PDF comp. correlation of w and r_t for 1st normal        [-]
        corr_w_rt_2_zm,    & ! PDF comp. correlation of w and r_t for 2nd normal        [-]
        corr_w_rt_2_zt,    & ! PDF comp. correlation of w and r_t for 2nd normal        [-]
        corr_w_thl_1_zm,   & ! PDF comp. correlation of w and th_l for 1st normal       [-]
        corr_w_thl_1_zt,   & ! PDF comp. correlation of w and th_l for 1st normal       [-]
        corr_w_thl_2_zm,   & ! PDF comp. correlation of w and th_l for 2nd normal       [-]
        corr_w_thl_2_zt,   & ! PDF comp. correlation of w and th_l for 2nd normal       [-]
        corr_rt_thl_1_zm,  & ! PDF comp. correlation of r_t and th_l for 1st normal     [-]
        corr_rt_thl_1_zt,  & ! PDF comp. correlation of r_t and th_l for 1st normal     [-]
        corr_rt_thl_2_zm,  & ! PDF comp. correlation of r_t and th_l for 2nd normal     [-]
        corr_rt_thl_2_zt,  & ! PDF comp. correlation of r_t and th_l for 2nd normal     [-]
        alpha_thl_zm,      & ! Factor relating to normalized variance for th_l          [-]
        alpha_thl_zt,      & ! Factor relating to normalized variance for th_l          [-]
        alpha_rt_zm,       & ! Factor relating to normalized variance for r_t           [-]
        alpha_rt_zt          ! Factor relating to normalized variance for r_t           [-]

      integer :: i

      !----------------------- Begin Code -----------------------------

      ! Store components of pdf_params in the locally declared variables
      ! We only apply the trapezoidal rule to these when
      ! l_apply_rule_to_pdf_params is true.  This is because when we apply the
      ! rule to the final result of pdf_closure rather than the intermediate
      ! results it can lead to an inconsistency in how we determine which
      ! PDF component a point is in and whether the point is in or out of cloud,
      ! which is turn will break the latin hypercube code that samples
      ! preferentially in cloud. -dschanen 13 Feb 2012

      if ( l_apply_rule_to_pdf_params ) then
        w_1_zt          = pdf_params%w_1
        w_2_zt          = pdf_params%w_2
        varnce_w_1_zt   = pdf_params%varnce_w_1
        varnce_w_2_zt   = pdf_params%varnce_w_2
        rt_1_zt         = pdf_params%rt_1
        rt_2_zt         = pdf_params%rt_2
        varnce_rt_1_zt  = pdf_params%varnce_rt_1
        varnce_rt_2_zt  = pdf_params%varnce_rt_2
        crt_1_zt        = pdf_params%crt_1
        crt_2_zt        = pdf_params%crt_2
        cthl_1_zt       = pdf_params%cthl_1
        cthl_2_zt       = pdf_params%cthl_2
        thl_1_zt        = pdf_params%thl_1
        thl_2_zt        = pdf_params%thl_2
        varnce_thl_1_zt = pdf_params%varnce_thl_1
        varnce_thl_2_zt = pdf_params%varnce_thl_2
        mixt_frac_zt   = pdf_params%mixt_frac
        rc_1_zt         = pdf_params%rc_1
        rc_2_zt         = pdf_params%rc_2
        rsatl_1_zt        = pdf_params%rsatl_1
        rsatl_2_zt        = pdf_params%rsatl_2
        cloud_frac_1_zt = pdf_params%cloud_frac_1
        cloud_frac_2_zt = pdf_params%cloud_frac_2
        chi_1_zt          = pdf_params%chi_1
        chi_2_zt          = pdf_params%chi_2
        stdev_chi_1_zt    = pdf_params%stdev_chi_1
        stdev_chi_2_zt    = pdf_params%stdev_chi_2
        stdev_eta_1_zt    = pdf_params%stdev_eta_1
        stdev_eta_2_zt    = pdf_params%stdev_eta_2
        corr_w_rt_1_zt    = pdf_params%corr_w_rt_1
        corr_w_rt_2_zt    = pdf_params%corr_w_rt_2
        corr_w_thl_1_zt   = pdf_params%corr_w_thl_1
        corr_w_thl_2_zt   = pdf_params%corr_w_thl_2
        corr_rt_thl_1_zt  = pdf_params%corr_rt_thl_1
        corr_rt_thl_2_zt  = pdf_params%corr_rt_thl_2
        alpha_thl_zt   = pdf_params%alpha_thl
        alpha_rt_zt    = pdf_params%alpha_rt
      end if

      ! If l_call_pdf_closure_twice is true, the _zm variables already have
      ! values from the second call to pdf_closure in advance_clubb_core.
      ! If it is false, the variables are interpolated to the _zm levels.
      if ( l_call_pdf_closure_twice ) then

        ! Store, in locally declared variables, the pdf_params output
        ! from the second call to pdf_closure
        if ( l_apply_rule_to_pdf_params ) then
          w_1_zm          = pdf_params_zm%w_1
          w_2_zm          = pdf_params_zm%w_2
          varnce_w_1_zm   = pdf_params_zm%varnce_w_1
          varnce_w_2_zm   = pdf_params_zm%varnce_w_2
          rt_1_zm         = pdf_params_zm%rt_1
          rt_2_zm         = pdf_params_zm%rt_2
          varnce_rt_1_zm  = pdf_params_zm%varnce_rt_1
          varnce_rt_2_zm  = pdf_params_zm%varnce_rt_2
          crt_1_zm        = pdf_params_zm%crt_1
          crt_2_zm        = pdf_params_zm%crt_2
          cthl_1_zm       = pdf_params_zm%cthl_1
          cthl_2_zm       = pdf_params_zm%cthl_2
          thl_1_zm        = pdf_params_zm%thl_1
          thl_2_zm        = pdf_params_zm%thl_2
          varnce_thl_1_zm = pdf_params_zm%varnce_thl_1
          varnce_thl_2_zm = pdf_params_zm%varnce_thl_2
          mixt_frac_zm   = pdf_params_zm%mixt_frac
          rc_1_zm         = pdf_params_zm%rc_1
          rc_2_zm         = pdf_params_zm%rc_2
          rsatl_1_zm        = pdf_params_zm%rsatl_1
          rsatl_2_zm        = pdf_params_zm%rsatl_2
          cloud_frac_1_zm = pdf_params_zm%cloud_frac_1
          cloud_frac_2_zm = pdf_params_zm%cloud_frac_2
          chi_1_zm          = pdf_params_zm%chi_1
          chi_2_zm          = pdf_params_zm%chi_2
          stdev_chi_1_zm    = pdf_params_zm%stdev_chi_1
          stdev_chi_2_zm    = pdf_params_zm%stdev_chi_2
          stdev_eta_1_zm    = pdf_params_zm%stdev_eta_1
          stdev_eta_2_zm    = pdf_params_zm%stdev_eta_2
          corr_w_rt_1_zm    = pdf_params_zm%corr_w_rt_1
          corr_w_rt_2_zm    = pdf_params_zm%corr_w_rt_2
          corr_w_thl_1_zm   = pdf_params_zm%corr_w_thl_1
          corr_w_thl_2_zm   = pdf_params_zm%corr_w_thl_2
          corr_rt_thl_1_zm  = pdf_params_zm%corr_rt_thl_1
          corr_rt_thl_2_zm  = pdf_params_zm%corr_rt_thl_2
          alpha_thl_zm      = pdf_params_zm%alpha_thl
          alpha_rt_zm       = pdf_params_zm%alpha_rt
        end if

      else

        ! Interpolate thermodynamic variables to the momentum grid.
        ! Since top momentum level is higher than top thermo. level,
        ! set variables at top momentum level to 0.
        wprtp2_zm           = zt2zm( wprtp2 )
        wprtp2_zm(gr%nz) = 0.0_core_rknd
        wpthlp2_zm           = zt2zm( wpthlp2 )
        wpthlp2_zm(gr%nz) = 0.0_core_rknd
        wprtpthlp_zm           = zt2zm( wprtpthlp )
        wprtpthlp_zm(gr%nz)  = 0.0_core_rknd
        cloud_frac_zm          = zt2zm( cloud_frac )
        cloud_frac_zm(gr%nz) = 0.0_core_rknd
        ice_supersat_frac_zm   = zt2zm( ice_supersat_frac )
        ice_supersat_frac_zm(gr%nz) = 0.0_core_rknd
        rcm_zm                 = zt2zm( rcm )
        rcm_zm(gr%nz)        = 0.0_core_rknd
        wp2thvp_zm             = zt2zm( wp2thvp )
        wp2thvp_zm(gr%nz)    = 0.0_core_rknd

        do i = 1, sclr_dim
          wpsclrprtp_zm(:,i)        = zt2zm( wpsclrprtp(:,i) )
          wpsclrprtp_zm(gr%nz,i)  = 0.0_core_rknd
          wpsclrp2_zm(:,i)          = zt2zm( wpsclrp2(:,i) )
          wpsclrp2_zm(gr%nz,i)    = 0.0_core_rknd
          wpsclrpthlp_zm(:,i)       = zt2zm( wpsclrpthlp(:,i) )
          wpsclrpthlp_zm(gr%nz,i) = 0.0_core_rknd
        end do ! i = 1, sclr_dim

        if ( l_apply_rule_to_pdf_params ) then
          w_1_zm                 = zt2zm( pdf_params%w_1 )
          w_1_zm(gr%nz)          = 0.0_core_rknd
          w_2_zm                 = zt2zm( pdf_params%w_2 )
          w_2_zm(gr%nz)          = 0.0_core_rknd
          varnce_w_1_zm          = zt2zm( pdf_params%varnce_w_1 )
          varnce_w_1_zm(gr%nz)   = 0.0_core_rknd
          varnce_w_2_zm          = zt2zm( pdf_params%varnce_w_2 )
          varnce_w_2_zm(gr%nz)   = 0.0_core_rknd
          rt_1_zm                = zt2zm( pdf_params%rt_1 )
          rt_1_zm(gr%nz)         = 0.0_core_rknd
          rt_2_zm                = zt2zm( pdf_params%rt_2 )
          rt_2_zm(gr%nz)         = 0.0_core_rknd
          varnce_rt_1_zm         = zt2zm( pdf_params%varnce_rt_1 )
          varnce_rt_1_zm(gr%nz)  = 0.0_core_rknd
          varnce_rt_2_zm         = zt2zm( pdf_params%varnce_rt_2 )
          varnce_rt_2_zm(gr%nz)  = 0.0_core_rknd
          crt_1_zm               = zt2zm( pdf_params%crt_1 )
          crt_1_zm(gr%nz)        = 0.0_core_rknd
          crt_2_zm               = zt2zm( pdf_params%crt_2 )
          crt_2_zm(gr%nz)        = 0.0_core_rknd
          cthl_1_zm              = zt2zm( pdf_params%cthl_1 )
          cthl_1_zm(gr%nz)       = 0.0_core_rknd
          cthl_2_zm              = zt2zm( pdf_params%cthl_2 )
          cthl_2_zm(gr%nz)       = 0.0_core_rknd
          thl_1_zm               = zt2zm( pdf_params%thl_1 )
          thl_1_zm(gr%nz)        = 0.0_core_rknd
          thl_2_zm               = zt2zm( pdf_params%thl_2 )
          thl_2_zm(gr%nz)        = 0.0_core_rknd
          varnce_thl_1_zm        = zt2zm( pdf_params%varnce_thl_1 )
          varnce_thl_1_zm(gr%nz) = 0.0_core_rknd
          varnce_thl_2_zm        = zt2zm( pdf_params%varnce_thl_2 )
          varnce_thl_2_zm(gr%nz) = 0.0_core_rknd
          mixt_frac_zm          = zt2zm( pdf_params%mixt_frac )
          mixt_frac_zm(gr%nz)   = 0.0_core_rknd
          rc_1_zm                = zt2zm( pdf_params%rc_1 )
          rc_1_zm(gr%nz)         = 0.0_core_rknd
          rc_2_zm                = zt2zm( pdf_params%rc_2 )
          rc_2_zm(gr%nz)         = 0.0_core_rknd
          rsatl_1_zm               = zt2zm( pdf_params%rsatl_1 )
          rsatl_1_zm(gr%nz)        = 0.0_core_rknd
          rsatl_2_zm               = zt2zm( pdf_params%rsatl_2 )
          rsatl_2_zm(gr%nz)        = 0.0_core_rknd
          cloud_frac_1_zm        = zt2zm( pdf_params%cloud_frac_1 )
          cloud_frac_1_zm(gr%nz) = 0.0_core_rknd
          cloud_frac_2_zm        = zt2zm( pdf_params%cloud_frac_2 )
          cloud_frac_2_zm(gr%nz) = 0.0_core_rknd
          chi_1_zm                 = zt2zm( pdf_params%chi_1 )
          chi_1_zm(gr%nz)          = 0.0_core_rknd
          chi_2_zm                 = zt2zm( pdf_params%chi_2 )
          chi_2_zm(gr%nz)          = 0.0_core_rknd
          stdev_chi_1_zm           = zt2zm( pdf_params%stdev_chi_1 )
          stdev_chi_1_zm(gr%nz)    = 0.0_core_rknd
          stdev_chi_2_zm           = zt2zm( pdf_params%stdev_chi_2 )
          stdev_chi_2_zm(gr%nz)    = 0.0_core_rknd
          stdev_eta_1_zm           = zt2zm( pdf_params%stdev_eta_1 )
          stdev_eta_1_zm(gr%nz)    = 0.0_core_rknd
          stdev_eta_2_zm           = zt2zm( pdf_params%stdev_eta_2 )
          stdev_eta_2_zm(gr%nz)    = 0.0_core_rknd
          corr_w_rt_1_zm           = zt2zm( pdf_params%corr_w_rt_1 )
          corr_w_rt_1_zm(gr%nz)    = 0.0_core_rknd
          corr_w_rt_2_zm           = zt2zm( pdf_params%corr_w_rt_2 )
          corr_w_rt_2_zm(gr%nz)    = 0.0_core_rknd
          corr_w_thl_1_zm          = zt2zm( pdf_params%corr_w_thl_1 )
          corr_w_thl_1_zm(gr%nz)   = 0.0_core_rknd
          corr_w_thl_2_zm          = zt2zm( pdf_params%corr_w_thl_2 )
          corr_w_thl_2_zm(gr%nz)   = 0.0_core_rknd
          corr_rt_thl_1_zm         = zt2zm( pdf_params%corr_rt_thl_1 )
          corr_rt_thl_1_zm(gr%nz)  = 0.0_core_rknd
          corr_rt_thl_2_zm         = zt2zm( pdf_params%corr_rt_thl_2 )
          corr_rt_thl_2_zm(gr%nz)  = 0.0_core_rknd
          alpha_thl_zm             = zt2zm( pdf_params%alpha_thl )
          alpha_thl_zm(gr%nz)      = 0.0_core_rknd
          alpha_rt_zm              = zt2zm( pdf_params%alpha_rt )
          alpha_rt_zm(gr%nz)       = 0.0_core_rknd
        end if
      end if ! l_call_pdf_closure_twice

      if ( l_stats ) then
        ! Use the trapezoidal rule to recompute the variables on the stats_zt level
        if ( iwprtp2 > 0 ) then
          wprtp2     = trapezoid_zt( wprtp2, wprtp2_zm )
        end if
        if ( iwpthlp2 > 0 ) then
          wpthlp2    = trapezoid_zt( wpthlp2, wpthlp2_zm )
        end if
        if ( iwprtpthlp > 0 ) then
          wprtpthlp  = trapezoid_zt( wprtpthlp, wprtpthlp_zm )
        end if

        do i = 1, sclr_dim
          if ( iwpsclrprtp(i) > 0 ) then
            wpsclrprtp(:,i)  = trapezoid_zt( wpsclrprtp(:,i), wpsclrprtp_zm(:,i) )
          end if
          if ( iwpsclrpthlp(i) > 0 ) then
            wpsclrpthlp(:,i) = trapezoid_zt( wpsclrpthlp(:,i), wpsclrpthlp_zm(:,i) )
          end if
          if ( iwpsclrp2(i) > 0 ) then
            wpsclrp2(:,i)    = trapezoid_zt( wpsclrp2(:,i), wpsclrp2_zm(:,i) )
          end if
        end do ! i = 1, sclr_dim
      end if ! l_stats

      cloud_frac = trapezoid_zt( cloud_frac, cloud_frac_zm )
      ice_supersat_frac = trapezoid_zt( ice_supersat_frac, ice_supersat_frac_zm )
      rcm        = trapezoid_zt( rcm, rcm_zm )

      wp2thvp    = trapezoid_zt( wp2thvp, wp2thvp_zm )

      if ( l_apply_rule_to_pdf_params ) then
        ! Note: this code makes PDF component cloud water mixing ratios and
        !       cloud fractions inconsistent with the PDF.  Other parts of
        !       CLUBB require PDF component cloud fractions to remain
        !       consistent with the PDF.  This code needs to be refactored
        !       so that cloud_frac_1 and cloud_frac_2 are preserved.
        write(fstderr,*) "The code in l_apply_rule_to_pdf_params does not " &
                         // "preserve cloud_frac_1 and cloud_frac_2 in a " &
                         // "manner consistent with the PDF as required " &
                         // "by other parts of CLUBB."
        write(fstderr,*) "Please refactor before continuing."
        return
        pdf_params%w_1          = trapezoid_zt( w_1_zt, w_1_zm )
        pdf_params%w_2          = trapezoid_zt( w_2_zt, w_2_zm )
        pdf_params%varnce_w_1   = trapezoid_zt( varnce_w_1_zt, varnce_w_1_zm )
        pdf_params%varnce_w_2   = trapezoid_zt( varnce_w_2_zt, varnce_w_2_zm )
        pdf_params%rt_1         = trapezoid_zt( rt_1_zt, rt_1_zm )
        pdf_params%rt_2         = trapezoid_zt( rt_2_zt, rt_2_zm )
        pdf_params%varnce_rt_1  = trapezoid_zt( varnce_rt_1_zt, varnce_rt_1_zm )
        pdf_params%varnce_rt_2  = trapezoid_zt( varnce_rt_2_zt, varnce_rt_2_zm )
        pdf_params%crt_1        = trapezoid_zt( crt_1_zt, crt_1_zm )
        pdf_params%crt_2        = trapezoid_zt( crt_2_zt, crt_2_zm )
        pdf_params%cthl_1       = trapezoid_zt( cthl_1_zt, cthl_1_zm )
        pdf_params%cthl_2       = trapezoid_zt( cthl_2_zt, cthl_2_zm )
        pdf_params%thl_1        = trapezoid_zt( thl_1_zt, thl_1_zm )
        pdf_params%thl_2        = trapezoid_zt( thl_2_zt, thl_2_zm )
        pdf_params%varnce_thl_1 = trapezoid_zt( varnce_thl_1_zt, varnce_thl_1_zm )
        pdf_params%varnce_thl_2 = trapezoid_zt( varnce_thl_2_zt, varnce_thl_2_zm )
        pdf_params%mixt_frac   = trapezoid_zt( mixt_frac_zt, mixt_frac_zm )
        pdf_params%rc_1         = trapezoid_zt( rc_1_zt, rc_1_zm )
        pdf_params%rc_2         = trapezoid_zt( rc_2_zt, rc_2_zm )
        pdf_params%rsatl_1        = trapezoid_zt( rsatl_1_zt, rsatl_1_zm )
        pdf_params%rsatl_2        = trapezoid_zt( rsatl_2_zt, rsatl_2_zm )
        pdf_params%cloud_frac_1 = trapezoid_zt( cloud_frac_1_zt, cloud_frac_1_zm )
        pdf_params%cloud_frac_2 = trapezoid_zt( cloud_frac_2_zt, cloud_frac_2_zm )
        pdf_params%chi_1          = trapezoid_zt( chi_1_zt, chi_1_zm )
        pdf_params%chi_2          = trapezoid_zt( chi_2_zt, chi_2_zm )
        pdf_params%corr_w_rt_1    = trapezoid_zt( corr_w_rt_1_zt, corr_w_rt_1_zm )
        pdf_params%corr_w_rt_2    = trapezoid_zt( corr_w_rt_2_zt, corr_w_rt_2_zm )
        pdf_params%corr_w_thl_1   = trapezoid_zt( corr_w_thl_1_zt, corr_w_thl_1_zm )
        pdf_params%corr_w_thl_2   = trapezoid_zt( corr_w_thl_2_zt, corr_w_thl_2_zm )
        pdf_params%corr_rt_thl_1  = trapezoid_zt( corr_rt_thl_1_zt, corr_rt_thl_1_zm )
        pdf_params%corr_rt_thl_2  = trapezoid_zt( corr_rt_thl_2_zt, corr_rt_thl_2_zm )
        pdf_params%alpha_thl      = trapezoid_zt( alpha_thl_zt, alpha_thl_zm )
        pdf_params%alpha_rt       = trapezoid_zt( alpha_rt_zt, alpha_rt_zm )
        pdf_params%stdev_chi_1    = trapezoid_zt( stdev_chi_1_zt, stdev_chi_1_zm )
        pdf_params%stdev_chi_2    = trapezoid_zt( stdev_chi_2_zt, stdev_chi_2_zm )
        pdf_params%stdev_eta_1    = trapezoid_zt( stdev_eta_1_zt, stdev_eta_1_zm )
        pdf_params%stdev_eta_2    = trapezoid_zt( stdev_eta_2_zt, stdev_eta_2_zm )
      end if

      ! End of trapezoidal rule

      return
    end subroutine trapezoidal_rule_zt

    !-----------------------------------------------------------------------
    subroutine trapezoidal_rule_zm &
               ( wpthvp_zt, thlpthvp_zt, rtpthvp_zt, & ! intent(in)
                 wpthvp, thlpthvp, rtpthvp )           ! intent(inout)
      !
      ! Description:
      !   This subroutine recomputes three variables on the
      !   momentum grid from pdf_closure -- wpthvp, thlpthvp, and
      !   rtpthvp -- by calling the function trapezoid_zm.  Only these three
      !   variables are used in this subroutine because they are the only
      !   pdf_closure momentum variables used elsewhere in CLUBB.
      !
      !   The _zt variables are output from the first call to pdf_closure.
      !   The _zm variables are output from the second call to pdf_closure
      !   on the momentum levels.
      !   This is done before the call to this subroutine.
      !
      !   ldgrant Feb. 2010
      !
      !  References:
      !    None
      !-----------------------------------------------------------------------

      use grid_class, only: gr ! Variable

      use clubb_precision, only: &
        core_rknd ! variable(s)

      implicit none

      ! Input variables
      real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
        wpthvp_zt,   & ! Buoyancy flux (on thermo. grid)  [(K m)/s]
        thlpthvp_zt, & ! th_l' th_v' (on thermo. grid)    [K^2]
        rtpthvp_zt     ! r_t' th_v' (on thermo. grid)     [(kg K)/kg]

      ! Input/Output variables
      real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
        wpthvp,   & ! Buoyancy flux   [(K m)/s]
        thlpthvp, & ! th_l' th_v'     [K^2]
        rtpthvp     ! r_t' th_v'      [(kg K)/kg]

      !----------------------- Begin Code -----------------------------

      ! Use the trapezoidal rule to recompute the variables on the zm level
      wpthvp     = trapezoid_zm( wpthvp, wpthvp_zt )
      thlpthvp   = trapezoid_zm( thlpthvp, thlpthvp_zt )
      rtpthvp    = trapezoid_zm( rtpthvp, rtpthvp_zt )

      return
    end subroutine trapezoidal_rule_zm

    !-----------------------------------------------------------------------
    pure function trapezoid_zt( variable_zt, variable_zm )
      !
      ! Description:
      !   Function which uses the trapezoidal rule from calculus
      !   to recompute the values for the variables on the thermo. grid which
      !   are output from the first call to pdf_closure in module clubb_core.
      !
      !   ldgrant June 2009
      !--------------------------------------------------------------------

      use grid_class, only: gr ! Variable

      use clubb_precision, only: &
        core_rknd ! Variable(s)

      implicit none

      ! Input Variables
      real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
        variable_zt, & ! Variable on the zt grid
        variable_zm    ! Variable on the zm grid

      ! Result
      real( kind = core_rknd ), dimension(gr%nz) :: trapezoid_zt

      ! Local Variable
      integer :: k ! Loop index

      !------------ Begin Code --------------

      ! Boundary condition: trapezoidal rule not valid at zt level 1
      trapezoid_zt(1) = variable_zt(1)

      do k = 2, gr%nz
        ! Trapezoidal rule from calculus
        trapezoid_zt(k) =  0.5_core_rknd * ( variable_zm(k) + variable_zt(k) ) &
                               * ( gr%zm(k) - gr%zt(k) ) * gr%invrs_dzt(k) &
                         + 0.5_core_rknd * ( variable_zt(k) + variable_zm(k-1) ) &
                               * ( gr%zt(k) - gr%zm(k-1) ) * gr%invrs_dzt(k)
      end do ! k = 2, gr%nz

      return
    end function trapezoid_zt

    !-----------------------------------------------------------------------
    pure function trapezoid_zm( variable_zm, variable_zt )
      !
      ! Description:
      !   Function which uses the trapezoidal rule from calculus
      !   to recompute the values for the important variables on the momentum
      !   grid which are output from pdf_closure in module clubb_core.
      !   These momentum variables only include wpthvp, thlpthvp, and rtpthvp.
      !
      !   ldgrant Feb. 2010
      !--------------------------------------------------------------------

      use grid_class, only: gr ! Variable

      use clubb_precision, only: &
        core_rknd ! Variable(s)

      implicit none

      ! Input Variables
      real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
        variable_zm, & ! Variable on the zm grid
        variable_zt    ! Variable on the zt grid

      ! Result
      real( kind = core_rknd ), dimension(gr%nz) :: trapezoid_zm

      ! Local Variable
      integer :: k ! Loop index

      !------------ Begin Code --------------

      ! Boundary conditions: trapezoidal rule not valid at top zm level, nzmax.
      ! Trapezoidal rule also not used at zm level 1.
      trapezoid_zm(1)       = variable_zm(1)
      trapezoid_zm(gr%nz) = variable_zm(gr%nz)

      do k = 2, gr%nz-1
        ! Trapezoidal rule from calculus
        trapezoid_zm(k) =  0.5_core_rknd * ( variable_zt(k+1) + variable_zm(k) ) &
                               * ( gr%zt(k+1) - gr%zm(k) ) * gr%invrs_dzm(k) &
                         + 0.5_core_rknd * ( variable_zm(k) + variable_zt(k) ) &
                               * ( gr%zm(k) - gr%zt(k) ) * gr%invrs_dzm(k)
      end do ! k = 2, gr%nz-1

      return
    end function trapezoid_zm

    !-----------------------------------------------------------------------
    subroutine compute_cloud_cover &
             ( pdf_params, cloud_frac, rcm, & ! intent(in)
               cloud_cover, rcm_in_layer )    ! intent(out)
      !
      ! Description:
      !   Subroutine to compute cloud cover (the amount of sky
      !   covered by cloud) and rcm in layer (liquid water mixing ratio in
      !   the portion of the grid box filled by cloud).
      !
      ! References:
      !   Definition of 's' comes from:
      !   ``The Gaussian Cloud Model Relations'' G. L. Mellor (1977)
      !   JAS, Vol. 34, pp. 356--358.
      !
      ! Notes:
      !   Added July 2009
      !---------------------------------------------------------------------

      use constants_clubb, only: &
          rc_tol, & ! Variable(s)
          fstderr

      use grid_class, only: gr ! Variable

      use pdf_parameter_module, only: &
          pdf_parameter ! Derived data type

      use clubb_precision, only: &
          core_rknd ! Variable(s)

      use error_code, only: &
          clubb_at_least_debug_level  ! Procedure

      implicit none

      ! External functions
      intrinsic :: abs, min, max

      ! Input variables
      real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
        cloud_frac, & ! Cloud fraction             [-]
        rcm           ! Liquid water mixing ratio  [kg/kg]

      type (pdf_parameter), intent(in) :: &
        pdf_params    ! PDF Parameters  [units vary]

      ! Output variables
      real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
        cloud_cover,  & ! Cloud cover                               [-]
        rcm_in_layer    ! Liquid water mixing ratio in cloud layer  [kg/kg]

      ! Local variables
      real( kind = core_rknd ), dimension(gr%nz) :: &
        chi_mean,                & ! Mean extended cloud water mixing ratio of the
    !                            two Gaussian distributions
        vert_cloud_frac_upper, & ! Fraction of cloud in top half of grid box
        vert_cloud_frac_lower, & ! Fraction of cloud in bottom half of grid box
        vert_cloud_frac          ! Fraction of cloud filling the grid box in the vertical

      integer :: k

      ! ------------ Begin code ---------------

      do k = 1, gr%nz

        chi_mean(k) =      pdf_params%mixt_frac(k)  * pdf_params%chi_1(k) + &
                    (1.0_core_rknd-pdf_params%mixt_frac(k)) * pdf_params%chi_2(k)

      end do

      do k = 2, gr%nz-1, 1

        if ( rcm(k) < rc_tol ) then ! No cloud at this level

          cloud_cover(k)  = cloud_frac(k)
          rcm_in_layer(k) = rcm(k)

        else if ( ( rcm(k+1) >= rc_tol ) .and. ( rcm(k-1) >= rc_tol ) ) then
          ! There is cloud above and below,
          !   so assume cloud fills grid box from top to bottom

          cloud_cover(k) = cloud_frac(k)
          rcm_in_layer(k) = rcm(k)

        else if ( ( rcm(k+1) < rc_tol ) .or. ( rcm(k-1) < rc_tol) ) then
          ! Cloud may fail to reach gridbox top or base or both

          ! First let the cloud fill the entire grid box, then overwrite
          ! vert_cloud_frac_upper(k) and/or vert_cloud_frac_lower(k)
          ! for a cloud top, cloud base, or one-point cloud.
          vert_cloud_frac_upper(k) = 0.5_core_rknd
          vert_cloud_frac_lower(k) = 0.5_core_rknd

          if ( rcm(k+1) < rc_tol ) then ! Cloud top

            vert_cloud_frac_upper(k) = &
                     ( ( 0.5_core_rknd / gr%invrs_dzm(k) ) / ( gr%zm(k) - gr%zt(k) ) ) &
                     * ( rcm(k) / ( rcm(k) + abs( chi_mean(k+1) ) ) )

            vert_cloud_frac_upper(k) = min( 0.5_core_rknd, vert_cloud_frac_upper(k) )

            ! Make the transition in cloudiness more gradual than using
            ! the above min statement alone.
            vert_cloud_frac_upper(k) = vert_cloud_frac_upper(k) + &
              ( ( rcm(k+1)/rc_tol )*( 0.5_core_rknd -vert_cloud_frac_upper(k) ) )

          else

            vert_cloud_frac_upper(k) = 0.5_core_rknd

          end if

          if ( rcm(k-1) < rc_tol ) then ! Cloud base

            vert_cloud_frac_lower(k) = &
                     ( ( 0.5_core_rknd / gr%invrs_dzm(k-1) ) / ( gr%zt(k) - gr%zm(k-1) ) ) &
                     * ( rcm(k) / ( rcm(k) + abs( chi_mean(k-1) ) ) )

            vert_cloud_frac_lower(k) = min( 0.5_core_rknd, vert_cloud_frac_lower(k) )

            ! Make the transition in cloudiness more gradual than using
            ! the above min statement alone.
            vert_cloud_frac_lower(k) = vert_cloud_frac_lower(k) + &
              ( ( rcm(k-1)/rc_tol )*( 0.5_core_rknd -vert_cloud_frac_lower(k) ) )

          else

            vert_cloud_frac_lower(k) = 0.5_core_rknd

          end if

          vert_cloud_frac(k) = &
            vert_cloud_frac_upper(k) + vert_cloud_frac_lower(k)

          vert_cloud_frac(k) = &
            max( cloud_frac(k), min( 1.0_core_rknd, vert_cloud_frac(k) ) )

          cloud_cover(k)  = cloud_frac(k) / vert_cloud_frac(k)
          rcm_in_layer(k) = rcm(k) / vert_cloud_frac(k)

        else

          if ( clubb_at_least_debug_level( 0 ) ) then

            write(fstderr,*)  &
               "Error: Should not arrive here in computation of cloud_cover"

            write(fstderr,*) "At grid level k = ", k
            write(fstderr,*) "pdf_params(k)%mixt_frac = ", pdf_params%mixt_frac(k)
            write(fstderr,*) "pdf_params(k)%chi_1 = ", pdf_params%chi_1(k)
            write(fstderr,*) "pdf_params(k)%chi_2 = ", pdf_params%chi_2(k)
            write(fstderr,*) "cloud_frac(k) = ", cloud_frac(k)
            write(fstderr,*) "rcm(k) = ", rcm(k)
            write(fstderr,*) "rcm(k+1) = ", rcm(k+1)
            write(fstderr,*) "rcm(k-1) = ", rcm(k-1)

            return

          end if

        end if ! rcm(k) < rc_tol

      end do ! k = 2, gr%nz-1, 1

      cloud_cover(1)       = cloud_frac(1)
      cloud_cover(gr%nz) = cloud_frac(gr%nz)

      rcm_in_layer(1)       = rcm(1)
      rcm_in_layer(gr%nz) = rcm(gr%nz)

      return
    end subroutine compute_cloud_cover
    !-----------------------------------------------------------------------
    subroutine clip_rcm &
             ( rtm, message, & ! intent(in)
               rcm )    ! intent(inout)
      !
      ! Description:
      !   Subroutine that reduces cloud water (rcm) whenever
      !   it exceeds total water (rtm = vapor + liquid).
      !   This avoids negative values of rvm = water vapor mixing ratio.
      !   However, it will not ensure that rcm <= rtm if rtm <= 0.
      !
      ! References:
      !   None
      !---------------------------------------------------------------------


      use grid_class, only: gr ! Variable

      use error_code, only: &
        clubb_at_least_debug_level  ! Procedure

      use constants_clubb, only: &
        fstderr, & ! Variable(s)
        zero_threshold

      use clubb_precision, only: &
        core_rknd ! Variable(s)

      implicit none

      ! External functions
      intrinsic :: max, epsilon

      ! Input variables
      real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
        rtm           ! Total water mixing ratio             [kg/kg]

      character(len= * ), intent(in) :: message

      real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
        rcm           ! Cloud water mixing ratio  [kg/kg]

      integer :: k

      ! ------------ Begin code ---------------

      ! Vince Larson clipped rcm in order to prevent rvm < 0.  5 Apr 2008.
      ! This code won't work unless rtm >= 0 !!!
      ! We do not clip rcm_in_layer because rcm_in_layer only influences
      ! radiation, and we do not want to bother recomputing it.  6 Aug 2009
      do k = 1, gr%nz
        if ( rtm(k) < rcm(k) ) then

          if ( clubb_at_least_debug_level( 1 ) ) then
            write(fstderr,*) message, ' at k=', k, 'rcm(k) = ', rcm(k), &
              'rtm(k) = ', rtm(k), '.',  '  Clipping rcm.'

          end if ! clubb_at_least_debug_level( 1 )

          rcm(k) = max( zero_threshold, rtm(k) - epsilon( rtm(k) ) )

        end if ! rtm(k) < rcm(k)

      end do ! k=1..gr%nz

      return
    end subroutine clip_rcm

    !-----------------------------------------------------------------------------
    subroutine set_Lscale_max( l_implemented, host_dx, host_dy, &
                               Lscale_max )

      ! Description:
      !   This subroutine sets the value of Lscale_max, which is the maximum
      !   allowable value of Lscale.  For standard CLUBB, it is set to a very large
      !   value so that Lscale will not be limited.  However, when CLUBB is running
      !   as part of a host model, the value of Lscale_max is dependent on the size
      !   of the host model's horizontal grid spacing.  The smaller the host model's
      !   horizontal grid spacing, the smaller the value of Lscale_max.  When Lscale
      !   is limited to a small value, the value of time-scale Tau is reduced, which
      !   in turn produces greater damping on CLUBB's turbulent parameters.  This
      !   is the desired effect on turbulent parameters for a host model with small
      !   horizontal grid spacing, for small areas usually contain much less
      !   variation in meteorological quantities than large areas.

      ! References:
      !   None
      !-----------------------------------------------------------------------

      use clubb_precision, only: &
        core_rknd ! Variable(s)

      implicit none

      ! Input Variables
      logical, intent(in) :: &
        l_implemented    ! Flag to see if CLUBB is running on it's own,
      !                    or if it's implemented as part of a host model.

      real( kind = core_rknd ), intent(in) :: &
        host_dx, & ! Host model's east-west horizontal grid spacing     [m]
        host_dy    ! Host model's north-south horizontal grid spacing   [m]

      ! Output Variable
      real( kind = core_rknd ), intent(out) :: &
        Lscale_max    ! Maximum allowable value for Lscale   [m]

      ! ---- Begin Code ----

      ! Determine the maximum allowable value for Lscale (in meters).
      if ( l_implemented ) then
        Lscale_max = 0.25_core_rknd * min( host_dx, host_dy )
      else
        Lscale_max = 1.0e5_core_rknd
      end if

      return
    end subroutine set_Lscale_max

!===============================================================================
  pure subroutine calculate_thlp2_rad &
                  ( nz, rcm_zm, thlprcp, radht_zm, &      ! Intent(in)
                    thlp2_forcing )                       ! Intent(inout)

  ! Description:
  !   Computes the contribution of radiative cooling to thlp2

  ! References:
  !   See clubb:ticket:632
  !----------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd                     ! Constant(s)

  use grid_class, only:  &
    zt2zm                         ! Procedure

  use constants_clubb, only: &
    two, &
    rc_tol

  use parameters_tunable, only: &
    thlp2_rad_coef                 ! Variable(s)

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

  ! Local Variables
  integer :: &
    k                     ! Loop iterator                                  [-]

  !----------------------------------------------------------------------


    do k = 1, nz

       if ( rcm_zm(k) > rc_tol ) then

          thlp2_forcing(k) = thlp2_forcing(k) + &
                    thlp2_rad_coef * ( two ) * radht_zm(k) / rcm_zm(k) * thlprcp(k)

       end if

    end do


    return
  end subroutine calculate_thlp2_rad


    !-----------------------------------------------------------------------

end module advance_clubb_core_module
