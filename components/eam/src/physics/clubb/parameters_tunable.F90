!-----------------------------------------------------------------------
! $Id$ 
!===============================================================================
module parameters_tunable
 
  ! Description:
  !   This module contains tunable model parameters.  The purpose of the module is to make it
  !   easier for the clubb_tuner code to use the params vector without "knowing" any information
  !   about the individual parameters contained in the vector itself.  It makes it easier to add
  !   new parameters to be tuned for, but does not make the CLUBB_core code itself any simpler.
  !   The parameters within the vector do not need to be the same variables used in the rest of
  !   CLUBB_core (see for e.g. nu1_vert_res_dep or lmin_coef).
  !   The parameters in the params vector only need to be those parameters for which we're not
  !   sure the correct value and we'd like to tune for.
  !
  ! References:
  !   None
  ! 
  ! Notes:
  !   To make it easier to verify of code correctness, please keep the omp threadprivate
  !   directives just after the variable declaration.  All parameters in this
  !   module should be declared threadprivate because of the CLUBB tuner.
  !-----------------------------------------------------------------------

  use constants_clubb, only: eps ! Epsilon

  use parameter_indices, only: nparams ! Variable(s)

  use clubb_precision, only: &
      core_rknd ! Variable(s)

  implicit none

  ! Default to private
  private

  public :: set_default_parameters, setup_parameters, read_parameters, &
            read_param_minmax, read_param_constraints, &
            adj_low_res_nu, nu_vertical_res_dep

  type nu_vertical_res_dep
    real( kind = core_rknd ) :: & 
      nu1,   & ! Background Coefficient of Eddy Diffusion: wp2      [m^2/s]
      nu2,   & ! Background Coefficient of Eddy Diffusion: xp2      [m^2/s]
      nu6,   & ! Background Coefficient of Eddy Diffusion: wpxp     [m^2/s]
      nu8,   & ! Background Coefficient of Eddy Diffusion: wp3      [m^2/s]
      nu9,   & ! Background Coefficient of Eddy Diffusion: up2/vp2  [m^2/s]
      nu10,  & ! Background Coefficient of Eddy Diffusion: edsclrm  [m^2/s]
      nu_hm    ! Background Coefficient of Eddy Diffusion: hydromet [m^2/s]
  end type nu_vertical_res_dep

  ! These are referenced together often enough that it made sense to
  ! make a list of them.  Note that lmin_coef is the input parameter,
  ! while the actual lmin model constant is computed from this.
  !***************************************************************
  !                    ***** IMPORTANT *****
  ! If you change the order of the parameters in the parameter_indices,
  ! you will need to change the order of this list as well or the
  ! tuner will break!
  !                    ***** IMPORTANT *****
  !***************************************************************
  character(len=28), dimension(nparams), parameter, public ::  & 
  params_list = & 
     (/"C1                          ", "C1b                         ", &
       "C1c                         ", &
       "C2rt                        ", "C2thl                       ", &
       "C2rtthl                     ", "C4                          ", &
       "C_uu_shr                    ", "C_uu_buoy                   ", &
       "C6rt                        ", &
       "C6rtb                       ", "C6rtc                       ", &
       "C6thl                       ", "C6thlb                      ", &
       "C6thlc                      ", "C7                          ", &
       "C7b                         ", "C7c                         ", &
       "C8                          ", "C8b                         ", &
       "C10                         ", "C11                         ", &
       "C11b                        ", "C11c                        ", &
       "C12                         ", "C13                         ", &
       "C14                         ", "C_wp2_pr_dfsn               ", &
       "C_wp3_pr_tp                 ", "C_wp3_pr_turb               ", &
       "C_wp3_pr_dfsn               ", "C_wp2_splat                 ", &
       "C6rt_Lscale0                ", "C6thl_Lscale0               ", &
       "C7_Lscale0                  ", "wpxp_L_thresh               ", &
       "c_K                         ", "c_K1                        ", &
       "nu1                         ", "c_K2                        ", &
       "nu2                         ", "c_K6                        ", &
       "nu6                         ", "c_K8                        ", &
       "nu8                         ", "c_K9                        ", &
       "nu9                         ", "nu10                        ", &
       "c_K_hm                      ", "c_K_hmb                     ", &
       "K_hm_min_coef               ", "nu_hm                       ", &
       "slope_coef_spread_DG_means_w", "pdf_component_stdev_factor_w", &
       "coef_spread_DG_means_rt     ", "coef_spread_DG_means_thl    ", &
       "gamma_coef                  ", "gamma_coefb                 ", &
       "gamma_coefc                 ", "mu                          ", &
       "beta                        ", "lmin_coef                   ", &
       "omicron                     ", "zeta_vrnce_rat              ", &
       "upsilon_precip_frac_rat     ", "lambda0_stability_coef      ", &
       "mult_coef                   ", "taumin                      ", &
       "taumax                      ", "Lscale_mu_coef              ", &
       "Lscale_pert_coef            ", "alpha_corr                  ", &
       "Skw_denom_coef              ", "c_K10                       ", &
       "c_K10h                      ", "thlp2_rad_coef              ", &
       "thlp2_rad_cloud_frac_thresh ", "up2_sfc_coef                ", &
       "Skw_max_mag                 ", "C_invrs_tau_bkgnd           ", &
       "C_invrs_tau_sfc             ", "C_invrs_tau_shear           ", &
       "C_invrs_tau_N2              ", "C_invrs_tau_N2_wp2          ", &
       "C_invrs_tau_N2_xp2          ", "C_invrs_tau_N2_wpxp         ", &
       "C_invrs_tau_N2_clear_wp3    ", "C_invrs_tau_wpxp_Ri         ", &
       "C_invrs_tau_wpxp_N2_thresh  ", "xp3_coef_base               ", &
       "xp3_coef_slope              ", "altitude_threshold          ", &
       "rtp2_clip_coef              ", "Cx_min                      ", &
       "Cx_max                      ", "Richardson_num_min          ", &
       "Richardson_num_max          "/)

  real( kind = core_rknd ), parameter, private :: &
    init_value = -999._core_rknd ! Initial value for the parameters, used to detect missing values

#ifdef E3SM
  public :: clubb_param_readnl

  ! The parameters below have the same meaning as those without prefix 'clubb_'
  ! They can be specified via namelist to ease parameter tuning
  ! If adding more parameters for tuning via namelist, need to insert blocks
  ! accordingly in subroutines read_parameters and clubb_paramm_readnl
  ! (including updating list of nml variables)
  real( kind = core_rknd ) ::           &
    clubb_C1,                           &
    clubb_C1b,                          &
    clubb_C1c,                          &
    clubb_C2rt,                         &
    clubb_C2thl,                        &
    clubb_C2rtthl,                      &
    clubb_C4,                           &
    clubb_C_uu_shr,                     &
    clubb_C_uu_buoy,                    &
    clubb_C6rt,                         &
    clubb_C6rtb,                        &
    clubb_C6rtc,                        &
    clubb_C6thl,                        &
    clubb_C6thlb,                       &
    clubb_C6thlc,                       &
    clubb_C7,                           &
    clubb_C7b,                          &
    clubb_C8,                           &
    clubb_C8b,                          &
    clubb_C11,                          &
    clubb_C11b,                         &
    clubb_C11c,                         &
    clubb_C14,                          &
    clubb_C_wp2_pr_dfsn,                &
    clubb_C_wp3_pr_tp,                  &
    clubb_C_wp3_pr_turb,                &
    clubb_C_wp3_pr_dfsn,                &
    clubb_beta,                         &
    clubb_gamma_coef,                   &
    clubb_gamma_coefb,                  &
    clubb_gamma_coefc,                  &
    clubb_pdf_component_stdev_factor_w, &
    clubb_mu,                           &
    clubb_c_K1,                         &
    clubb_nu1,                          &
    clubb_c_K2,                         &
    clubb_nu2,                          &
    clubb_c_K8,                         &
    clubb_nu8,                          &
    clubb_c_K9,                         &
    clubb_nu9,                          &
    clubb_c_K10,                        &
    clubb_c_K10h,                       &
    clubb_c_K_hmb,                      &
    clubb_wpxp_L_thresh,                &
    clubb_lmin_coef,                    &
    clubb_mult_coef,                    &
    clubb_Skw_denom_coef,               &
    clubb_up2_sfc_coef,                 &
    clubb_Skw_max_mag,                  &
    clubb_C_invrs_tau_bkgnd,            &
    clubb_C_invrs_tau_sfc,              &
    clubb_C_invrs_tau_shear,            &
    clubb_C_invrs_tau_N2,               &
    clubb_C_invrs_tau_N2_wp2,           &
    clubb_C_invrs_tau_N2_xp2,           &
    clubb_C_invrs_tau_N2_wpxp,          &
    clubb_C_invrs_tau_N2_clear_wp3,     &
    clubb_C_invrs_tau_wpxp_Ri,          &
    clubb_C_invrs_tau_wpxp_N2_thresh,   &
    clubb_C_wp2_splat,                  &
    clubb_xp3_coef_base,                &
    clubb_xp3_coef_slope,               &
    clubb_altitude_threshold,           &
    clubb_rtp2_clip_coef,               &
    clubb_Cx_min,                       &
    clubb_Cx_max,                       &
    clubb_Richardson_num_min,           & 
    clubb_Richardson_num_max 
    
!$omp threadprivate(clubb_C1, clubb_C1b, clubb_C1c, &
!$omp   clubb_C2rt, clubb_C2thl, clubb_C2rtthl, clubb_C4, &
!$omp   clubb_C_uu_shr, clubb_C_uu_buoy, clubb_C6rt, clubb_C6rtb, clubb_C6rtc, &
!$omp   clubb_C6thl, clubb_C6thlb, clubb_C6thlc, &
!$omp   clubb_C7, clubb_C7b, clubb_C8, clubb_C8b, clubb_C11, clubb_C11b, &
!$omp   clubb_C11c, clubb_C14, clubb_C_wp2_pr_dfsn, & 
!$omp   clubb_C_wp3_pr_tp, clubb_C_wp3_pr_turb, clubb_C_wp3_pr_dfsn, &
!$omp   clubb_beta, clubb_gamma_coef, clubb_gamma_coefb, clubb_gamma_coefc, &
!$omp   clubb_pdf_component_stdev_factor_w, clubb_mu, clubb_c_K1, clubb_nu1, &
!$omp   clubb_c_K2, clubb_nu2, clubb_c_K8, clubb_nu8, clubb_c_K9, clubb_nu9, &
!$omp   clubb_c_K10, clubb_c_K10h, clubb_c_K_hmb, clubb_wpxp_L_thresh, &
!$omp   clubb_lmin_coef, clubb_mult_coef, clubb_Skw_denom_coef, &
!$omp   clubb_up2_sfc_coef, clubb_Skw_max_mag, clubb_C_invrs_tau_bkgnd, &
!$omp   clubb_C_invrs_tau_sfc, clubb_C_invrs_tau_shear, clubb_C_invrs_tau_N2, &
!$omp   clubb_C_invrs_tau_N2_wp2, clubb_C_invrs_tau_N2_xp2, &
!$omp   clubb_C_invrs_tau_N2_wpxp, clubb_C_invrs_tau_N2_clear_wp3, &
!$omp   clubb_C_invrs_tau_wpxp_Ri, clubb_C_invrs_tau_wpxp_N2_thresh, &
!$omp   clubb_C_wp2_splat, clubb_xp3_coef_base, clubb_xp3_coef_slope, &
!$omp   clubb_altitude_threshold, clubb_rtp2_clip_coef, clubb_Cx_min, &
!$omp   clubb_Cx_max, clubb_Richardson_num_min, clubb_Richardson_num_max)
    
#endif /*E3SM*/

  contains

  !=============================================================================
  subroutine set_default_parameters( &
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
               Cx_min, Cx_max, Richardson_num_min, Richardson_num_max )

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
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max


    ! NOTE: In CLUBB standalone, as well as some host models, the hardcoded
    !       default values of some or all of the parameters below have no
    !       effect, as the values are simply read in using a namelist or set in
    !       host model specific code.

    ! Model tunable parameters
    ! WARNING: THE DEFAULT VALUES OF THE PARAMETERS BELOW MAY BE OVERWRITTEN
    !    BY NAMELIST VALUES FROM, E.G., clubb_params_nl!!!
    C1          = 1.000000_core_rknd ! Low Skewness in C1 Skw. Function    [-]
    C1b         = 1.000000_core_rknd ! High Skewness in C1 Skw. Function   [-]
    C1c         = 1.000000_core_rknd ! Degree of Slope of C1 Skw. Function [-]
    C2rt        = 2.000000_core_rknd ! C2 coef. for the rtp2_dp1 term      [-]
    C2thl       = 2.000000_core_rknd ! C2 coef. for the thlp2_dp1 term     [-]
    C2rtthl     = 2.000000_core_rknd ! C2 coef. for the rtpthlp_dp1 term   [-]
    C4          = 2.000000_core_rknd ! Used only when l_tke_aniso is true  [-]
    C_uu_shr    = 0.400000_core_rknd ! Coef. in pressure terms (shear): w'^2 eqn   [-]
    C_uu_buoy   = 0.300000_core_rknd ! Coef. in pressure terms (buoyancy): w'^2 eqn [-]
    C6rt        = 2.000000_core_rknd ! Low Skewness in C6rt Skw. Function  [-]
    C6rtb       = 2.000000_core_rknd ! High Skewness in C6rt Skw. Function [-]
    C6rtc       = 1.000000_core_rknd ! Degree of Slope of C6rt Skw. Fnct.  [-]
    C6thl       = 2.000000_core_rknd ! Low Skewness in C6thl Skw. Function [-]
    C6thlb      = 2.000000_core_rknd ! High Skewness in C6thl Skw. Fnct.   [-]
    C6thlc      = 1.000000_core_rknd ! Degree of Slope of C6thl Skw. Fnct. [-]
    C7          = 0.500000_core_rknd ! Low Skewness in C7 Skw. Function    [-]
    C7b         = 0.500000_core_rknd ! High Skewness in C7 Skw. Function   [-]
    C7c         = 0.500000_core_rknd ! Degree of Slope of C7 Skw. Function [-]
    C8          = 0.500000_core_rknd ! Coef. #1 in C8 Skewness Equation    [-]
    C8b         = 0.020000_core_rknd ! Coef. #2 in C8 Skewness Equation    [-]
    C10         = 3.300000_core_rknd ! Currently Not Used in the Model     [-]
    C11         = 0.400000_core_rknd ! Low Skewness in C11 Skw. Function   [-]
    C11b        = 0.400000_core_rknd ! High Skewness in C11 Skw. Function  [-]
    C11c        = 0.500000_core_rknd ! Degree of Slope of C11 Skw. Fnct.   [-]
    C12         = 1.000000_core_rknd ! Constant in w'^3 Crank-Nich. diff.  [-]
    C13         = 0.100000_core_rknd ! Not currently used in model         [-]
    C14         = 1.000000_core_rknd ! Constant for u'^2 and v'^2 terms    [-]
    C_wp2_pr_dfsn = 0.000000_core_rknd ! Coefficient for the wp2_pr_dfsn term [-]
    C_wp3_pr_tp   = 0.000000_core_rknd ! Coefficient for the wp3_pr_tp term [-]
    C_wp3_pr_turb = 0.000000_core_rknd ! Coefficient for the wp3_pr_turb term [-]
    C_wp3_pr_dfsn = 0.000000_core_rknd ! Coefficient for the wp3_pr_dfsn term [-]
    C_wp2_splat   = 2.000000_core_rknd ! Coefficient for gustiness near ground [-]
    C6rt_Lscale0  = 14.0_core_rknd     ! Damp C6rt as a fnct. of Lscale  [-]
    C6thl_Lscale0 = 14.0_core_rknd     ! Damp C6thl as a fnct. of Lscale [-]
    C7_Lscale0    = 0.8500000_core_rknd ! Damp C7 as a fnct. of Lscale    [-]
    wpxp_L_thresh = 60.0_core_rknd      ! Lscale threshold: damp C6 & C7  [m]
    c_K         = 0.200000_core_rknd ! Constant C_mu^(1/4) in DD 1987 [m^2/s]
    c_K1        = 0.200000_core_rknd ! Coef. of Eddy Diffusion: wp2   [m^2/s]
    c_K2        = 0.025000_core_rknd ! Coef. of Eddy Diffusion: xp2   [m^2/s]
    c_K6        = 0.375000_core_rknd ! Coef. of Eddy Diffusion: wpxp  [m^2/s]
    c_K8        = 5.000000_core_rknd ! Coef. of Eddy Diffusion: wp3   [m^2/s]
    c_K9        = 0.100000_core_rknd ! Coef. of Eddy Diff.: up2/vp2   [m^2/s]
    c_K_hm      = 0.750000_core_rknd ! Coef. of Eddy Diffusion: hmm   [m^2/s]
    c_K_hmb     = 0.750000_core_rknd ! Coef. of Non-Local Factor, Eddy Diffusion: hmm   [m^2/s]
    K_hm_min_coef = 0.10000_core_rknd ! Min. of Non-Local Factor, Eddy Diffusion: hmm   [m^2/s]
    gamma_coef  = 0.250000_core_rknd ! Low Skw.: gamma coef. Skw. Fnct.   [-]
    gamma_coefb = 0.250000_core_rknd ! High Skw.: gamma coef. Skw. Fnct.  [-]
    gamma_coefc = 5.000000_core_rknd ! Deg. Slope: gamma coef. Skw. Fnct. [-]
    mu          = 1.000E-3_core_rknd ! Fract entrain rate per unit alt  [1/m]
    mult_coef   = 0.500000_core_rknd ! Coef. applied to log(avg dz/thresh)[-]
    taumin      = 90.00000_core_rknd ! Min. allow. value: time-scale tau  [s]
    taumax      = 3600.000_core_rknd ! Max. allow. value: time-scale tau  [s]
    Lscale_mu_coef   = 2.0_core_rknd ! Coef perturb mu: av calc Lscale    [-]
    Lscale_pert_coef = 0.1_core_rknd ! Coef pert thlm/rtm: av calc Lscale [-]
    alpha_corr = 0.15_core_rknd   ! Coef. for the corr. diagnosis algorithm  [-]
    nu1   = 20.00000_core_rknd ! Bg. Coef. Eddy Diffusion: wp2        [m^2/s]
    nu2   = 1.000000_core_rknd ! Bg. Coef. Eddy Diffusion: xp2        [m^2/s]
    nu6   = 5.000000_core_rknd ! Bg. Coef. Eddy Diffusion: wpxp       [m^2/s]
    nu8   = 20.00000_core_rknd ! Bg. Coef. Eddy Diffusion: wp3        [m^2/s]
    nu9   = 10.00000_core_rknd ! Bg. Coef. Eddy Diffusion: up2/vp2    [m^2/s]
    nu10  = 0.000000_core_rknd ! Bg. Coef. Eddy Diffusion: edsclrm    [m^2/s]
    nu_hm = 1.500000_core_rknd ! Bg. Coef. Eddy Diffusion: hmm        [m^2/s]

    ! Vince Larson added a constant to set plume widths for theta_l and rt
    ! beta should vary between 0 and 3.
    beta = 1.000000_core_rknd    ! Beta coefficient     [-]
    lmin_coef = 0.500000_core_rknd   ! Coefficient of lmin    [-]
    Skw_max_mag = 10.0_core_rknd     ! Max magnitude of skewness [-]

    C_invrs_tau_bkgnd          = 1.1_core_rknd 
    C_invrs_tau_sfc            = 0.1_core_rknd
    C_invrs_tau_shear          = 0.15_core_rknd
    C_invrs_tau_N2             = 0.4_core_rknd 
    C_invrs_tau_N2_wp2         = 0.2_core_rknd
    C_invrs_tau_N2_xp2         = 0.05_core_rknd
    C_invrs_tau_N2_wpxp        = 0.0_core_rknd
    C_invrs_tau_N2_clear_wp3   = 1.0_core_rknd
    C_invrs_tau_wpxp_Ri        = 0.35_core_rknd
    C_invrs_tau_wpxp_N2_thresh = 3.3e-4_core_rknd

    ! Parameters for the new PDF (w, rt, and theta-l).
    !
    ! Brian Griffin added a tunable parameter for the PDF of w,
    ! slope_coef_spread_DG_means_w, to increase or decrease the spread between
    ! the two PDF component means of w.  When the value of this slope parameter
    ! is larger, F_w is smaller and the PDF component means of w are closer
    ! together.
    ! Valid values are slope_coef_spread_DG_means_w > 0.
    !
    ! A second parameter for the PDF of w, pdf_component_stdev_factor_w, is used
    ! to adjust the standard deviations of the 1st PDF component against the 2nd
    ! PDF component for w.  This parameter is related to zeta_w, where:
    !
    ! 1 + zeta_w = ( mixt_frac * sigma_w_1^2 )
    !              / ( ( 1 - mixt_frac ) * sigma_w_2^2 );
    !
    ! The pdf_component_stdev_factor_w is set such that:
    !
    ! pdf_component_stdev_factor_w = zeta_w + 1.
    !
    ! Valid values are pdf_component_stdev_factor_w > 0.
    !
    ! The parameter for the PDF of rt is coef_spread_DG_means_rt.  Valid values
    ! are 0 <= coef_spread_DG_means_rt < 1.  When coef_spread_DG_means_rt
    ! approaches 0, F_rt approaches min_F_rt, and the two PDF component means
    ! become closer together.  When coef_spread_DG_means_rt approaches 1, F_rt
    ! approaches max_F_rt, and the two PDF component means are spread farther
    ! apart.
    !
    ! The parameter for the PDF of theta-l is coef_spread_DG_means_thl.
    ! Valid values are 0 <= coef_spread_DG_means_thl < 1.  When
    ! coef_spread_DG_means_thl approaches 0, F_thl approaches min_F_thl, and the
    ! two PDF component means become closer together.  When
    ! coef_spread_DG_means_thl approaches 1, F_thl approaches max_F_thl, and the
    ! two PDF component means are spread farther apart.
    ! Slope coefficient for the spread between the PDF component means of w.
    slope_coef_spread_DG_means_w = 21.0_core_rknd
    ! Parameter to adjust the PDF component standard deviations of w.
    pdf_component_stdev_factor_w = 1.0_core_rknd
    ! Coefficient for the spread between the PDF component means of rt.
    coef_spread_DG_means_rt = 0.8_core_rknd
    ! Coefficient for the spread between the PDF component means of thl.
    coef_spread_DG_means_thl = 0.8_core_rknd

    ! Parameters for the hydrometeor portion of the PDF.
    !
    ! Brian Griffin added a parameter for hydrometeors, omicron, to increase the
    ! standard deviation of each component and decrease the spread between the
    ! component means as the value of omicron inreases.  Valid value are
    ! 0 < omicron <= 1.
    ! A second parameter for hydrometeors, zeta, increases the standard
    ! deviation of component 1 at the expense of the standard deviation of
    ! component 2 when the value of zeta > 0 (and increasingly so as zeta
    ! increases).  Valid values are zeta > -1.
    omicron        = 0.5_core_rknd ! Hydromet width/spread-of-means param [-]
    zeta_vrnce_rat = 0.0_core_rknd ! Ratio sigma^2/mu^2 comp. 1 / comp. 2 [-]
    ! ratio mixt_frac*precip_frac_1/precip_frac (precip_frac_calc_type=2)    [-]
    upsilon_precip_frac_rat = 0.55_core_rknd

    ! Intensity of stability correction applied to C1 and C6 [-]
    lambda0_stability_coef = 0.03_core_rknd
    ! Factor to decrease sensitivity in the denominator of Skw calculation
    Skw_denom_coef = 4.0_core_rknd

    ! Momentum coefficient of Kh_zm
    c_K10 = 1.0_core_rknd
    ! Thermodynamic coefficient of Kh_zm
    c_K10h = 1.0_core_rknd

    thlp2_rad_coef = 1.0_core_rknd ! Coefficient of thlp2_rad               [-]
    thlp2_rad_cloud_frac_thresh = 0.1_core_rknd ! Minimum cloud fraction for
                                                ! computation of thlp2_rad  [-]

    up2_sfc_coef = 4.0_core_rknd ! Coefficients of up2 and vp2    [-]

    xp3_coef_base  = 0.25_core_rknd ! "Base" value of xp3_coef in simple eqn
    xp3_coef_slope = 0.01_core_rknd ! Slope in regards to Brunt-Vaisla freq.

    altitude_threshold = 100.0_core_rknd ! Altitude above which damping should occur for wpxp
    rtp2_clip_coef = 0.5_core_rknd       ! Coef. appled the clipping threshold on rtp2

    Cx_min = 0.33_core_rknd              ! Threshold on Cx_fnc_Richardson
    Cx_max = 0.95_core_rknd              ! Threshold on Cx_fnc_Richardson
    Richardson_num_min = 0.25_core_rknd  ! Threshold on Richardson number
    Richardson_num_max = 400.0_core_rknd ! Threshold on Richardson number


    return

  end subroutine set_default_parameters

  !=============================================================================
  subroutine setup_parameters & 
            ( deltaz, params, nzmax, &
              grid_type, momentum_heights, thermodynamic_heights, &
              l_prescribed_avg_deltaz, &
              lmin, nu_vert_res_dep, err_code_out )

    ! Description:
    ! Subroutine to setup model parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------


    use constants_clubb, only: &
        three,   & ! Variable(s)
        one,     &
        zero,    &
        fstderr

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    use parameter_indices, only: &
        izeta_vrnce_rat

    implicit none

    ! Constant Parameters
    real( kind = core_rknd ), parameter :: &
      lmin_deltaz = 40.0_core_rknd ! Fixed value for minimum value for the length scale.

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

    real( kind = core_rknd ), intent(out) :: &
      lmin    ! Min. value for the length scale    [m]

    type(nu_vertical_res_dep), intent(out) :: &
      nu_vert_res_dep    ! Vertical resolution dependent nu values

    integer, intent(out) :: &
      err_code_out  ! Error code indicator

    integer :: k    ! loop variable

    real( kind = core_rknd ) :: & 
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
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max

    !-------------------- Begin code --------------------

    ! Ensure all variables are greater than 0, and zeta_vrnce_rat is greater than -1
    do k = 1, nparams

        if ( k /= izeta_vrnce_rat .and. params(k) < zero ) then

            write(fstderr,*) params_list(k), " = ", params(k)
            write(fstderr,*) params_list(k), " must satisfy 0.0 <= ", params_list(k)
            err_code = clubb_fatal_error

        else if ( params(k) < -one ) then

            write(fstderr,*) "zeta_vrnce_rat = ", zeta_vrnce_rat
            write(fstderr,*) "zeta_vrnce_rat must satisfy -1.0 <= zeta_vrnce_rat"
            err_code = clubb_fatal_error

        end if

    end do

    call unpack_parameters & 
             ( params, & ! intent(in)
               C1, C1b, C1c, C2rt, C2thl, C2rtthl, & ! intent(out)
               C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & ! intent(out)
               C7, C7b, C7c, C8, C8b, C10, & ! intent(out)
               C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, & ! intent(out)
               C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, & ! intent(out)
               C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, & ! intent(out)
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, & ! intent(out)
               c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, & ! intent(out)
               slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, & ! intent(out)
               coef_spread_DG_means_rt, coef_spread_DG_means_thl, & ! intent(out)
               gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, & ! intent(out)
               omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, & ! intent(out)
               lambda0_stability_coef, mult_coef, taumin, taumax, & ! intent(out)
               Lscale_mu_coef, Lscale_pert_coef, alpha_corr, & ! intent(out)
               Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, & ! intent(out)
               thlp2_rad_cloud_frac_thresh, up2_sfc_coef, & ! intent(out)
               Skw_max_mag, xp3_coef_base, xp3_coef_slope, & ! intent(out)
               altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, & ! intent(out)
               C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, & ! intent(out)
               C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, & ! intent(out)
               C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, & ! intent(out)
               C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, & ! intent(out)
               Cx_min, Cx_max, Richardson_num_min, Richardson_num_max ) ! intent(out)


    ! It was decided after some experimentation, that the best
    ! way to produce grid independent results is to set lmin to be
    ! some fixed value. -dschanen 21 May 2007
    !lmin = lmin_coef * deltaz  ! Old
    lmin = lmin_coef * lmin_deltaz ! New fixed value

    ! ### Adjust Constant Diffusivity Coefficients Based On Grid Spacing ###
    call adj_low_res_nu &
           ( nzmax, grid_type, deltaz,  & ! Intent(in)
             momentum_heights, thermodynamic_heights, & ! Intent(in)
             l_prescribed_avg_deltaz, mult_coef, &  ! Intent(in)
             nu1, nu2, nu6, nu8, nu9, nu10, nu_hm, &  ! Intent(in)
             nu_vert_res_dep )  ! Intent(out)

    if ( beta < zero .or. beta > three ) then

       ! Constraints on beta
       write(fstderr,*) "beta = ", beta
       write(fstderr,*) "beta cannot be < 0 or > 3"
       err_code = clubb_fatal_error

    endif ! beta < 0 or beta > 3

    if ( slope_coef_spread_DG_means_w <= zero ) then

       ! Constraint on slope_coef_spread_DG_means_w
       write(fstderr,*) "slope_coef_spread_DG_means_w = ", &
                        slope_coef_spread_DG_means_w
       write(fstderr,*) "slope_coef_spread_DG_means_w cannot be <= 0"
       err_code = clubb_fatal_error

    endif ! slope_coef_spread_DG_means_w <= 0

    if ( pdf_component_stdev_factor_w <= zero ) then

       ! Constraint on pdf_component_stdev_factor_w
       write(fstderr,*) "pdf_component_stdev_factor_w = ", &
                        pdf_component_stdev_factor_w
       write(fstderr,*) "pdf_component_stdev_factor_w cannot be <= 0"
       err_code = clubb_fatal_error

    endif ! pdf_component_stdev_factor_w <= 0

    if ( coef_spread_DG_means_rt < zero &
         .or. coef_spread_DG_means_rt >= one ) then

       ! Constraint on coef_spread_DG_means_rt
       write(fstderr,*) "coef_spread_DG_means_rt = ", coef_spread_DG_means_rt
       write(fstderr,*) "coef_spread_DG_means_rt cannot be < 0 or >= 1"
       err_code = clubb_fatal_error

    endif ! coef_spread_DG_means_rt < 0 or coef_spread_DG_means_rt >= 1

    if ( coef_spread_DG_means_thl < zero &
         .or. coef_spread_DG_means_thl >= one ) then

       ! Constraint on coef_spread_DG_means_thl
       write(fstderr,*) "coef_spread_DG_means_thl = ", coef_spread_DG_means_thl
       write(fstderr,*) "coef_spread_DG_means_thl cannot be < 0 or >= 1"
       err_code = clubb_fatal_error

    endif ! coef_spread_DG_means_thl < 0 or coef_spread_DG_means_thl >= 1

    if ( omicron <= zero .or. omicron > one ) then

       ! Constraints on omicron
       write(fstderr,*) "omicron = ", omicron
       write(fstderr,*) "omicron cannot be <= 0 or > 1"
       err_code = clubb_fatal_error

    endif ! omicron <= 0 or omicron > 1

    if ( zeta_vrnce_rat <= -one ) then

       ! Constraints on zeta_vrnce_rat
       write(fstderr,*) "zeta_vrnce_rat = ", zeta_vrnce_rat
       write(fstderr,*) "zeta_vrnce_rat cannot be <= -1"
       err_code = clubb_fatal_error

    endif ! zeta_vrnce_rat <= -1

    if ( upsilon_precip_frac_rat < zero &
         .or. upsilon_precip_frac_rat > one ) then

       ! Constraints on upsilon_precip_frac_rat
       write(fstderr,*) "upsilon_precip_frac_rat = ", upsilon_precip_frac_rat
       write(fstderr,*) "upsilon_precip_frac_rat cannot be < 0 or > 1"
       err_code = clubb_fatal_error

    endif ! upsilon_precip_frac_rat < 0 or upsilon_precip_frac_rat > 1

    if ( mu < zero ) then

       ! Constraints on entrainment rate, mu.
       write(fstderr,*) "mu = ", mu
       write(fstderr,*) "mu cannot be < 0"
       err_code = clubb_fatal_error

    endif ! mu < 0.0

    if ( lmin < 1.0_core_rknd ) then

       ! Constraints on mixing length
       write(fstderr,*) "lmin = ", lmin
       write(fstderr,*) "lmin is < 1.0_core_rknd"
       err_code = clubb_fatal_error

    endif ! lmin < 1.0

     ! The C6rt parameters must be set equal to the C6thl parameters.
     ! Otherwise, the wpthlp pr1 term will be calculated inconsistently.

     if ( abs(C6rt - C6thl) > abs(C6rt + C6thl) / 2 * eps ) then
        write(fstderr,*) "C6rt = ", C6rt
        write(fstderr,*) "C6thl = ", C6thl
        write(fstderr,*) "C6rt and C6thl must be equal."
        err_code = clubb_fatal_error
     endif ! C6rt /= C6thl

     if ( abs(C6rtb - C6thlb) > abs(C6rtb + C6thlb) / 2 * eps ) then
        write(fstderr,*) "C6rtb = ", C6rtb
        write(fstderr,*) "C6thlb = ", C6thlb
        write(fstderr,*) "C6rtb and C6thlb must be equal."
        err_code = clubb_fatal_error
     endif ! C6rtb /= C6thlb

     if ( abs(C6rtc - C6thlc) > abs(C6rtc + C6thlc) / 2 * eps ) then
        write(fstderr,*) "C6rtc = ", C6rtc
        write(fstderr,*) "C6thlc = ", C6thlc
        write(fstderr,*) "C6rtc and C6thlc must be equal."
        err_code = clubb_fatal_error
     endif ! C6rtc /= C6thlc

     if ( abs(C6rt_Lscale0 - C6thl_Lscale0) > abs(C6rt_Lscale0 + C6thl_Lscale0) / 2 * eps ) then
        write(fstderr,*) "C6rt_Lscale0 = ", C6rt_Lscale0
        write(fstderr,*) "C6thl_Lscale0 = ", C6thl_Lscale0
        write(fstderr,*) "C6rt_Lscale0 and C6thl_Lscale0 must be equal."
        err_code = clubb_fatal_error
     endif ! C6rt_Lscale0 /= C6thl_Lscale0




    if ( C1 < zero ) then
        write(fstderr,*) "C1 = ", C1
        write(fstderr,*) "C1 must satisfy 0.0 <= C1"
        err_code = clubb_fatal_error
    end if

    if ( C7 > one .or. C7 < zero ) then
        write(fstderr,*) "C7 = ", C7
        write(fstderr,*) "C7 must satisfy 0.0 <= C7 <= 1.0"
        err_code = clubb_fatal_error
    end if

    if ( C7b > one .or. C7b < zero ) then
        write(fstderr,*) "C7b = ", C7b
        write(fstderr,*) "C7b must satisfy 0.0 <= C7b <= 1.0"
        err_code = clubb_fatal_error
    end if

    if ( C11 > one .or. C11 < zero ) then
        write(fstderr,*) "C11 = ", C11
        write(fstderr,*) "C11 must satisfy 0.0 <= C11 <= 1.0"
        err_code = clubb_fatal_error
    end if

    if ( C11b > one .or. C11b < zero ) then
        write(fstderr,*) "C11b = ", C11b
        write(fstderr,*) "C11b must satisfy 0.0 <= C11b <= 1.0"
        err_code = clubb_fatal_error
    end if

    if ( C_wp2_splat < zero ) then
        write(fstderr,*) "C_wp2_splat = ", C_wp2_splat
        write(fstderr,*) "C_wp2_splat must satisfy C_wp2_splat >= 0"
        err_code = clubb_fatal_error
    end if
    
    err_code_out = err_code

    return

  end subroutine setup_parameters

  !=============================================================================
  subroutine adj_low_res_nu &
               ( nzmax, grid_type, deltaz, & ! Intent(in)
                 momentum_heights, thermodynamic_heights, & ! Intent(in)
                 l_prescribed_avg_deltaz, mult_coef, &  ! Intent(in)
                 nu1, nu2, nu6, nu8, nu9, nu10, nu_hm, & ! Intent(out)
                 nu_vert_res_dep )  ! Intent(out)

    ! Description:
    !   Adjust the values of background eddy diffusivity based on
    !   vertical grid spacing.
    !   This code was made into a public subroutine so that it may be
    !   called multiple times per model run in scenarios where grid
    !   altitudes, and hence average grid spacing, change through space
    !   and/or time.  This occurs, for example, when CLUBB is
    !   implemented in WRF.  --ldgrant Jul 2010
    !----------------------------------------------------------------------

    use constants_clubb, only: &
        fstderr ! Constant(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Flag for adjusting the values of the constant background eddy diffusivity
    ! coefficients based on the average vertical grid spacing.  If this flag is
    ! turned off, the values of the various nu coefficients will remain as they
    ! are declared in the tunable_parameters.in file.
    logical, parameter :: l_adj_low_res_nu = .true.

    ! The size of the average vertical grid spacing that serves as a threshold
    ! for when to increase the size of the background eddy diffusivity
    ! coefficients (nus) by a certain factor above what the background
    ! coefficients are specified to be in tunable_parameters.in.  At any average
    ! grid spacing at or below this value, the values of the background
    ! diffusivities remain the same.  However, at any average vertical grid
    ! spacing above this value, the values of the background eddy diffusivities
    ! are increased.  Traditionally, the threshold grid spacing has been set to
    ! 40.0 meters.  This is only relevant if l_adj_low_res_nu is turned on.
    real( kind = core_rknd ), parameter :: &
      grid_spacing_thresh = 40.0_core_rknd  ! grid spacing threshold  [m]

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

    ! Local Variables
    real( kind = core_rknd ) :: avg_deltaz  ! Average grid box height   [m]

    ! The factor by which to multiply the coefficients of background eddy
    ! diffusivity if the grid spacing threshold is exceeded and l_adj_low_res_nu
    ! is turned on.
    real( kind = core_rknd ) :: &
      mult_factor_zt, &  ! Uses gr%dzt for nu values on zt levels
      mult_factor_zm     ! Uses gr%dzm for nu values on zm levels

    !--------------- Begin code -------------------------

    ! Flag for adjusting the values of the constant diffusivity coefficients
    ! based on the grid spacing.  If this flag is turned off, the values of the
    ! various nu coefficients will remain as they are declared in the
    ! parameters.in file.
    if ( l_adj_low_res_nu ) then

      ! ### Adjust Constant Diffusivity Coefficients Based On Grid Spacing ###

      ! All of the background coefficients of eddy diffusivity, as well as the
      ! constant coefficient for 4th-order hyper-diffusion, must be adjusted
      ! based on the size of the grid spacing.  For a case that uses an
      ! evenly-spaced grid, the adjustment is based on the constant grid
      ! spacing deltaz.  For a case that uses a stretched grid, the adjustment
      ! is based on avg_deltaz, which is the average grid spacing over the
      ! vertical domain.
 
      if ( l_prescribed_avg_deltaz ) then
        
        avg_deltaz = deltaz

      else if ( grid_type == 3 ) then

        ! CLUBB is implemented in a host model, or is using grid_type = 3

        ! Find the average deltaz over the grid based on momentum level
        ! inputs.

        avg_deltaz  &
           = ( momentum_heights(nzmax) - momentum_heights(1) )  &
             / real( nzmax - 1, kind = core_rknd )

      else if ( grid_type == 1 ) then

        ! Evenly-spaced grid.

        avg_deltaz = deltaz

      else if ( grid_type == 2 ) then

        ! Stretched (unevenly-spaced) grid:  stretched thermodynamic level
        ! input.

        ! Find the average deltaz over the stretched grid based on
        ! thermodynamic level inputs.

        avg_deltaz  &
          = ( thermodynamic_heights(nzmax) - thermodynamic_heights(1) )  &
             / real( nzmax - 1, kind = core_rknd )
      else
        ! Eric Raut added to remove compiler warning. (Obviously, this value is not used)
        avg_deltaz = 0.0_core_rknd
        write(fstderr,*) "Invalid grid_type:", grid_type
        error stop "Fatal error"

      end if ! grid_type

      ! The nu's are chosen for deltaz <= 40 m. Looks like they must
      ! be adjusted for larger grid spacings (Vince Larson)

      ! Use a constant mult_factor so nu does not depend on grid spacing
      if( avg_deltaz > grid_spacing_thresh ) then
        mult_factor_zt = 1.0_core_rknd + mult_coef * log( avg_deltaz / grid_spacing_thresh )
        mult_factor_zm = mult_factor_zt
      else
        mult_factor_zt = 1.0_core_rknd
        mult_factor_zm = 1.0_core_rknd
      end if

      !mult_factor = 1.0_core_rknd + mult_coef * log( avg_deltaz / grid_spacing_thresh )
      nu_vert_res_dep%nu1   =  nu1 * mult_factor_zm
      nu_vert_res_dep%nu2   =  nu2 * mult_factor_zm
      nu_vert_res_dep%nu6   =  nu6 * mult_factor_zm
      nu_vert_res_dep%nu8   =  nu8 * mult_factor_zt
      nu_vert_res_dep%nu9   =  nu9 * mult_factor_zm
      nu_vert_res_dep%nu10  =  nu10 * mult_factor_zt !We're unsure of the grid
      nu_vert_res_dep%nu_hm =  nu_hm * mult_factor_zt

    else ! nu values are not adjusted

      nu_vert_res_dep%nu1   =  nu1
      nu_vert_res_dep%nu2   =  nu2
      nu_vert_res_dep%nu6   =  nu6
      nu_vert_res_dep%nu8   =  nu8
      nu_vert_res_dep%nu9   =  nu9
      nu_vert_res_dep%nu10  =  nu10
      nu_vert_res_dep%nu_hm =  nu_hm

    end if  ! l_adj_low_res_nu

    return
  end subroutine adj_low_res_nu

#ifdef E3SM
  !=============================================================================
  subroutine clubb_param_readnl(filename)

    ! Description:
    ! Read clubb tunable parameters from namelist to override preset values
    ! To be called by clubb_readnl in clubb_intr.F90
    ! Author: Wuyin Lin

    !-----------------------------------------------------------------------
    use spmd_utils,      only: masterproc
    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use cam_abortutils,  only: endrun
    use mpishorthand, only: mpicom, mpir8

    implicit none

    ! Input variables

    character(len=*), intent(in) :: filename

    namelist /clubb_param_nl/      &
    clubb_C1,                           &
    clubb_C1b,                          &
    clubb_C1c,                          &
    clubb_C2rt,                         &
    clubb_C2thl,                        &
    clubb_C2rtthl,                      &
    clubb_C4,                           &
    clubb_C_uu_shr,                     &
    clubb_C_uu_buoy,                    &
    clubb_C6rt,                         &
    clubb_C6rtb,                        &
    clubb_C6rtc,                        &
    clubb_C6thl,                        &
    clubb_C6thlb,                       &
    clubb_C6thlc,                       &
    clubb_C7,                           &
    clubb_C7b,                          &
    clubb_C8,                           &
    clubb_C8b,                          &
    clubb_C11,                          &
    clubb_C11b,                         &
    clubb_C11c,                         &
    clubb_C14,                          &
    clubb_C_wp2_pr_dfsn,                &
    clubb_C_wp3_pr_tp,                  &
    clubb_C_wp3_pr_turb,                &
    clubb_C_wp3_pr_dfsn,                &
    clubb_beta,                         &
    clubb_gamma_coef,                   &
    clubb_gamma_coefb,                  &
    clubb_gamma_coefc,                  &
    clubb_pdf_component_stdev_factor_w, &
    clubb_mu,                           &
    clubb_c_K1,                         &
    clubb_nu1,                          &
    clubb_c_K2,                         &
    clubb_nu2,                          &
    clubb_c_K8,                         &
    clubb_nu8,                          &
    clubb_c_K9,                         &
    clubb_nu9,                          &
    clubb_c_K10,                        &
    clubb_c_K10h,                       &
    clubb_c_K_hmb,                      &
    clubb_wpxp_L_thresh,                &
    clubb_lmin_coef,                    &
    clubb_mult_coef,                    &
    clubb_Skw_denom_coef,               &
    clubb_up2_sfc_coef,                 &
    clubb_Skw_max_mag,                  &
    clubb_C_invrs_tau_bkgnd,            &
    clubb_C_invrs_tau_sfc,              &
    clubb_C_invrs_tau_shear,            &
    clubb_C_invrs_tau_N2,               &
    clubb_C_invrs_tau_N2_wp2,           &
    clubb_C_invrs_tau_N2_xp2,           &
    clubb_C_invrs_tau_N2_wpxp,          &
    clubb_C_invrs_tau_N2_clear_wp3,     &
    clubb_C_invrs_tau_wpxp_Ri,          &
    clubb_C_invrs_tau_wpxp_N2_thresh,   &
    clubb_C_wp2_splat,                  &
    clubb_xp3_coef_base,                &
    clubb_xp3_coef_slope,               &
    clubb_altitude_threshold,           &
    clubb_rtp2_clip_coef,               &
    clubb_Cx_min,                       &
    clubb_Cx_max,                       &
    clubb_Richardson_num_min,           &
    clubb_Richardson_num_max

    integer :: read_status
    integer :: iunit

    ! ---- Begin Code ----

    ! If clubb_tunables_nl present, read in to replace preset values
    ! This is made available for tuning 
     
    clubb_C1 = init_value
    clubb_C1b = init_value
    clubb_C1c = init_value
    clubb_C2rt = init_value
    clubb_C2thl = init_value
    clubb_C2rtthl = init_value
    clubb_C4 = init_value
    clubb_C_uu_shr = init_value
    clubb_C_uu_buoy = init_value
    clubb_C6rt = init_value
    clubb_C6rtb = init_value
    clubb_C6rtc = init_value
    clubb_C6thl = init_value
    clubb_C6thlb = init_value
    clubb_C6thlc = init_value
    clubb_C7 = init_value
    clubb_C7b = init_value
    clubb_C8 = init_value
    clubb_C8b = init_value
    clubb_C11 = init_value
    clubb_C11b = init_value
    clubb_C11c = init_value
    clubb_C14 = init_value
    clubb_C_wp2_pr_dfsn = init_value
    clubb_C_wp3_pr_tp = init_value
    clubb_C_wp3_pr_turb = init_value
    clubb_C_wp3_pr_dfsn = init_value
    clubb_beta = init_value
    clubb_gamma_coef = init_value
    clubb_gamma_coefb = init_value
    clubb_gamma_coefc = init_value
    clubb_pdf_component_stdev_factor_w = init_value
    clubb_mu = init_value
    clubb_c_K1 = init_value
    clubb_nu1 = init_value
    clubb_c_K2 = init_value
    clubb_nu2 = init_value
    clubb_c_K8 = init_value
    clubb_nu8 = init_value
    clubb_c_K9 = init_value
    clubb_nu9 = init_value
    clubb_c_K10 = init_value
    clubb_c_K10h = init_value
    clubb_c_K_hmb = init_value
    clubb_wpxp_L_thresh = init_value
    clubb_lmin_coef = init_value
    clubb_mult_coef = init_value
    clubb_Skw_denom_coef = init_value
    clubb_up2_sfc_coef = init_value
    clubb_Skw_max_mag = init_value
    clubb_C_invrs_tau_bkgnd = init_value
    clubb_C_invrs_tau_sfc = init_value
    clubb_C_invrs_tau_shear = init_value
    clubb_C_invrs_tau_N2 = init_value
    clubb_C_invrs_tau_N2_wp2 = init_value
    clubb_C_invrs_tau_N2_xp2 = init_value
    clubb_C_invrs_tau_N2_wpxp = init_value
    clubb_C_invrs_tau_N2_clear_wp3 = init_value
    clubb_C_invrs_tau_wpxp_Ri = init_value
    clubb_C_invrs_tau_wpxp_N2_thresh = init_value
    clubb_C_wp2_splat = init_value
    clubb_xp3_coef_base = init_value
    clubb_xp3_coef_slope = init_value
    clubb_altitude_threshold = init_value
    clubb_rtp2_clip_coef = init_value
    clubb_Cx_min = init_value
    clubb_Cx_max = init_value
    clubb_Richardson_num_min = init_value
    clubb_Richardson_num_max = init_value

    if (masterproc) then
      iunit = getunit()
      open( iunit, file=trim(filename), status='old' )
      call find_group_name(iunit, 'clubb_param_nl', status=read_status)
      if (read_status == 0) then
         read(unit=iunit, nml=clubb_param_nl,iostat=read_status)
         if (read_status /= 0) then
            call endrun('clubb_param_readnl:  error reading namelist')
         end if
      endif
      close(unit=iunit)
      call freeunit(iunit)
    end if
#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(clubb_C1,         1, mpir8,  0, mpicom)
   call mpibcast(clubb_C1b,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_C1c,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_C2rt,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_C2thl,      1, mpir8,  0, mpicom)
   call mpibcast(clubb_C2rtthl,    1, mpir8,  0, mpicom)
   call mpibcast(clubb_C4,         1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_uu_shr,   1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_uu_buoy,  1, mpir8,  0, mpicom)
   call mpibcast(clubb_C6rt,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_C6rtb,      1, mpir8,  0, mpicom)
   call mpibcast(clubb_C6rtc,      1, mpir8,  0, mpicom)
   call mpibcast(clubb_C6thl,      1, mpir8,  0, mpicom)
   call mpibcast(clubb_C6thlb,     1, mpir8,  0, mpicom)
   call mpibcast(clubb_C6thlc,     1, mpir8,  0, mpicom)
   call mpibcast(clubb_C7,         1, mpir8,  0, mpicom)
   call mpibcast(clubb_C7b,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_C8,         1, mpir8,  0, mpicom)
   call mpibcast(clubb_C8b,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_C11,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_C11b,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_C11c,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_C14,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_wp2_pr_dfsn, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_wp3_pr_tp,   1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_wp3_pr_turb, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_wp3_pr_dfsn, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_beta,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_gamma_coef, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_gamma_coefb,1, mpir8,  0, mpicom)
   call mpibcast(clubb_gamma_coefc,1, mpir8,  0, mpicom)
   call mpibcast(clubb_pdf_component_stdev_factor_w, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_mu,         1, mpir8,  0, mpicom)
   call mpibcast(clubb_c_K1,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_nu1,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_c_K2,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_nu2,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_c_K8,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_nu8,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_c_K9,       1, mpir8,  0, mpicom)
   call mpibcast(clubb_nu9,        1, mpir8,  0, mpicom)
   call mpibcast(clubb_c_K10,      1, mpir8,  0, mpicom)
   call mpibcast(clubb_c_K10h,     1, mpir8,  0, mpicom)
   call mpibcast(clubb_c_K_hmb,    1, mpir8,  0, mpicom)
   call mpibcast(clubb_wpxp_L_thresh, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_lmin_coef,  1, mpir8,  0, mpicom)
   call mpibcast(clubb_mult_coef,  1, mpir8,  0, mpicom)
   call mpibcast(clubb_Skw_denom_coef, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_up2_sfc_coef, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_Skw_max_mag, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_invrs_tau_bkgnd, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_invrs_tau_sfc, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_invrs_tau_shear, 1, mpir8,  0, mpicom) 
   call mpibcast(clubb_C_invrs_tau_N2, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_invrs_tau_N2_wp2, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_invrs_tau_N2_xp2, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_invrs_tau_N2_wpxp, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_C_invrs_tau_N2_clear_wp3, 1, mpir8, 0, mpicom)
   call mpibcast(clubb_C_invrs_tau_wpxp_Ri, 1, mpir8, 0, mpicom)
   call mpibcast(clubb_C_invrs_tau_wpxp_N2_thresh, 1, mpir8, 0, mpicom)
   call mpibcast(clubb_C_wp2_splat, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_xp3_coef_base, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_xp3_coef_slope, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_altitude_threshold, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_rtp2_clip_coef, 1, mpir8,  0, mpicom)
   call mpibcast(clubb_Cx_min, 1, mpir8, 0, mpicom)
   call mpibcast(clubb_Cx_max, 1, mpir8, 0, mpicom)
   call mpibcast(clubb_Richardson_num_min, 1, mpir8, 0, mpicom)
   call mpibcast(clubb_Richardson_num_max, 1, mpir8, 0, mpicom)
#endif


  end subroutine clubb_param_readnl
#endif /*E3SM*/
  !=============================================================================
  subroutine read_parameters( iunit, filename, &
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
                              Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, &
                              params )

    ! Description:
    ! Read a namelist containing the model parameters

    ! References:
    ! None
    !-----------------------------------------------------------------------
!    use constants_clubb, only: fstderr ! Constant

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
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(nparams) :: params

    ! Local variables
!    integer :: i

!    logical :: l_error

    ! Since we lack a devious way to do this just once, this namelist
    ! must be changed as well when a new parameter is added.
    namelist /clubb_params_nl/  & 
      C1, C1b, C1c, & 
      C2rt, C2thl, C2rtthl, C4, C_uu_shr, C_uu_buoy, & 
      C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & 
      C7, C7b, C7c, C8, C8b, C10, C11, C11b, C11c, & 
      C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, &
      C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, & 
      C6rt_Lscale0, C6thl_Lscale0, &
      C7_Lscale0, wpxp_L_thresh, c_K, c_K1, nu1, c_K2, nu2, & 
      c_K6, nu6, c_K8, nu8, c_K9, nu9, nu10, &
      c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, &
      slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, &
      coef_spread_DG_means_rt, coef_spread_DG_means_thl, &
      beta, gamma_coef, gamma_coefb, gamma_coefc, lmin_coef, &
      omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, &
      lambda0_stability_coef, mult_coef, taumin, taumax, mu, Lscale_mu_coef, &
      Lscale_pert_coef, alpha_corr, Skw_denom_coef, c_K10, c_K10h, &
      thlp2_rad_coef, thlp2_rad_cloud_frac_thresh, up2_sfc_coef, &
      Skw_max_mag, xp3_coef_base, xp3_coef_slope, altitude_threshold, &
      rtp2_clip_coef, C_invrs_tau_bkgnd, C_invrs_tau_sfc, &
      C_invrs_tau_shear, C_invrs_tau_N2, C_invrs_tau_N2_wp2, &
      C_invrs_tau_N2_xp2, C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, &
      C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, &
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max

    ! ---- Begin Code ----

    ! If the filename is empty, assume we're using a `working' set of
    ! parameters that are set statically here (handy for host models).
    ! Read the namelist
    if ( filename /= "" ) then
      ! Read the namelist
      open(unit=iunit, file=filename, status='old', action='read')

      read(unit=iunit, nml=clubb_params_nl)

      close(unit=iunit)

    end if

#ifdef E3SM
    if (clubb_C1 /= init_value) then
       C1 = clubb_C1
    endif
    if (clubb_C1b /= init_value) then
       C1b = clubb_C1b
    end if
    if (clubb_C1c /= init_value) then
       C1c = clubb_C1c
    end if
    ! if clubb_C2thl and clubb_C2rtthl not specified, continue to use C2thl=C2rt, C2rtthl = 1.3*C2rt
    ! to preserve existing compsets that have assumed so and only vary C2rt
    if (clubb_C2rt /= init_value) then
       C2rt = clubb_C2rt
       if (clubb_C2thl == init_value) C2thl = C2rt
       if (clubb_C2rtthl == init_value) C2rtthl = C2rt*2.0_core_rknd
    end if
    ! Allows C2thl and C2rtthl to vary separately
    if (clubb_C2thl /= init_value) C2thl = clubb_C2thl
    if (clubb_C2rtthl /= init_value) C2rtthl = clubb_C2rtthl
    if (clubb_C4 /= init_value) C4 = clubb_C4
    if (clubb_C_uu_shr /= init_value) C_uu_shr = clubb_C_uu_shr
    if (clubb_C_uu_buoy /= init_value) C_uu_buoy = clubb_C_uu_buoy
    if (clubb_C6rt /= init_value) then
       C6rt = clubb_C6rt
       if (clubb_C6thl == init_value) C6thl = C6rt
    end if
    if (clubb_C6rtb /= init_value) C6rtb = clubb_C6rtb
    if (clubb_C6rtc /= init_value) C6rtc = clubb_C6rtc
    if (clubb_C6thl /= init_value) C6thl = clubb_C6thl
    if (clubb_C6thlb /= init_value) C6thlb = clubb_C6thlb
    if (clubb_C6thlc /= init_value) C6thlc = clubb_C6thlc
    if (clubb_C7 /= init_value) C7 = clubb_C7
    if (clubb_C7b /= init_value) C7b = clubb_C7b
    if (clubb_C8 /= init_value) C8 = clubb_C8
    if (clubb_C8b /= init_value) C8b = clubb_C8b
    if (clubb_C11 /= init_value) C11 = clubb_C11
    if (clubb_C11b /= init_value) C11b = clubb_C11b
    if (clubb_C11c /= init_value) C11c = clubb_C11c
    if (clubb_C14 /= init_value) C14 = clubb_C14
    if (clubb_C_wp2_pr_dfsn /= init_value) C_wp2_pr_dfsn = clubb_C_wp2_pr_dfsn
    if (clubb_C_wp3_pr_tp /= init_value) C_wp3_pr_tp = clubb_C_wp3_pr_tp
    if (clubb_C_wp3_pr_turb /= init_value) C_wp3_pr_turb = clubb_C_wp3_pr_turb
    if (clubb_C_wp3_pr_dfsn /= init_value) C_wp3_pr_dfsn = clubb_C_wp3_pr_dfsn
    if (clubb_beta /= init_value) beta = clubb_beta
    ! if clubb_gamma_coefb not specified, continue to use gamma_coefb=gamma_coef
    ! to preserve existing compsets that have assumed so  and only vary gamma_coef
    if (clubb_gamma_coef /= init_value) then
       gamma_coef = clubb_gamma_coef
       if (clubb_gamma_coefb == init_value) gamma_coefb = gamma_coef
    end if
    ! Allows gamma_coefb to vary separately
    if (clubb_gamma_coefb /= init_value) gamma_coefb = clubb_gamma_coefb
    if (clubb_gamma_coefc /= init_value) gamma_coefc = clubb_gamma_coefc
    if (clubb_pdf_component_stdev_factor_w /= init_value) &
       pdf_component_stdev_factor_w = clubb_pdf_component_stdev_factor_w
    if (clubb_mu /= init_value) mu = clubb_mu
    if (clubb_c_K1 /= init_value) c_K1 = clubb_c_K1
    if (clubb_nu1 /= init_value) nu1 = clubb_nu1
    if (clubb_c_K2 /= init_value) c_K2 = clubb_c_K2
    if (clubb_nu2 /= init_value) nu2 = clubb_nu2
    if (clubb_c_K8 /= init_value) c_K8 = clubb_c_K8
    if (clubb_nu8 /= init_value) nu8 = clubb_nu8
    if (clubb_c_K9 /= init_value) c_K9 = clubb_c_K9
    if (clubb_nu9 /= init_value) nu9 = clubb_nu9
    if (clubb_c_K10 /= init_value) c_K10 = clubb_c_K10
    if (clubb_c_K10h /= init_value) c_K10h = clubb_c_K10h
    if (clubb_c_K_hmb /= init_value) c_K_hmb = clubb_c_K_hmb
    if (clubb_wpxp_L_thresh /= init_value) wpxp_L_thresh = clubb_wpxp_L_thresh
    if (clubb_lmin_coef /= init_value) lmin_coef = clubb_lmin_coef
    if (clubb_mult_coef /= init_value) mult_coef = clubb_mult_coef
    if (clubb_Skw_denom_coef /= init_value) &
       Skw_denom_coef = clubb_Skw_denom_coef
    if (clubb_up2_sfc_coef /= init_value) &
       up2_sfc_coef = clubb_up2_sfc_coef
    if (clubb_Skw_max_mag /= init_value) &
       Skw_max_mag = clubb_Skw_max_mag
    if (clubb_C_invrs_tau_bkgnd /= init_value) &
       C_invrs_tau_bkgnd = clubb_C_invrs_tau_bkgnd
    if (clubb_C_invrs_tau_sfc /= init_value) &
       C_invrs_tau_sfc = clubb_C_invrs_tau_sfc
    if (clubb_C_invrs_tau_shear /= init_value) &
       C_invrs_tau_shear = clubb_C_invrs_tau_shear
    if (clubb_C_invrs_tau_N2 /= init_value) &
       C_invrs_tau_N2 = clubb_C_invrs_tau_N2
    if (clubb_C_invrs_tau_N2_wp2 /= init_value) &
       C_invrs_tau_N2_wp2 = clubb_C_invrs_tau_N2_wp2
    if (clubb_C_invrs_tau_N2_xp2 /= init_value) &
       C_invrs_tau_N2_xp2 = clubb_C_invrs_tau_N2_xp2
    if (clubb_C_invrs_tau_N2_wpxp /= init_value) &
       C_invrs_tau_N2_wpxp = clubb_C_invrs_tau_N2_wpxp
    if (clubb_C_invrs_tau_N2_clear_wp3 /= init_value) &
       C_invrs_tau_N2_clear_wp3 = clubb_C_invrs_tau_N2_clear_wp3
    if (clubb_C_invrs_tau_wpxp_Ri /= init_value) &
       C_invrs_tau_wpxp_Ri = clubb_C_invrs_tau_wpxp_Ri
    if (clubb_C_invrs_tau_wpxp_N2_thresh /= init_value) &
       C_invrs_tau_wpxp_N2_thresh = clubb_C_invrs_tau_wpxp_N2_thresh
    if (clubb_C_wp2_splat  /= init_value ) C_wp2_splat = clubb_C_wp2_splat
    if (clubb_xp3_coef_base /= init_value) xp3_coef_base = clubb_xp3_coef_base
    if (clubb_xp3_coef_slope /= init_value) &
       xp3_coef_slope = clubb_xp3_coef_slope
    if (clubb_altitude_threshold /= init_value) &
       altitude_threshold = clubb_altitude_threshold
    if (clubb_rtp2_clip_coef /= init_value) &
       rtp2_clip_coef = clubb_rtp2_clip_coef
    if (clubb_Cx_min /= init_value) Cx_min = clubb_Cx_min
    if (clubb_Cx_max /= init_value) Cx_max = clubb_Cx_max
    if (clubb_Richardson_num_min /= init_value) &
       Richardson_num_min = clubb_Richardson_num_min
    if (clubb_Richardson_num_max /= init_value) &
       Richardson_num_max = clubb_Richardson_num_max
#endif /*E3SM*/

    ! Put the variables in the output array
    call pack_parameters &
             ( C1, C1b, C1c, C2rt, C2thl, C2rtthl, & ! intent(in)
               C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, & ! intent(in)
               C7, C7b, C7c, C8, C8b, C10, & ! intent(in)
               C11, C11b, C11c, C12, C13, C14, C_wp2_pr_dfsn, C_wp3_pr_tp, & ! intent(in)
               C_wp3_pr_turb, C_wp3_pr_dfsn, C_wp2_splat, & ! intent(in)
               C6rt_Lscale0, C6thl_Lscale0, C7_Lscale0, wpxp_L_thresh, & ! intent(in)
               c_K, c_K1, nu1, c_K2, nu2, c_K6, nu6, c_K8, nu8, & ! intent(in)
               c_K9, nu9, nu10, c_K_hm, c_K_hmb, K_hm_min_coef, nu_hm, & ! intent(in)
               slope_coef_spread_DG_means_w, pdf_component_stdev_factor_w, & ! intent(in)
               coef_spread_DG_means_rt, coef_spread_DG_means_thl, & ! intent(in)
               gamma_coef, gamma_coefb, gamma_coefc, mu, beta, lmin_coef, & ! intent(in)
               omicron, zeta_vrnce_rat, upsilon_precip_frac_rat, & ! intent(in)
               lambda0_stability_coef, mult_coef, taumin, taumax, & ! intent(in)
               Lscale_mu_coef, Lscale_pert_coef, alpha_corr, & ! intent(in)
               Skw_denom_coef, c_K10, c_K10h, thlp2_rad_coef, & ! intent(in)
               thlp2_rad_cloud_frac_thresh, up2_sfc_coef, & ! intent(in)
               Skw_max_mag, xp3_coef_base, xp3_coef_slope, & ! intent(in)
               altitude_threshold, rtp2_clip_coef, C_invrs_tau_bkgnd, & ! intent(in)
               C_invrs_tau_sfc, C_invrs_tau_shear, C_invrs_tau_N2, & ! intent(in)
               C_invrs_tau_N2_wp2, C_invrs_tau_N2_xp2, &   ! intent(in)
               C_invrs_tau_N2_wpxp, C_invrs_tau_N2_clear_wp3, & ! intent(in)
               C_invrs_tau_wpxp_Ri, C_invrs_tau_wpxp_N2_thresh, & ! intent(in)
               Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, & ! intent(in)
               params ) ! intent(out)

!    l_error = .false.

!    This error check is currently broken since we are not initializing the
!    parameters to -999 ( = init_value).
!    do i = 1, nparams
!      if ( abs(params(i)-init_value) < abs(params(i)+init_value) / 2 * eps) then
!        write(fstderr,*) "Tuning parameter "//trim( params_list(i) )// &
!          " was missing from "//trim( filename )
!        l_error = .true.
!      end if
!    end do

!    if ( l_error ) error stop "Fatal error."

    return

  end subroutine read_parameters

  !=============================================================================
  subroutine read_param_minmax & 
           ( iunit, filename, nindex, params_minmax, ndim )

    ! Description:
    ! Read a namelist containing the amount to vary model parameters.
    ! Used by the downhill simplex / simulated annealing algorithm.

    ! References:
    ! None
    !-----------------------------------------------------------------------
    use constants_clubb, only: fstderr ! Constant

    use clubb_precision, only: &
        core_rknd ! Variable(s)

   use parameter_indices, only: &
        iC1,  & ! Variable(s)
        iC1b, &
        iC1c, &
        iC2rt, &
        iC2thl, &
        iC2rtthl, &
        iC4, &
        iC_uu_shr, &
        iC_uu_buoy, &
        iC6rt, &
        iC6rtb, &
        iC6rtc, &
        iC6thl, &
        iC6thlb, &
        iC6thlc, &
        iC7, &
        iC7b, &
        iC7c, &
        iC8, &
        iC8b, &
        iC10, &
        iC11, &
        iC11b, &
        iC11c, &
        iC12, &
        iC13, &
        iC14, &
        iC_wp2_pr_dfsn, &
        iC_wp3_pr_tp, &
        iC_wp3_pr_turb, &
        iC_wp3_pr_dfsn, &
        iC_wp2_splat

    use parameter_indices, only: &
        iC6rt_Lscale0, &
        iC6thl_Lscale0, &
        iC7_Lscale0, &
        iwpxp_L_thresh

    use parameter_indices, only: &
        ic_K,  &
        ic_K1, &
        inu1, &
        ic_K2, &
        inu2, &
        ic_K6, &
        inu6, &
        ic_K8, &
        inu8, &
        ic_K9, &
        inu9, &
        inu10, &
        ic_K_hm, &
        ic_K_hmb, &
        iK_hm_min_coef, &
        inu_hm, &
        islope_coef_spread_DG_means_w, &
        ipdf_component_stdev_factor_w, &
        icoef_spread_DG_means_rt, &
        icoef_spread_DG_means_thl, &
        igamma_coef, &
        igamma_coefb, &
        igamma_coefc, &
        imu, &
        ibeta, &
        ilmin_coef, &
        iomicron, &
        izeta_vrnce_rat, &
        iupsilon_precip_frac_rat, &
        ilambda0_stability_coef, &
        imult_coef, &
        itaumin, &
        itaumax, &
        iLscale_mu_coef, &
        iLscale_pert_coef, &
        ialpha_corr, &
        iSkw_denom_coef, &
        ic_K10, &
        ic_K10h, &
        ithlp2_rad_coef, &
        ithlp2_rad_cloud_frac_thresh, &
        iup2_sfc_coef, &
        iSkw_max_mag, &
        ixp3_coef_base, &
        ixp3_coef_slope, &
        ialtitude_threshold, &
        irtp2_clip_coef, &
        iC_invrs_tau_bkgnd, &
        iC_invrs_tau_sfc, &
        iC_invrs_tau_shear, &
        iC_invrs_tau_N2, &
        iC_invrs_tau_N2_wp2, &
        iC_invrs_tau_N2_xp2, &
        iC_invrs_tau_N2_wpxp, &
        iC_invrs_tau_N2_clear_wp3, &
        iC_invrs_tau_wpxp_Ri, &
        iC_invrs_tau_wpxp_N2_thresh, &
        iCx_min, &
        iCx_max, &
        iRichardson_num_min, &
        iRichardson_num_max

    implicit none

    ! Input variables
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Output variables

    ! An array of array indices (i.e. which elements of the array `params'
    ! are contained within the simplex and the max variable)
    integer, intent(out), dimension(nparams) :: nindex

    real( kind = core_rknd ), intent(out), dimension(2,nparams) ::  & 
      params_minmax  ! Amount to vary the parameter in the initial simplex

    integer, intent(out) :: &
      ndim  ! Number of variables, e.g. rcm, to be tuned. Dimension of the init simplex

    ! Local variables
    integer :: i

    ! Amount to change each parameter for the initial simplex
    ! This MUST be changed to match the clubb_params_nl namelist if parameters are added.
    real ( kind = core_rknd ), dimension(2) :: &
      C1_minmax, C1b_minmax, C1c_minmax, C2rt_minmax, C2thl_minmax, C2rtthl_minmax, C4_minmax, &
      C_uu_shr_minmax, C_uu_buoy_minmax, C6rt_minmax, C6rtb_minmax, C6rtc_minmax, C6thl_minmax, &
      C6thlb_minmax, C6thlc_minmax, C7_minmax, C7b_minmax, C7c_minmax, C8_minmax, C8b_minmax, &
      C10_minmax, C11_minmax, C11b_minmax, C11c_minmax, C12_minmax, C13_minmax, C14_minmax, &
      C_wp2_pr_dfsn_minmax, C_wp3_pr_tp_minmax, C_wp3_pr_turb_minmax, C_wp3_pr_dfsn_minmax, &
      C_wp2_splat_minmax, C6rt_Lscale0_minmax, C6thl_Lscale0_minmax, C7_Lscale0_minmax, &
      wpxp_L_thresh_minmax, c_K_minmax, c_K1_minmax, nu1_minmax, c_K2_minmax, nu2_minmax, &
      c_K6_minmax, nu6_minmax, c_K8_minmax, nu8_minmax, c_K9_minmax, nu9_minmax, nu10_minmax, &
      c_K_hm_minmax, c_K_hmb_minmax, K_hm_min_coef_minmax, nu_hm_minmax, &
      slope_coef_spread_DG_means_w_minmax, pdf_component_stdev_factor_w_minmax, &
      coef_spread_DG_means_rt_minmax, coef_spread_DG_means_thl_minmax, &
      beta_minmax, gamma_coef_minmax, gamma_coefb_minmax, gamma_coefc_minmax, lmin_coef_minmax, &
      omicron_minmax, zeta_vrnce_rat_minmax, upsilon_precip_frac_rat_minmax, &
      lambda0_stability_coef_minmax, mult_coef_minmax, taumin_minmax, taumax_minmax, &
      mu_minmax, Lscale_mu_coef_minmax, Lscale_pert_coef_minmax, alpha_corr_minmax, &
      Skw_denom_coef_minmax, c_K10_minmax, c_K10h_minmax, thlp2_rad_coef_minmax, &
      thlp2_rad_cloud_frac_thresh_minmax, up2_sfc_coef_minmax, Skw_max_mag_minmax, &
      xp3_coef_base_minmax, xp3_coef_slope_minmax, altitude_threshold_minmax, &
      rtp2_clip_coef_minmax, C_invrs_tau_bkgnd_minmax, C_invrs_tau_sfc_minmax, &
      C_invrs_tau_shear_minmax, C_invrs_tau_N2_minmax, C_invrs_tau_N2_wp2_minmax, &
      C_invrs_tau_N2_xp2_minmax, C_invrs_tau_N2_wpxp_minmax, C_invrs_tau_N2_clear_wp3_minmax, &
      C_invrs_tau_wpxp_Ri_minmax, C_invrs_tau_wpxp_N2_thresh_minmax, Cx_min_minmax, &
      Cx_max_minmax, Richardson_num_min_minmax, Richardson_num_max_minmax

    namelist /init_minmax/  & 
      C1_minmax, C1b_minmax, C1c_minmax, C2rt_minmax, C2thl_minmax, C2rtthl_minmax, C4_minmax, &
      C_uu_shr_minmax, C_uu_buoy_minmax, C6rt_minmax, C6rtb_minmax, C6rtc_minmax, C6thl_minmax, &
      C6thlb_minmax, C6thlc_minmax, C7_minmax, C7b_minmax, C7c_minmax, C8_minmax, C8b_minmax, &
      C10_minmax, C11_minmax, C11b_minmax, C11c_minmax, C12_minmax, C13_minmax, C14_minmax, &
      C_wp2_pr_dfsn_minmax, C_wp3_pr_tp_minmax, C_wp3_pr_turb_minmax, C_wp3_pr_dfsn_minmax, &
      C_wp2_splat_minmax, C6rt_Lscale0_minmax, C6thl_Lscale0_minmax, C7_Lscale0_minmax, &
      wpxp_L_thresh_minmax, c_K_minmax, c_K1_minmax, nu1_minmax, c_K2_minmax, nu2_minmax, &
      c_K6_minmax, nu6_minmax, c_K8_minmax, nu8_minmax, c_K9_minmax, nu9_minmax, nu10_minmax, &
      c_K_hm_minmax, c_K_hmb_minmax, K_hm_min_coef_minmax, nu_hm_minmax, &
      slope_coef_spread_DG_means_w_minmax, pdf_component_stdev_factor_w_minmax, &
      coef_spread_DG_means_rt_minmax, coef_spread_DG_means_thl_minmax, &
      beta_minmax, gamma_coef_minmax, gamma_coefb_minmax, gamma_coefc_minmax, lmin_coef_minmax, &
      omicron_minmax, zeta_vrnce_rat_minmax, upsilon_precip_frac_rat_minmax, &
      lambda0_stability_coef_minmax, mult_coef_minmax, taumin_minmax, taumax_minmax, &
      mu_minmax, Lscale_mu_coef_minmax, Lscale_pert_coef_minmax, alpha_corr_minmax, &
      Skw_denom_coef_minmax, c_K10_minmax, c_K10h_minmax, thlp2_rad_coef_minmax, &
      thlp2_rad_cloud_frac_thresh_minmax, up2_sfc_coef_minmax, Skw_max_mag_minmax, &
      xp3_coef_base_minmax, xp3_coef_slope_minmax, altitude_threshold_minmax, &
      rtp2_clip_coef_minmax, C_invrs_tau_bkgnd_minmax, C_invrs_tau_sfc_minmax, &
      C_invrs_tau_shear_minmax, C_invrs_tau_N2_minmax, C_invrs_tau_N2_wp2_minmax, &
      C_invrs_tau_N2_xp2_minmax, C_invrs_tau_N2_wpxp_minmax, C_invrs_tau_N2_clear_wp3_minmax, &
      C_invrs_tau_wpxp_Ri_minmax, C_invrs_tau_wpxp_N2_thresh_minmax, Cx_min_minmax, &
      Cx_max_minmax, Richardson_num_min_minmax, Richardson_num_max_minmax


! ----- Begin code -------------

    ! Read the namelist
    open(unit=iunit, file=filename, status='old', action='read')
    read(unit=iunit, nml=init_minmax)
    close(unit=iunit)

    ! Put min/max values into params_minmax output array
    params_minmax(:,iC1) = C1_minmax
    params_minmax(:,iC1b) = C1b_minmax
    params_minmax(:,iC1c) = C1c_minmax
    params_minmax(:,iC2rt) = C2rt_minmax
    params_minmax(:,iC2thl) = C2thl_minmax
    params_minmax(:,iC2rtthl) = C2rtthl_minmax
    params_minmax(:,iC4) = C4_minmax
    params_minmax(:,iC_uu_shr) = C_uu_shr_minmax
    params_minmax(:,iC_uu_buoy) = C_uu_buoy_minmax
    params_minmax(:,iC6rt) = C6rt_minmax
    params_minmax(:,iC6rtb) = C6rtb_minmax
    params_minmax(:,iC6rtc) = C6rtc_minmax
    params_minmax(:,iC6thl) = C6thl_minmax
    params_minmax(:,iC6thlb) = C6thlb_minmax
    params_minmax(:,iC6thlc) = C6thlc_minmax
    params_minmax(:,iC7) = C7_minmax
    params_minmax(:,iC7b) = C7b_minmax
    params_minmax(:,iC7c) = C7c_minmax
    params_minmax(:,iC8) = C8_minmax
    params_minmax(:,iC8b) = C8b_minmax
    params_minmax(:,iC10) = C10_minmax
    params_minmax(:,iC11) = C11_minmax
    params_minmax(:,iC11b) = C11b_minmax
    params_minmax(:,iC11c) = C11c_minmax
    params_minmax(:,iC12) = C12_minmax
    params_minmax(:,iC13) = C13_minmax
    params_minmax(:,iC14) = C14_minmax
    params_minmax(:,iC_wp2_pr_dfsn) = C_wp2_pr_dfsn_minmax
    params_minmax(:,iC_wp3_pr_tp) = C_wp3_pr_tp_minmax
    params_minmax(:,iC_wp3_pr_turb) = C_wp3_pr_turb_minmax
    params_minmax(:,iC_wp3_pr_dfsn) = C_wp3_pr_dfsn_minmax
    params_minmax(:,iC_wp2_splat) = C_wp2_splat_minmax
    params_minmax(:,iC6rt_Lscale0) = C6rt_Lscale0_minmax
    params_minmax(:,iC6thl_Lscale0) = C6thl_Lscale0_minmax
    params_minmax(:,iC7_Lscale0) = C7_Lscale0_minmax
    params_minmax(:,iwpxp_L_thresh) = wpxp_L_thresh_minmax
    params_minmax(:,ic_K) = c_K_minmax
    params_minmax(:,ic_K1) = c_K1_minmax
    params_minmax(:,inu1) = nu1_minmax
    params_minmax(:,ic_K2) = c_K2_minmax
    params_minmax(:,inu2) = nu2_minmax
    params_minmax(:,ic_K6) = c_K6_minmax
    params_minmax(:,inu6) = nu6_minmax
    params_minmax(:,ic_K8) = c_K8_minmax
    params_minmax(:,inu8) = nu8_minmax
    params_minmax(:,ic_K9) = c_K9_minmax
    params_minmax(:,inu9) = nu9_minmax
    params_minmax(:,inu10) = nu10_minmax
    params_minmax(:,ic_K_hm) = c_K_hm_minmax
    params_minmax(:,ic_K_hmb) = c_K_hmb_minmax
    params_minmax(:,iK_hm_min_coef) = K_hm_min_coef_minmax
    params_minmax(:,inu_hm) = nu_hm_minmax
    params_minmax(:,islope_coef_spread_DG_means_w) = slope_coef_spread_DG_means_w_minmax
    params_minmax(:,ipdf_component_stdev_factor_w) = pdf_component_stdev_factor_w_minmax
    params_minmax(:,icoef_spread_DG_means_rt) = coef_spread_DG_means_rt_minmax
    params_minmax(:,icoef_spread_DG_means_thl) = coef_spread_DG_means_thl_minmax
    params_minmax(:,igamma_coef) = gamma_coef_minmax
    params_minmax(:,igamma_coefb) = gamma_coefb_minmax
    params_minmax(:,igamma_coefc) = gamma_coefc_minmax
    params_minmax(:,imu) = mu_minmax
    params_minmax(:,ibeta) = beta_minmax
    params_minmax(:,ilmin_coef) = lmin_coef_minmax
    params_minmax(:,iomicron) = omicron_minmax
    params_minmax(:,izeta_vrnce_rat) = zeta_vrnce_rat_minmax
    params_minmax(:,iupsilon_precip_frac_rat) = upsilon_precip_frac_rat_minmax
    params_minmax(:,ilambda0_stability_coef) = lambda0_stability_coef_minmax
    params_minmax(:,imult_coef) = mult_coef_minmax
    params_minmax(:,itaumin) = taumin_minmax
    params_minmax(:,itaumax) = taumax_minmax
    params_minmax(:,iLscale_mu_coef) = Lscale_mu_coef_minmax
    params_minmax(:,iLscale_pert_coef) = Lscale_pert_coef_minmax
    params_minmax(:,ialpha_corr) = alpha_corr_minmax
    params_minmax(:,iSkw_denom_coef) = Skw_denom_coef_minmax
    params_minmax(:,ic_K10) = c_K10_minmax
    params_minmax(:,ic_K10h) = c_K10h_minmax
    params_minmax(:,ithlp2_rad_coef) = thlp2_rad_coef_minmax
    params_minmax(:,ithlp2_rad_cloud_frac_thresh) = thlp2_rad_cloud_frac_thresh_minmax
    params_minmax(:,iup2_sfc_coef) = up2_sfc_coef_minmax
    params_minmax(:,iSkw_max_mag) = Skw_max_mag_minmax
    params_minmax(:,ixp3_coef_base) = xp3_coef_base_minmax
    params_minmax(:,ixp3_coef_slope) = xp3_coef_slope_minmax
    params_minmax(:,ialtitude_threshold) = altitude_threshold_minmax
    params_minmax(:,irtp2_clip_coef) = rtp2_clip_coef_minmax
    params_minmax(:,iC_invrs_tau_bkgnd) = C_invrs_tau_bkgnd_minmax
    params_minmax(:,iC_invrs_tau_sfc) = C_invrs_tau_sfc_minmax
    params_minmax(:,iC_invrs_tau_shear) = C_invrs_tau_shear_minmax
    params_minmax(:,iC_invrs_tau_N2) = C_invrs_tau_N2_minmax
    params_minmax(:,iC_invrs_tau_N2_wp2) = C_invrs_tau_N2_wp2_minmax
    params_minmax(:,iC_invrs_tau_N2_xp2) = C_invrs_tau_N2_xp2_minmax
    params_minmax(:,iC_invrs_tau_N2_wpxp) = C_invrs_tau_N2_wpxp_minmax
    params_minmax(:,iC_invrs_tau_N2_clear_wp3) = C_invrs_tau_N2_clear_wp3_minmax
    params_minmax(:,iC_invrs_tau_wpxp_Ri) = C_invrs_tau_wpxp_Ri_minmax
    params_minmax(:,iC_invrs_tau_wpxp_N2_thresh) = C_invrs_tau_wpxp_N2_thresh_minmax
    params_minmax(:,iRichardson_num_min) = Richardson_num_min_minmax
    params_minmax(:,iRichardson_num_max) = Richardson_num_max_minmax
    params_minmax(:,iCx_min) = Cx_min_minmax
    params_minmax(:,iCx_max) = Cx_max_minmax

    ! Error checks:  if a minimum value is entered, it must have a
    ! corresponding maximum value of greater value; the min and max values
    ! should not be equal or too close together (current threshold 0.01,
    ! although there may be circumstances where the user might want to adjust
    ! this); and neither min or max should be less than zero.
    do i = 1, nparams, 1
      if ( params_minmax(1,i) > params_minmax(2,i) ) then
        write(fstderr,*) "Check init_minmax namelist: " // trim(params_list(i)) // &
                   " has a minimum value greater than its maximum value."
        error stop
      elseif ( params_minmax(1,i) > 0.0_core_rknd .and. &
               abs( params_minmax(2,i) - params_minmax(1,i) ) < 0.01_core_rknd ) then
        write(fstderr,*) "Check init_minmax namelist: " // trim(params_list(i)) // &
                   " has a minimum value too close to its maximum value."
        error stop
      elseif (params_minmax(1,i) < 0.0_core_rknd .or. params_minmax(2,i) < 0.0_core_rknd ) then
        write(fstderr,*) "Check init_minmax namelist: " // trim(params_list(i)) // &
                   " has a minimum and/or maximum value less than zero."
        error stop
      end if
    end do

    ! Initialize to zero
    nindex(1:nparams) = 0
    ndim = 0

    ! Determine how many variables are being changed
    do i = 1, nparams, 1
      if ( abs(params_minmax(2,i)) > eps) then
        ndim = ndim + 1   ! Increase the total
        nindex(ndim) = i  ! Set the next array index
      endif
    enddo

    return

  end subroutine read_param_minmax

  !=============================================================================
  subroutine read_param_constraints &
           ( iunit, filename, param_constraints )

    ! Description:
    ! For the tuner.  While tuning it can be useful to keep certain parameters
    ! equal to each other.  This subroutine reads in a namelist that specifies
    ! which parameters should be kept equal to another.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants_clubb, only: fstderr ! Constant

    use parameter_indices, only: &
      iC1, &
      iC1b, &
      iC6rt, &
      iC6rtb, &
      iC6rtc, &
      iC6thl, &
      iC6thlb, &
      iC6thlc, &
      iC7, &
      iC7b, &
      iC11, &
      iC11b, &
      iC14, &
      iC6rt_Lscale0, &
      iC6thl_Lscale0, &
      igamma_coef, &
      igamma_coefb

    implicit none

    ! Input variables
    integer, intent(in) :: iunit

    character(len=*), intent(in) :: filename

    ! Output variables
    character(len=28), dimension(nparams), intent(out) ::  &
      param_constraints  ! Which variables should be kept equal to others

    ! Local variables
    character(len=28) ::  &
      C1_equals, C1b_equals, C6rt_equals, C6rtb_equals, C6rtc_equals, &
      C6thl_equals, C6thlb_equals, C6thlc_equals, C7_equals, C7b_equals, &
      C11_equals, C11b_equals, C14_equals, C6rt_Lscale0_equals, &
      C6thl_Lscale0_equals, gamma_coef_equals, gamma_coefb_equals

    integer :: i

    ! This namelist specifies if a variable should be kept equal to another.
    namelist /parameter_constraints/  &
      C1_equals, C1b_equals, C6rt_equals, C6rtb_equals, C6rtc_equals, &
      C6thl_equals, C6thlb_equals, C6thlc_equals, C7_equals, C7b_equals, &
      C11_equals, C11b_equals, C14_equals, C6rt_Lscale0_equals, &
      C6thl_Lscale0_equals, gamma_coef_equals, gamma_coefb_equals

    ! Initialize output array
    param_constraints = ""

    ! Read the namelist
    open(unit=iunit, file=filename, status='old', action='read')

    read(unit=iunit, nml=parameter_constraints)

    close(unit=iunit)

    ! Put the variables in the output array
    param_constraints(iC1)            = C1_equals
    param_constraints(iC1b)           = C1b_equals
    param_constraints(iC6rt)          = C6rt_equals
    param_constraints(iC6rtb)         = C6rtb_equals
    param_constraints(iC6rtc)         = C6rtc_equals
    param_constraints(iC6thl)         = C6thl_equals
    param_constraints(iC6thlb)        = C6thlb_equals
    param_constraints(iC6thlc)        = C6thlc_equals
    param_constraints(iC7)            = C7_equals
    param_constraints(iC7b)           = C7b_equals
    param_constraints(iC11)           = C11_equals
    param_constraints(iC11b)          = C11b_equals
    param_constraints(iC14)           = C14_equals
    param_constraints(iC6rt_Lscale0)  = C6rt_Lscale0_equals
    param_constraints(iC6thl_Lscale0) = C6thl_Lscale0_equals
    param_constraints(igamma_coef)    = gamma_coef_equals
    param_constraints(igamma_coefb)   = gamma_coefb_equals

    ! Error check to make sure no constraint parameter is set equal to itself.
    do i = 1, nparams, 1
      if ( params_list(i) == param_constraints(i) ) then
        write(fstderr,*) "Check parameter_constraints namelist: "//trim( params_list(i) )// &
          " was set equal to itself."
        error stop
      end if
    end do

    return

  end subroutine read_param_constraints

  !=============================================================================
  subroutine pack_parameters &
             ( C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
               C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
               C7, C7b, C7c, C8, C8b, C10, &
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
               Cx_min, Cx_max, Richardson_num_min, Richardson_num_max, params )

    ! Description:
    ! Takes the list of scalar variables and puts them into a 1D vector.
    ! It is here for the purpose of keeping the code generalized
    ! when new variables are added.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameter_indices, only: & 
        iC1,  & ! Variable(s)
        iC1b, & 
        iC1c, & 
        iC2rt, & 
        iC2thl, & 
        iC2rtthl, & 
        iC4, & 
        iC_uu_shr, &
        iC_uu_buoy, & 
        iC6rt, & 
        iC6rtb, & 
        iC6rtc, & 
        iC6thl, & 
        iC6thlb, & 
        iC6thlc, & 
        iC7, & 
        iC7b, & 
        iC7c, & 
        iC8, & 
        iC8b, & 
        iC10, & 
        iC11, & 
        iC11b, & 
        iC11c, & 
        iC12, & 
        iC13, & 
        iC14, &
        iC_wp2_pr_dfsn, &
        iC_wp3_pr_tp, &
        iC_wp3_pr_turb, &
        iC_wp3_pr_dfsn, &
        iC_wp2_splat

    use parameter_indices, only: &
        iC6rt_Lscale0, &
        iC6thl_Lscale0, &
        iC7_Lscale0, &
        iwpxp_L_thresh

    use parameter_indices, only: & 
        ic_K,  & 
        ic_K1, & 
        inu1, & 
        ic_K2, & 
        inu2, & 
        ic_K6, & 
        inu6, & 
        ic_K8, & 
        inu8, & 
        ic_K9, & 
        inu9, & 
        inu10, &
        ic_K_hm, & 
        ic_K_hmb, & 
        iK_hm_min_coef, &
        inu_hm, & 
        islope_coef_spread_DG_means_w, &
        ipdf_component_stdev_factor_w, &
        icoef_spread_DG_means_rt, &
        icoef_spread_DG_means_thl, &
        igamma_coef, & 
        igamma_coefb, & 
        igamma_coefc, & 
        imu, & 
        ibeta, & 
        ilmin_coef, &
        iomicron, &
        izeta_vrnce_rat, &
        iupsilon_precip_frac_rat, &
        ilambda0_stability_coef, &
        imult_coef, &
        itaumin, & 
        itaumax, & 
        iLscale_mu_coef, &
        iLscale_pert_coef, &
        ialpha_corr, &
        iSkw_denom_coef, &
        ic_K10, &
        ic_K10h, &
        ithlp2_rad_coef, &
        ithlp2_rad_cloud_frac_thresh, &
        iup2_sfc_coef, &
        iSkw_max_mag, &
        ixp3_coef_base, &
        ixp3_coef_slope, &
        ialtitude_threshold, &
        irtp2_clip_coef, &
        iC_invrs_tau_bkgnd, &
        iC_invrs_tau_sfc, &
        iC_invrs_tau_shear, &
        iC_invrs_tau_N2, &
        iC_invrs_tau_N2_wp2, &
        iC_invrs_tau_N2_xp2, &
        iC_invrs_tau_N2_wpxp, &
        iC_invrs_tau_N2_clear_wp3, &
        iC_invrs_tau_wpxp_Ri, &
        iC_invrs_tau_wpxp_N2_thresh, &
        iCx_min, &
        iCx_max, &
        iRichardson_num_min, &
        iRichardson_num_max

    implicit none

    ! Input variables
    real( kind = core_rknd ), intent(in) :: & 
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
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max

    ! Output variables
    real( kind = core_rknd ), intent(out), dimension(nparams) :: params

    params(iC1)      = C1
    params(iC1b)     = C1b
    params(iC1c)     = C1c
    params(iC2rt)    = C2rt
    params(iC2thl)   = C2thl
    params(iC2rtthl) = C2rtthl
    params(iC4)      = C4
    params(iC_uu_shr) = C_uu_shr
    params(iC_uu_buoy) = C_uu_buoy
    params(iC6rt)    = C6rt
    params(iC6rtb)   = C6rtb
    params(iC6rtc)   = C6rtc
    params(iC6thl)   = C6thl
    params(iC6thlb)  = C6thlb
    params(iC6thlc)  = C6thlc
    params(iC7)      = C7
    params(iC7b)     = C7b
    params(iC7c)     = C7c
    params(iC8)      = C8
    params(iC8b)     = C8b
    params(iC10)     = C10
    params(iC11)     = C11
    params(iC11b)    = C11b
    params(iC11c)    = C11c
    params(iC12)     = C12
    params(iC13)     = C13
    params(iC14)     = C14
    params(iC_wp2_pr_dfsn)      = C_wp2_pr_dfsn
    params(iC_wp3_pr_tp)        = C_wp3_pr_tp
    params(iC_wp3_pr_turb)      = C_wp3_pr_turb
    params(iC_wp3_pr_dfsn)      = C_wp3_pr_dfsn
    params(iC_wp2_splat)        = C_wp2_splat
    params(iC6rt_Lscale0)       = C6rt_Lscale0
    params(iC6thl_Lscale0)      = C6thl_Lscale0
    params(iC7_Lscale0)         = C7_Lscale0
    params(iwpxp_L_thresh)    = wpxp_L_thresh
    params(ic_K)       = c_K
    params(ic_K1)      = c_K1
    params(inu1)       = nu1
    params(ic_K2)      = c_K2
    params(inu2)       = nu2
    params(ic_K6)      = c_K6
    params(inu6)       = nu6
    params(ic_K8)      = c_K8
    params(inu8)       = nu8
    params(ic_K9)      = c_K9
    params(inu9)       = nu9
    params(inu10)      = nu10
    params(ic_K_hm)    = c_K_hm
    params(ic_K_hmb)   = c_K_hmb
    params(iK_hm_min_coef)   = K_hm_min_coef
    params(inu_hm)     = nu_hm
    params(islope_coef_spread_DG_means_w) = slope_coef_spread_DG_means_w
    params(ipdf_component_stdev_factor_w) = pdf_component_stdev_factor_w
    params(icoef_spread_DG_means_rt) = coef_spread_DG_means_rt
    params(icoef_spread_DG_means_thl) = coef_spread_DG_means_thl
    params(igamma_coef)  = gamma_coef
    params(igamma_coefb) = gamma_coefb
    params(igamma_coefc) = gamma_coefc
    params(imu) = mu
    params(ibeta) = beta
    params(ilmin_coef) = lmin_coef
    params(iomicron) = omicron
    params(izeta_vrnce_rat) = zeta_vrnce_rat
    params(iupsilon_precip_frac_rat) = upsilon_precip_frac_rat
    params(ilambda0_stability_coef) = lambda0_stability_coef
    params(imult_coef) = mult_coef
    params(itaumin) = taumin
    params(itaumax) = taumax
    params(iLscale_mu_coef) = Lscale_mu_coef
    params(iLscale_pert_coef) = Lscale_pert_coef
    params(ialpha_corr) = alpha_corr
    params(iSkw_denom_coef) = Skw_denom_coef
    params(ic_K10) = c_K10
    params(ic_K10h) = c_K10h
    params(ithlp2_rad_coef) = thlp2_rad_coef
    params(ithlp2_rad_cloud_frac_thresh) = thlp2_rad_cloud_frac_thresh
    params(iup2_sfc_coef) = up2_sfc_coef
    params(iSkw_max_mag) = Skw_max_mag
    params(ixp3_coef_base) = xp3_coef_base
    params(ixp3_coef_slope) = xp3_coef_slope
    params(ialtitude_threshold) = altitude_threshold
    params(irtp2_clip_coef) = rtp2_clip_coef
    params(iC_invrs_tau_bkgnd)          = C_invrs_tau_bkgnd
    params(iC_invrs_tau_sfc)            = C_invrs_tau_sfc
    params(iC_invrs_tau_shear)          = C_invrs_tau_shear
    params(iC_invrs_tau_N2)             = C_invrs_tau_N2
    params(iC_invrs_tau_N2_wp2)         = C_invrs_tau_N2_wp2
    params(iC_invrs_tau_N2_xp2)         = C_invrs_tau_N2_xp2
    params(iC_invrs_tau_N2_wpxp)        = C_invrs_tau_N2_wpxp
    params(iC_invrs_tau_N2_clear_wp3)   = C_invrs_tau_N2_clear_wp3
    params(iC_invrs_tau_wpxp_Ri)        = C_invrs_tau_wpxp_Ri
    params(iC_invrs_tau_wpxp_N2_thresh) = C_invrs_tau_wpxp_N2_thresh
    params(iCx_min) = Cx_min
    params(iCx_max) = Cx_max
    params(iRichardson_num_min) = Richardson_num_min
    params(iRichardson_num_max) = Richardson_num_max


    return
  end subroutine pack_parameters

  !=============================================================================
  subroutine unpack_parameters & 
             ( params, & 
               C1, C1b, C1c, C2rt, C2thl, C2rtthl, &
               C4, C_uu_shr, C_uu_buoy, C6rt, C6rtb, C6rtc, C6thl, C6thlb, C6thlc, &
               C7, C7b, C7c, C8, C8b, C10, &
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
               Cx_min, Cx_max, Richardson_num_min, Richardson_num_max )

    ! Description:
    ! Takes the 1D vector and returns the list of scalar variables.
    ! Here for the purposes of keeping the code generalized
    ! when new variables are added.

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use parameter_indices, only: & 
        iC1,  & ! Variable(s)
        iC1b, & 
        iC1c, & 
        iC2rt, & 
        iC2thl, & 
        iC2rtthl, & 
        iC4, & 
        iC_uu_shr, &
        iC_uu_buoy, & 
        iC6rt, & 
        iC6rtb, & 
        iC6rtc, & 
        iC6thl, & 
        iC6thlb, & 
        iC6thlc, & 
        iC7, & 
        iC7b, & 
        iC7c, & 
        iC8, & 
        iC8b, & 
        iC10, & 
        iC11, & 
        iC11b, & 
        iC11c, & 
        iC12, & 
        iC13, & 
        iC14, &
        iC_wp2_pr_dfsn, &
        iC_wp3_pr_tp, &
        iC_wp3_pr_turb, &
        iC_wp3_pr_dfsn, &
        iC_wp2_splat

    use parameter_indices, only: &
        iC6rt_Lscale0, &
        iC6thl_Lscale0, &
        iC7_Lscale0, &
        iwpxp_L_thresh

    use parameter_indices, only: & 
        ic_K,  & 
        ic_K1, & 
        inu1, & 
        ic_K2, & 
        inu2, & 
        ic_K6, & 
        inu6, & 
        ic_K8, & 
        inu8, & 
        ic_K9, & 
        inu9, & 
        inu10, &
        ic_K_hm, & 
        ic_K_hmb, & 
        iK_hm_min_coef, & 
        inu_hm, & 
        islope_coef_spread_DG_means_w, &
        ipdf_component_stdev_factor_w, &
        icoef_spread_DG_means_rt, &
        icoef_spread_DG_means_thl, &
        igamma_coef, & 
        igamma_coefb, & 
        igamma_coefc, & 
        imu, & 
        ibeta, & 
        ilmin_coef, &
        iomicron, &
        izeta_vrnce_rat, &
        iupsilon_precip_frac_rat, &
        ilambda0_stability_coef, &
        imult_coef, &
        itaumin, & 
        itaumax, & 
        iLscale_mu_coef, &
        iLscale_pert_coef, &
        ialpha_corr, &
        iSkw_denom_coef, &
        ic_K10, &
        ic_K10h, & 
        ithlp2_rad_coef, &
        ithlp2_rad_cloud_frac_thresh, &
        iup2_sfc_coef, &
        iSkw_max_mag, &
        ixp3_coef_base, &
        ixp3_coef_slope, &
        ialtitude_threshold, &
        irtp2_clip_coef, &
        iC_invrs_tau_bkgnd, &
        iC_invrs_tau_sfc, &
        iC_invrs_tau_shear, &
        iC_invrs_tau_N2, &
        iC_invrs_tau_N2_wp2, &
        iC_invrs_tau_N2_xp2, &
        iC_invrs_tau_N2_wpxp, &
        iC_invrs_tau_N2_clear_wp3, &
        iC_invrs_tau_wpxp_Ri, &
        iC_invrs_tau_wpxp_N2_thresh, &
        iCx_min, &
        iCx_max, &
        iRichardson_num_min, &
        iRichardson_num_max, &
        nparams

    implicit none

    ! Input variables
    real( kind = core_rknd ), intent(in), dimension(nparams) :: params

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
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max

    C1      = params(iC1)
    C1b     = params(iC1b)
    C1c     = params(iC1c)
    C2rt    = params(iC2rt)
    C2thl   = params(iC2thl)
    C2rtthl = params(iC2rtthl)
    C4      = params(iC4)
    C_uu_shr = params(iC_uu_shr)
    C_uu_buoy = params(iC_uu_buoy)
    C6rt    = params(iC6rt)
    C6rtb   = params(iC6rtb)
    C6rtc   = params(iC6rtc)
    C6thl   = params(iC6thl)
    C6thlb  = params(iC6thlb)
    C6thlc  = params(iC6thlc)
    C7      = params(iC7)
    C7b     = params(iC7b)
    C7c     = params(iC7c)
    C8      = params(iC8)
    C8b     = params(iC8b)
    C10     = params(iC10)
    C11     = params(iC11)
    C11b    = params(iC11b)
    C11c    = params(iC11c)
    C12     = params(iC12)
    C13     = params(iC13)
    C14     = params(iC14)

    C_wp2_pr_dfsn      = params(iC_wp2_pr_dfsn)
    C_wp3_pr_tp        = params(iC_wp3_pr_tp)
    C_wp3_pr_turb      = params(iC_wp3_pr_turb)
    C_wp3_pr_dfsn      = params(iC_wp3_pr_dfsn)
    C_wp2_splat        = params(iC_wp2_splat)

    C6rt_Lscale0       = params(iC6rt_Lscale0)
    C6thl_Lscale0      = params(iC6thl_Lscale0)
    C7_Lscale0         = params(iC7_Lscale0)
    wpxp_L_thresh      = params(iwpxp_L_thresh)

    c_K       = params(ic_K)
    c_K1      = params(ic_K1)
    nu1       = params(inu1)
    c_K2      = params(ic_K2)
    nu2       = params(inu2)
    c_K6      = params(ic_K6)
    nu6       = params(inu6)
    c_K8      = params(ic_K8)
    nu8       = params(inu8)
    c_K9      = params(ic_K9)
    nu9       = params(inu9)
    nu10      = params(inu10)
    c_K_hm    = params(ic_K_hm)
    c_K_hmb   = params(ic_K_hmb)
    K_hm_min_coef   = params(iK_hm_min_coef)
    nu_hm     = params(inu_hm)

    slope_coef_spread_DG_means_w = params(islope_coef_spread_DG_means_w)
    pdf_component_stdev_factor_w = params(ipdf_component_stdev_factor_w)
    coef_spread_DG_means_rt = params(icoef_spread_DG_means_rt)
    coef_spread_DG_means_thl = params(icoef_spread_DG_means_thl)

    gamma_coef  = params(igamma_coef)
    gamma_coefb = params(igamma_coefb)
    gamma_coefc = params(igamma_coefc)

    mu = params(imu)

    beta = params(ibeta)

    lmin_coef = params(ilmin_coef)

    omicron = params(iomicron)
    zeta_vrnce_rat = params(izeta_vrnce_rat)

    upsilon_precip_frac_rat = params(iupsilon_precip_frac_rat)
    lambda0_stability_coef = params(ilambda0_stability_coef)
    mult_coef = params(imult_coef)

    taumin = params(itaumin)
    taumax = params(itaumax)

    Lscale_mu_coef = params(iLscale_mu_coef)
    Lscale_pert_coef = params(iLscale_pert_coef)
    alpha_corr = params(ialpha_corr)
    Skw_denom_coef = params(iSkw_denom_coef)
    c_K10 = params(ic_K10)
    c_K10h = params(ic_K10h)

    thlp2_rad_coef = params(ithlp2_rad_coef)
    thlp2_rad_cloud_frac_thresh = params(ithlp2_rad_cloud_frac_thresh)
    up2_sfc_coef = params(iup2_sfc_coef)
    Skw_max_mag = params(iSkw_max_mag)
    xp3_coef_base = params(ixp3_coef_base)
    xp3_coef_slope = params(ixp3_coef_slope)
    altitude_threshold = params(ialtitude_threshold)
    rtp2_clip_coef = params(irtp2_clip_coef)
    C_invrs_tau_bkgnd          = params(iC_invrs_tau_bkgnd)
    C_invrs_tau_sfc            = params(iC_invrs_tau_sfc )
    C_invrs_tau_shear          = params(iC_invrs_tau_shear)
    C_invrs_tau_N2             = params(iC_invrs_tau_N2)
    C_invrs_tau_N2_wp2         = params(iC_invrs_tau_N2_wp2)
    C_invrs_tau_N2_xp2         = params(iC_invrs_tau_N2_xp2)
    C_invrs_tau_N2_wpxp        = params(iC_invrs_tau_N2_wpxp)
    C_invrs_tau_N2_clear_wp3   = params(iC_invrs_tau_N2_clear_wp3)
    C_invrs_tau_wpxp_Ri        = params(iC_invrs_tau_wpxp_Ri)
    C_invrs_tau_wpxp_N2_thresh = params(iC_invrs_tau_wpxp_N2_thresh)
    Cx_min = params(iCx_min)
    Cx_max = params(iCx_max)
    Richardson_num_min = params(iRichardson_num_min)
    Richardson_num_max = params(iRichardson_num_max)

    return
  end subroutine unpack_parameters

  !=============================================================================
  subroutine init_parameters_999( &
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
               Cx_min, Cx_max, Richardson_num_min, Richardson_num_max )

    ! Description:
    ! Set all tunable parameters to NaN

    ! References:
    ! None
    !-----------------------------------------------------------------------

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
      Cx_min, Cx_max, Richardson_num_min, Richardson_num_max

    ! --- Begin Code ---

    C1                           = init_value
    C1b                          = init_value
    C1c                          = init_value
    C2rt                         = init_value
    C2thl                        = init_value
    C2rtthl                      = init_value
    C4                           = init_value
    C_uu_shr                     = init_value
    C_uu_buoy                    = init_value
    C6rt                         = init_value
    C6rtb                        = init_value
    C6rtc                        = init_value
    C6thl                        = init_value
    C6thlb                       = init_value
    C6thlc                       = init_value
    C7                           = init_value
    C7b                          = init_value
    C7c                          = init_value
    C8                           = init_value
    C8b                          = init_value
    C10                          = init_value
    C11                          = init_value
    C11b                         = init_value
    C11c                         = init_value
    C12                          = init_value
    C13                          = init_value
    C14                          = init_value
    C_wp2_pr_dfsn                = init_value
    C_wp3_pr_tp                  = init_value
    C_wp3_pr_turb                = init_value
    C_wp3_pr_dfsn                = init_value
    C_wp2_splat                  = init_value 
    C6rt_Lscale0                 = init_value
    C6thl_Lscale0                = init_value
    C7_Lscale0                   = init_value
    wpxp_L_thresh                = init_value
    c_K                          = init_value
    c_K1                         = init_value
    nu1                          = init_value
    c_K2                         = init_value
    nu2                          = init_value
    c_K6                         = init_value
    nu6                          = init_value
    c_K8                         = init_value
    nu8                          = init_value
    c_K9                         = init_value
    nu9                          = init_value
    nu10                         = init_value
    c_K_hm                       = init_value
    c_K_hmb                      = init_value
    K_hm_min_coef                = init_value
    nu_hm                        = init_value
    slope_coef_spread_DG_means_w = init_value
    pdf_component_stdev_factor_w = init_value
    coef_spread_DG_means_rt      = init_value
    coef_spread_DG_means_thl     = init_value
    beta                         = init_value
    gamma_coef                   = init_value
    gamma_coefb                  = init_value
    gamma_coefc                  = init_value
    mult_coef                    = init_value
    taumin                       = init_value
    taumax                       = init_value
    lmin_coef                    = init_value
    omicron                      = init_value
    zeta_vrnce_rat               = init_value
    upsilon_precip_frac_rat      = init_value
    lambda0_stability_coef       = init_value
    mu                           = init_value
    Lscale_mu_coef               = init_value
    Lscale_pert_coef             = init_value
    alpha_corr                   = init_value
    Skw_denom_coef               = init_value
    c_K10                        = init_value
    c_K10h                       = init_value
    thlp2_rad_coef               = init_value
    thlp2_rad_cloud_frac_thresh  = init_value
    up2_sfc_coef                 = init_value
    Skw_max_mag                  = init_value
    xp3_coef_base                = init_value
    xp3_coef_slope               = init_value
    altitude_threshold           = init_value
    rtp2_clip_coef               = init_value
    C_invrs_tau_bkgnd            = init_value 
    C_invrs_tau_sfc              = init_value 
    C_invrs_tau_shear            = init_value 
    C_invrs_tau_N2               = init_value 
    C_invrs_tau_N2_xp2           = init_value
    C_invrs_tau_N2_wp2           = init_value
    C_invrs_tau_N2_wpxp          = init_value
    C_invrs_tau_N2_clear_wp3     = init_value
    C_invrs_tau_wpxp_Ri          = init_value
    C_invrs_tau_wpxp_N2_thresh   = init_value
    Cx_min                       = init_value
    Cx_max                       = init_value
    Richardson_num_min           = init_value
    Richardson_num_max           = init_value

    return

  end subroutine init_parameters_999

!===============================================================================

end module parameters_tunable
