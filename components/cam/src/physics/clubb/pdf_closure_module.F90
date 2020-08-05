!---------------------------------------------------------------------------
! $Id$
!===============================================================================
module pdf_closure_module

  implicit none

  public :: pdf_closure, &
            calc_wp4_pdf, &
            calc_wp2xp_pdf, &
            calc_wpxp2_pdf, &
            calc_wpxpyp_pdf, &
            calc_vert_avg_cf_component

  ! Options for the two component normal (double Gaussian) PDF type to use for
  ! the w, rt, and theta-l (or w, chi, and eta) portion of CLUBB's multivariate,
  ! two-component PDF.
  integer, parameter, public :: &
    iiPDF_ADG1 = 1,     & ! ADG1 PDF
    iiPDF_ADG2 = 2,     & ! ADG2 PDF
    iiPDF_3D_Luhar = 3, & ! 3D Luhar PDF
    iiPDF_new = 4,      & ! new PDF
    iiPDF_TSDADG = 5,   & ! new TSDADG PDF
    iiPDF_LY93 = 6        ! Lewellen and Yoh (1993)

  ! The selected two component normal PDF for w, rt, and theta-l.
  integer, parameter, public :: &
    iiPDF_type = iiPDF_ADG1

  private ! Set Default Scope

  contains
!------------------------------------------------------------------------

  !#######################################################################
  !#######################################################################
  ! If you change the argument list of pdf_closure you also have to
  ! change the calls to this function in the host models CAM, WRF, SAM
  ! and GFDL.
  !#######################################################################
  !#######################################################################
  subroutine pdf_closure( hydromet_dim, p_in_Pa, exner, thv_ds,     &
                          wm, wp2, wp3, sigma_sqd_w,                &
                          Skw, Skthl_in, Skrt_in,                   &
                          rtm, rtp2, wprtp,                         &
                          thlm, thlp2, wpthlp,                      &
                          um, up2, upwp,                            &
                          vm, vp2, vpwp,                            &
                          rtpthlp,                                  &
                          sclrm, wpsclrp, sclrp2,                   &
                          sclrprtp, sclrpthlp,                      &
#ifdef GFDL
                          RH_crit, do_liquid_only_in_clubb,         & ! h1g, 2010-06-15
#endif
                          wphydrometp, wp2hmp,                      &
                          rtphmp, thlphmp,                          &
                          wp4, wprtp2, wp2rtp,                      &
                          wpthlp2, wp2thlp, wprtpthlp,              &
                          cloud_frac, ice_supersat_frac,            &
                          rcm, wpthvp, wp2thvp, rtpthvp,            &
                          thlpthvp, wprcp, wp2rcp, rtprcp,          &
                          thlprcp, rcp2,                            &
                          uprcp, vprcp,                             &
                          pdf_params, pdf_implicit_coefs_terms,     &
                          F_w, F_rt, F_thl,                         &
                          min_F_w, max_F_w,                         &
                          min_F_rt, max_F_rt,                       &
                          min_F_thl, max_F_thl,                     &
                          wpsclrprtp, wpsclrp2, sclrpthvp,          &
                          wpsclrpthlp, sclrprcp, wp2sclrp,          &
                          rc_coef                                   )


    ! Description:
    ! Subroutine that computes pdf parameters analytically.
    !
    ! Based of the original formulation, but with some tweaks
    ! to remove some of the less realistic assumptions and
    ! improve transport terms.

    !   Corrected version that should remove inconsistency

    ! References:
    !   The shape of CLUBB's PDF is given by the expression in
    !   https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:clubb_pdf

    !   Eqn. 29, 30, 31, 32 & 33  on p. 3547 of
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &  ! Constants
        three,          & ! 3
        two,            & ! 2
        one,            & ! 1
        one_half,       & ! 1/2
        zero,           & ! 0
        Cp,             & ! Dry air specific heat at constant p [J/kg/K]
        Lv,             & ! Latent heat of vaporization         [J/kg]
        Rd,             & ! Dry air gas constant                [J/kg/K]
        ep,             & ! Rd / Rv;     ep  = 0.622            [-]
        ep1,            & ! (1.0-ep)/ep; ep1 = 0.61             [-]
        ep2,            & ! 1.0/ep;      ep2 = 1.61             [-]
        rt_tol,         & ! Tolerance for r_t                   [kg/kg]
        thl_tol,        & ! Tolerance for th_l                  [K]
        T_freeze_K,     & ! Freezing point of water             [K]
        fstderr,        &
        zero_threshold, &
        chi_tol, &
        eta_tol, &
        max_mag_correlation, &
        eps, &
        w_tol

    use parameters_model, only: &
        mixt_frac_max_mag, & ! Variable(s)
        sclr_dim             ! Number of passive scalar variables

    use parameters_tunable, only: & 
        Skw_denom_coef ! Variable(s)

    use pdf_parameter_module, only:  &
        pdf_parameter,        & ! Variable Type
        implicit_coefs_terms

    use new_pdf_main, only: &
        new_pdf_driver    ! Procedure(s)

    use adg1_adg2_3d_luhar_pdf, only: &
        ADG1_pdf_driver,     & ! Procedure(s)
        ADG2_pdf_driver,     &
        Luhar_3D_pdf_driver

    use new_tsdadg_pdf, only: &
        tsdadg_pdf_driver    ! Procedure(s)

    use LY93_pdf, only: &
        LY93_driver    ! Procedure(s)

    use pdf_utilities, only: &
        calc_comp_corrs_binormal, & ! Procedure(s)
        calc_corr_chi_x,          &
        calc_corr_eta_x

    use array_index, only: &
        l_mix_rat_hm  ! Variable(s)

    use model_flags, only: &
        l_explicit_turbulent_adv_wp3,  & ! Variable(s)
        l_explicit_turbulent_adv_xpyp

    use numerical_check, only:  & 
        pdf_closure_check ! Procedure(s)

    use saturation, only:  & 
        sat_mixrat_liq, & ! Procedure(s)
        sat_mixrat_ice

    use stats_variables, only: &
        iwp4,       & ! Variables
        ircp2,      &
        iwprtp2,    &
        iwprtpthlp, &
        iwpthlp2

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use error_code, only: &
        clubb_at_least_debug_level,  & ! Procedure
        err_code,                    & ! Error Indicator
        clubb_fatal_error              ! Constant

    implicit none

    intrinsic :: sqrt, exp, min, max, abs, present

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim   ! Number of hydrometeor species              [#]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  & 
      p_in_Pa,     & ! Pressure                                   [Pa]
      exner,       & ! Exner function                             [-]
      thv_ds,      & ! Dry, base-state theta_v (ref. th_l here)   [K]
      wm,          & ! mean w-wind component (vertical velocity)  [m/s] 
      wp2,         & ! w'^2                                       [m^2/s^2] 
      wp3,         & ! w'^3                                       [m^3/s^3]
      Skw,         & ! Skewness of w                              [-]
      Skthl_in,    & ! Skewness of thl                            [-]
      Skrt_in,     & ! Skewness of rt                             [-]
      rtm,         & ! Mean total water mixing ratio              [kg/kg]
      rtp2,        & ! r_t'^2                                     [(kg/kg)^2]
      wprtp,       & ! w'r_t'                                     [(kg/kg)(m/s)]
      thlm,        & ! Mean liquid water potential temperature    [K]
      thlp2,       & ! th_l'^2                                    [K^2]
      wpthlp,      & ! w'th_l'                                    [K(m/s)]
      rtpthlp        ! r_t'th_l'                                  [K(kg/kg)]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      um,          & ! Grid-mean eastward wind     [m/s]
      up2,         & ! u'^2                        [(m/s)^2]
      upwp,        & ! u'w'                        [(m/s)^2]
      vm,          & ! Grid-mean northward wind    [m/s]
      vp2,         & ! v'^2                        [(m/s)^2]
      vpwp           ! v'w'                        [(m/s)^2]

    real( kind = core_rknd ), dimension(gr%nz, sclr_dim), intent(in) ::  & 
      sclrm,       & ! Mean passive scalar        [units vary]
      wpsclrp,     & ! w' sclr'                   [units vary]
      sclrp2,      & ! sclr'^2                    [units vary]
      sclrprtp,    & ! sclr' r_t'                 [units vary]
      sclrpthlp      ! sclr' th_l'                [units vary]

#ifdef  GFDL
    ! critial relative humidity for nucleation
    real( kind = core_rknd ), dimension( min(1,sclr_dim), 2 ), intent(in) ::  & ! h1g, 2010-06-15
       RH_crit     ! critical relative humidity for droplet and ice nucleation
! ---> h1g, 2012-06-14
    logical, intent(in)                 ::  do_liquid_only_in_clubb
! <--- h1g, 2012-06-14
#endif

    real( kind = core_rknd ), dimension(gr%nz, hydromet_dim), intent(in) :: &
      wphydrometp, & ! Covariance of w and a hydrometeor    [(m/s) <hm units>]
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' >    [(m/s)^2 <hm units>]
      rtphmp,      & ! Covariance of rt and a hydrometeor   [(kg/kg) <hm units>]
      thlphmp        ! Covariance of thl and a hydrometeor  [K <hm units>]

    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      ! If iiPDF_type == iiPDF_ADG2, this gets overwritten. Therefore,
      ! intent(inout). Otherwise it should be intent(in)
      sigma_sqd_w   ! Width of individual w plumes               [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  & 
      wp4,                & ! w'^4                  [m^4/s^4]
      wprtp2,             & ! w' r_t'               [(m kg)/(s kg)]
      wp2rtp,             & ! w'^2 r_t'             [(m^2 kg)/(s^2 kg)]
      wpthlp2,            & ! w' th_l'^2            [(m K^2)/s]
      wp2thlp,            & ! w'^2 th_l'            [(m^2 K)/s^2]
      cloud_frac,         & ! Cloud fraction        [-]
      ice_supersat_frac,  & ! Ice cloud fracion     [-]
      rcm,                & ! Mean liquid water     [kg/kg]
      wpthvp,             & ! Buoyancy flux         [(K m)/s] 
      wp2thvp,            & ! w'^2 th_v'            [(m^2 K)/s^2]
      rtpthvp,            & ! r_t' th_v'            [(kg K)/kg]
      thlpthvp,           & ! th_l' th_v'           [K^2]
      wprcp,              & ! w' r_c'               [(m kg)/(s kg)]
      wp2rcp,             & ! w'^2 r_c'             [(m^2 kg)/(s^2 kg)]
      rtprcp,             & ! r_t' r_c'             [(kg^2)/(kg^2)]
      thlprcp,            & ! th_l' r_c'            [(K kg)/kg]
      rcp2,               & ! r_c'^2                [(kg^2)/(kg^2)]
      wprtpthlp             ! w' r_t' th_l'         [(m kg K)/(s kg)]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      uprcp,              & ! u' r_c'               [(m kg)/(s kg)]
      vprcp                 ! v' r_c'               [(m kg)/(s kg)]

    type(pdf_parameter), intent(inout) :: & 
      pdf_params     ! pdf paramters         [units vary]

    type(implicit_coefs_terms), dimension(gr%nz), intent(out) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Parameters output only for recording statistics (new PDF).
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      F_w,   & ! Parameter for the spread of the PDF component means of w    [-]
      F_rt,  & ! Parameter for the spread of the PDF component means of rt   [-]
      F_thl    ! Parameter for the spread of the PDF component means of thl  [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      min_F_w,   & ! Minimum allowable value of parameter F_w      [-]
      max_F_w,   & ! Maximum allowable value of parameter F_w      [-]
      min_F_rt,  & ! Minimum allowable value of parameter F_rt     [-]
      max_F_rt,  & ! Maximum allowable value of parameter F_rt     [-]
      min_F_thl, & ! Minimum allowable value of parameter F_thl    [-]
      max_F_thl    ! Maximum allowable value of parameter F_thl    [-]

    ! Output (passive scalar variables)

    real( kind = core_rknd ), intent(out), dimension(gr%nz, sclr_dim) ::  & 
      sclrpthvp, & 
      sclrprcp, & 
      wpsclrp2, & 
      wpsclrprtp, & 
      wpsclrpthlp, & 
      wp2sclrp

    ! Local Variables

    ! Variables that are stored in derived data type pdf_params.
    real( kind = core_rknd ), dimension(gr%nz) ::  &
      w_1,           & ! Mean of w (1st PDF component)                     [m/s]
      w_2,           & ! Mean of w (2nd PDF component)                     [m/s]
      varnce_w_1,    & ! Variance of w (1st PDF component)             [m^2/s^2]
      varnce_w_2,    & ! Variance of w (2nd PDF component)             [m^2/s^2]
      rt_1,          & ! Mean of r_t (1st PDF component)                 [kg/kg]
      rt_2,          & ! Mean of r_t (2nd PDF component)                 [kg/kg]
      varnce_rt_1,   & ! Variance of r_t (1st PDF component)         [kg^2/kg^2]
      varnce_rt_2,   & ! Variance of r_t (2nd PDF component)         [kg^2/kg^2]
      thl_1,         & ! Mean of th_l (1st PDF component)                    [K]
      thl_2,         & ! Mean of th_l (2nd PDF component)                    [K]
      varnce_thl_1,  & ! Variance of th_l (1st PDF component)              [K^2]
      varnce_thl_2,  & ! Variance of th_l (2nd PDF component)              [K^2]
      u_1,           & ! Mean of eastward wind (1st PDF component)         [m/s]
      u_2,           & ! Mean of eastward wind (2nd PDF component)         [m/s]
      varnce_u_1,    & ! Variance of u (1st PDF component)             [m^2/s^2]
      varnce_u_2,    & ! Variance of u (2nd PDF component)             [m^2/s^2]
      v_1,           & ! Mean of northward wind (1st PDF component)        [m/s]
      v_2,           & ! Mean of northward wind (2nd PDF component)        [m/s]
      varnce_v_1,    & ! Variance of v (1st PDF component)             [m^2/s^2]
      varnce_v_2,    & ! Variance of v (2nd PDF component)             [m^2/s^2]
      corr_w_rt_1,   & ! Correlation of w and r_t (1st PDF component)        [-]
      corr_w_rt_2,   & ! Correlation of w and r_t (2nd PDF component)        [-]
      corr_w_thl_1,  & ! Correlation of w and th_l (1st PDF component)       [-]
      corr_w_thl_2,  & ! Correlation of w and th_l (2nd PDF component)       [-]
      corr_rt_thl_1, & ! Correlation of r_t and th_l (1st PDF component)     [-]
      corr_rt_thl_2, & ! Correlation of r_t and th_l (2nd PDF component)     [-]
      alpha_thl,     & ! Factor relating to normalized variance for th_l     [-]
      alpha_rt,      & ! Factor relating to normalized variance for r_t      [-]
      alpha_u,       & ! Factor relating to normalized variance for u        [-]
      alpha_v,       & ! Factor relating to normalized variance for v        [-]
      crt_1,         & ! Coef. on r_t in s/t eqns. (1st PDF comp.)           [-]
      crt_2,         & ! Coef. on r_t in s/t eqns. (2nd PDF comp.)           [-]
      cthl_1,        & ! Coef. on th_l in s/t eqns. (1st PDF comp.)  [(kg/kg)/K]
      cthl_2           ! Coef. on th_l in s/t eqns. (2nd PDF comp.)  [(kg/kg)/K]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      chi_1,           & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      chi_2,           & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      stdev_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      stdev_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      stdev_eta_1,     & ! Standard dev. of eta (old t) (1st PDF comp.)  [kg/kg]
      stdev_eta_2,     & ! Standard dev. of eta (old t) (2nd PDF comp.)  [kg/kg]
      covar_chi_eta_1, & ! Covariance of chi and eta (1st PDF comp.) [kg^2/kg^2]
      covar_chi_eta_2, & ! Covariance of chi and eta (2nd PDF comp.) [kg^2/kg^2]
      corr_w_chi_1,    & ! Correlation of w and chi (1st PDF component)      [-]
      corr_w_chi_2,    & ! Correlation of w and chi (2nd PDF component)      [-]
      corr_w_eta_1,    & ! Correlation of w and eta (1st PDF component)      [-]
      corr_w_eta_2,    & ! Correlation of w and eta (2nd PDF component)      [-]
      corr_u_w_1,      & ! Correlation of u and w   (1st PDF component)      [-]
      corr_u_w_2,      & ! Correlation of u and w   (2nd PDF component)      [-]
      corr_v_w_1,      & ! Correlation of v and w   (1st PDF component)      [-]
      corr_v_w_2,      & ! Correlation of v and w   (2nd PDF component)      [-]
      corr_chi_eta_1,  & ! Correlation of chi and eta (1st PDF component)    [-]
      corr_chi_eta_2,  & ! Correlation of chi and eta (2nd PDF component)    [-]
      rsatl_1,         & ! Mean of r_sl (1st PDF component)              [kg/kg]
      rsatl_2,         & ! Mean of r_sl (2nd PDF component)              [kg/kg]
      rc_1,            & ! Mean of r_c (1st PDF component)               [kg/kg]
      rc_2,            & ! Mean of r_c (2nd PDF component)               [kg/kg]
      cloud_frac_1,    & ! Cloud fraction (1st PDF component)                [-]
      cloud_frac_2,    & ! Cloud fraction (2nd PDF component)                [-]
      mixt_frac          ! Weight of 1st PDF component (Sk_w dependent)      [-]

    ! Note:  alpha coefficients = 0.5 * ( 1 - correlations^2 ).
    !        These are used to calculate the scalar widths
    !        varnce_thl_1, varnce_thl_2, varnce_rt_1, and varnce_rt_2 as in
    !        Eq. (34) of Larson and Golaz (2005)

    ! Passive scalar local variables

    real( kind = core_rknd ), dimension(gr%nz,sclr_dim) ::  & 
      sclr1, sclr2,  &
      varnce_sclr1, varnce_sclr2, & 
      alpha_sclr,  & 
      corr_sclr_thl_1, corr_sclr_thl_2, &
      corr_sclr_rt_1, corr_sclr_rt_2, &
      corr_w_sclr_1, corr_w_sclr_2

    logical :: &
      l_scalar_calc, &  ! True if sclr_dim > 0
      l_calc_ice_supersat_frac ! True if we should calculate ice_supersat_frac

    ! Quantities needed to predict higher order moments
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      tl1, tl2

    real( kind = core_rknd ), dimension(gr%nz) :: &
      sqrt_wp2, & ! Square root of wp2          [m/s]
      Skthl,    & ! Skewness of thl             [-]
      Skrt        ! Skewness of rt              [-]

    ! Thermodynamic quantity

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      rc_coef    ! Coefficient on X'r_c' in X'th_v' equation    [K/(kg/kg)]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      wprcp_contrib_comp_1,   & ! <w'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
      wprcp_contrib_comp_2,   & ! <w'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]
      wp2rcp_contrib_comp_1,  & ! <w'^2rc'> contrib. (1st comp) [m^2/s^2(kg/kg)]
      wp2rcp_contrib_comp_2,  & ! <w'^2rc'> contrib. (2nd comp) [m^2/s^2(kg/kg)]
      rtprcp_contrib_comp_1,  & ! <rt'rc'> contrib. (1st PDF comp.)  [kg^2/kg^2]
      rtprcp_contrib_comp_2,  & ! <rt'rc'> contrib. (2nd PDF comp.)  [kg^2/kg^2]
      thlprcp_contrib_comp_1, & ! <thl'rc'> contrib. (1st PDF comp.)  [K(kg/kg)]
      thlprcp_contrib_comp_2, & ! <thl'rc'> contrib. (2nd PDF comp.)  [K(kg/kg)]
      uprcp_contrib_comp_1,   & ! <u'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
      uprcp_contrib_comp_2,   & ! <u'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]
      vprcp_contrib_comp_1,   & ! <v'rc'> contrib. (1st PDF comp.)  [m/s(kg/kg)]
      vprcp_contrib_comp_2      ! <v'rc'> contrib. (2nd PDF comp.)  [m/s(kg/kg)]

    ! variables for computing ice cloud fraction
    real( kind = core_rknd), dimension(gr%nz) :: &
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      ice_supersat_frac_2, & ! Ice supersaturation fraction (2nd PDF comp.)  [-]
      rc_1_ice, rc_2_ice
    
    ! To test pdf parameters
    real( kind = core_rknd ), dimension(gr%nz) ::  &
    wm_clubb_pdf,    &
    rtm_clubb_pdf,   &
    thlm_clubb_pdf,  &
    wp2_clubb_pdf,   &
    rtp2_clubb_pdf,  &
    thlp2_clubb_pdf, &
    wp3_clubb_pdf,   &
    rtp3_clubb_pdf,  &
    thlp3_clubb_pdf, &
    Skw_clubb_pdf,   &
    Skrt_clubb_pdf,  &
    Skthl_clubb_pdf

    real( kind = core_rknd ) :: &
      chi_at_ice_sat1, chi_at_ice_sat2

    logical, parameter :: &
      l_liq_ice_loading_test = .false. ! Temp. flag liq./ice water loading test

    integer :: k, i, hm_idx   ! Indices

    ! Value of chi at saturation for liquid water; always 0
    real ( kind = core_rknd ), parameter :: chi_at_liq_sat = 0.0_core_rknd

#ifdef GFDL
    real ( kind = core_rknd ), parameter :: t1_combined = 273.16, &
                                            t2_combined = 268.16, &
                                            t3_combined = 238.16 
#endif

!------------------------ Code Begins ----------------------------------

    ! Check whether the passive scalars are present.
    if ( sclr_dim > 0 ) then
      l_scalar_calc = .true.
    else
      l_scalar_calc = .false.
    end if


    ! Initialize to default values to prevent a runtime error
    if ( ( iiPDF_type /= iiPDF_ADG1 ) .and. ( iiPDF_type /= iiPDF_ADG2 ) ) then

        alpha_thl = one_half
        alpha_rt = one_half

        if ( l_scalar_calc ) then
            alpha_sclr = one_half
        endif

        ! This allows for skewness to be clipped locally without passing the updated
        ! value back out.
        Skrt = Skrt_in
        Skthl = Skthl_in

    endif


    ! Initialize to 0 to prevent a runtime error
    if ( iiPDF_type /= iiPDF_new ) then

       do k = 1, gr%nz

           pdf_implicit_coefs_terms(k)%coef_wp4_implicit = zero
           pdf_implicit_coefs_terms(k)%coef_wprtp2_implicit = zero
           pdf_implicit_coefs_terms(k)%coef_wpthlp2_implicit = zero
           pdf_implicit_coefs_terms(k)%coef_wp2rtp_implicit = zero
           pdf_implicit_coefs_terms(k)%term_wp2rtp_explicit = zero
           pdf_implicit_coefs_terms(k)%coef_wp2thlp_implicit = zero
           pdf_implicit_coefs_terms(k)%term_wp2thlp_explicit = zero
           pdf_implicit_coefs_terms(k)%coef_wprtpthlp_implicit = zero
           pdf_implicit_coefs_terms(k)%term_wprtpthlp_explicit = zero

        end do

        F_w = zero
        F_rt = zero
        F_thl = zero
        min_F_w = zero
        max_F_w = zero
        min_F_rt = zero
        max_F_rt = zero
        min_F_thl = zero
        max_F_thl = zero

    endif


    ! To avoid recomputing
    sqrt_wp2 = sqrt( wp2 )

    ! Select the PDF closure method for the two-component PDF used by CLUBB for
    ! w, rt, theta-l, and passive scalar variables.
    ! Calculate the mixture fraction for the multivariate PDF, as well as both
    ! PDF component means and both PDF component variances for each of w, rt,
    ! theta-l, and passive scalar variables.
    if ( iiPDF_type == iiPDF_ADG1 ) then ! use ADG1

       call ADG1_pdf_driver( wm, rtm, thlm, um, vm,                   & ! In
                             wp2, rtp2, thlp2, up2, vp2,              & ! In
                             Skw, wprtp, wpthlp, upwp, vpwp, sqrt_wp2,& ! In
                             sigma_sqd_w, mixt_frac_max_mag,          & ! In
                             sclrm, sclrp2, wpsclrp, l_scalar_calc,   & ! In
                             w_1, w_2, rt_1, rt_2, thl_1, thl_2,      & ! Out
                             u_1, u_2, v_1, v_2,                      & ! Out
                             varnce_w_1, varnce_w_2, varnce_rt_1,     & ! Out
                             varnce_rt_2, varnce_thl_1, varnce_thl_2, & ! Out
                             varnce_u_1, varnce_u_2,                  & ! Out
                             varnce_v_1, varnce_v_2,                  & ! Out
                             mixt_frac, alpha_rt, alpha_thl,          & ! Out
                             alpha_u, alpha_v,                        & ! Out
                             sclr1, sclr2, varnce_sclr1,              & ! Out
                             varnce_sclr2, alpha_sclr )                 ! Out

    elseif ( iiPDF_type == iiPDF_ADG2 ) then ! use ADG2

       call ADG2_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2,         & ! In
                             Skw, wprtp, wpthlp, sqrt_wp2,            & ! In
                             sclrm, sclrp2, wpsclrp, l_scalar_calc,   & ! In
                             w_1, w_2, rt_1, rt_2, thl_1, thl_2,      & ! Out
                             varnce_w_1, varnce_w_2, varnce_rt_1,     & ! Out
                             varnce_rt_2, varnce_thl_1, varnce_thl_2, & ! Out
                             mixt_frac, alpha_rt, alpha_thl,          & ! Out
                             sigma_sqd_w, sclr1, sclr2,               & ! Out
                             varnce_sclr1, varnce_sclr2, alpha_sclr )   ! Out

    elseif ( iiPDF_type == iiPDF_3D_Luhar ) then ! use 3D Luhar

       call Luhar_3D_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2,      & ! In
                                 Skw, Skrt, Skthl, wprtp, wpthlp,      & ! In
                                 w_1, w_2, rt_1, rt_2, thl_1, thl_2,   & ! Out
                                 varnce_w_1, varnce_w_2, varnce_rt_1,  & ! Out
                                 varnce_rt_2, varnce_thl_1,            & ! Out
                                 varnce_thl_2, mixt_frac )               ! Out

    elseif ( iiPDF_type == iiPDF_new ) then ! use new PDF

       call new_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2, Skw,    & ! In
                            wprtp, wpthlp, rtpthlp,                  & ! In
                            Skrt, Skthl,                             & ! In/Out
                            w_1, w_2, rt_1, rt_2,                    & ! Out
                            thl_1, thl_2, varnce_w_1,                & ! Out
                            varnce_w_2, varnce_rt_1,                 & ! Out
                            varnce_rt_2, varnce_thl_1,               & ! Out
                            varnce_thl_2, mixt_frac,                 & ! Out
                            pdf_implicit_coefs_terms,                & ! Out
                            F_w, F_rt, F_thl, min_F_w, max_F_w,      & ! Out
                            min_F_rt, max_F_rt, min_F_thl, max_F_thl ) ! Out

    elseif ( iiPDF_type == iiPDF_TSDADG ) then

       call tsdadg_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2, & ! In
                               Skw, Skrt, Skthl, wprtp, wpthlp, & ! In
                               w_1, w_2, rt_1, rt_2,            & ! Out
                               thl_1, thl_2, varnce_w_1,        & ! Out
                               varnce_w_2, varnce_rt_1,         & ! Out
                               varnce_rt_2, varnce_thl_1,       & ! Out
                               varnce_thl_2, mixt_frac          ) ! Out

    elseif ( iiPDF_type == iiPDF_LY93 ) then ! use LY93

       call LY93_driver( wm, rtm, thlm, wp2, rtp2,  & ! In
                         thlp2, Skw, Skrt, Skthl,   & ! In
                         w_1, w_2, rt_1, rt_2,      & ! Out
                         thl_1, thl_2, varnce_w_1,  & ! Out
                         varnce_w_2, varnce_rt_1,   & ! Out
                         varnce_rt_2, varnce_thl_1, & ! Out
                         varnce_thl_2, mixt_frac    ) ! Out

    endif ! iiPDF_type


    ! Calculate the PDF component correlations of rt and thl.
    call calc_comp_corrs_binormal( rtpthlp, rtm, thlm, rt_1, rt_2, & ! In
                                   thl_1, thl_2, varnce_rt_1,      & ! In
                                   varnce_rt_2, varnce_thl_1,      & ! In
                                   varnce_thl_2, mixt_frac,        & ! In
                                   corr_rt_thl_1, corr_rt_thl_2    ) ! Out

    if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 ) then

       ! ADG1 and ADG2 define corr_w_rt_1, corr_w_rt_2, corr_w_thl_1, and
       ! corr_w_thl_2 to all have a value of 0, so skip the calculation.
       corr_w_rt_1  = zero
       corr_w_rt_2  = zero
       corr_w_thl_1 = zero
       corr_w_thl_2 = zero

    else

       ! Calculate the PDF component correlations of w and rt.
       call calc_comp_corrs_binormal( wprtp, wm, rtm, w_1, w_2, & ! In
                                      rt_1, rt_2, varnce_w_1,   & ! In
                                      varnce_w_2, varnce_rt_1,  & ! In
                                      varnce_rt_2, mixt_frac,   & ! In
                                      corr_w_rt_1, corr_w_rt_2  ) ! Out

       ! Calculate the PDF component correlations of w and thl.
       call calc_comp_corrs_binormal( wpthlp, wm, thlm, w_1, w_2, & ! In
                                      thl_1, thl_2, varnce_w_1,   & ! In
                                      varnce_w_2, varnce_thl_1,   & ! In
                                      varnce_thl_2, mixt_frac,    & ! In
                                      corr_w_thl_1, corr_w_thl_2  ) ! Out
    endif


    if ( l_scalar_calc ) then

       do i = 1, sclr_dim

          ! Calculate the PDF component correlations of a passive scalar and thl.
          call calc_comp_corrs_binormal( sclrpthlp(:,i), sclrm(:,i),      & ! In
                                         thlm, sclr1(:,i), sclr2(:,i),    & ! In
                                         thl_1, thl_2, varnce_sclr1(:,i), & ! In
                                         varnce_sclr2(:,i), varnce_thl_1, & ! In
                                         varnce_thl_2, mixt_frac,         & ! In
                                         corr_sclr_thl_1(:,i),            & ! Out
                                         corr_sclr_thl_2(:,i)             ) ! Out

          ! Calculate the PDF component correlations of a passive scalar and rt.
          call calc_comp_corrs_binormal( sclrprtp(:,i), sclrm(:,i), rtm, & ! In
                                         sclr1(:,i), sclr2(:,i),         & ! In
                                         rt_1, rt_2, varnce_sclr1(:,i),  & ! In
                                         varnce_sclr2(:,i), varnce_rt_1, & ! In
                                         varnce_rt_2, mixt_frac,         & ! In
                                         corr_sclr_rt_1(:,i),            & ! Out
                                         corr_sclr_rt_2(:,i)             ) ! Out

          if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 ) then

            ! ADG1 and ADG2 define all PDF component correlations involving w
            ! to have a value of 0, so skip the calculation.
            corr_w_sclr_1(:,i) = zero
            corr_w_sclr_2(:,i) = zero

          else

            ! Calculate the PDF component correlations of w and a passive
            ! scalar.
            call calc_comp_corrs_binormal( wpsclrp(:,i), wm, sclrm(:,i),  & ! In
                                           w_1, w_2, sclr1(:,i),          & ! In
                                           sclr2(:,i), varnce_w_1,        & ! In
                                           varnce_w_2, varnce_sclr1(:,i), & ! In
                                           varnce_sclr2(:,i), mixt_frac,  & ! In
                                           corr_w_sclr_1(:,i),            & ! Out
                                           corr_w_sclr_2(:,i)             ) ! Out
          endif

       enddo

    endif


    ! Compute higher order moments (these are interactive)
    wp2rtp = calc_wp2xp_pdf( wm, rtm, w_1, w_2, rt_1, rt_2, varnce_w_1, &
                             varnce_w_2, varnce_rt_1, varnce_rt_2, &
                             corr_w_rt_1, corr_w_rt_2, mixt_frac )

    wp2thlp = calc_wp2xp_pdf( wm, thlm, w_1, w_2, thl_1, thl_2, varnce_w_1, &
                              varnce_w_2, varnce_thl_1, varnce_thl_2, &
                              corr_w_thl_1, corr_w_thl_2, mixt_frac )


    ! Compute higher order moments (these may be interactive)
    if ( l_explicit_turbulent_adv_wp3 .or. iwp4 > 0 ) then
       wp4 = calc_wp4_pdf( wm, w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac )
    endif

    if ( l_explicit_turbulent_adv_xpyp .or. iwprtp2 > 0 ) then
       wprtp2 = calc_wpxp2_pdf( wm, rtm, w_1, w_2, rt_1, rt_2, varnce_w_1, &
                                varnce_w_2, varnce_rt_1, varnce_rt_2, &
                                corr_w_rt_1, corr_w_rt_2, mixt_frac )
    endif

    if ( l_explicit_turbulent_adv_xpyp .or. iwpthlp2 > 0 ) then
       wpthlp2 = calc_wpxp2_pdf( wm, thlm, w_1, w_2, thl_1, thl_2, varnce_w_1, &
                                 varnce_w_2, varnce_thl_1, varnce_thl_2, &
                                 corr_w_thl_1, corr_w_thl_2, mixt_frac )
    endif

    if ( l_explicit_turbulent_adv_xpyp .or. iwprtpthlp > 0 ) then
       wprtpthlp = calc_wpxpyp_pdf( wm, rtm, thlm, w_1, w_2, rt_1, rt_2, &
                                    thl_1, thl_2, varnce_w_1, varnce_w_2, &
                                    varnce_rt_1, varnce_rt_2, varnce_thl_1, &
                                    varnce_thl_2, corr_w_rt_1, corr_w_rt_2, &
                                    corr_w_thl_1, corr_w_thl_2, corr_rt_thl_1, &
                                    corr_rt_thl_2, mixt_frac )
    endif


    ! Scalar Addition to higher order moments
    if ( l_scalar_calc ) then

       do i = 1, sclr_dim

          wp2sclrp(:,i) &
          = calc_wp2xp_pdf( wm, sclrm(:,i), w_1, w_2, sclr1(:,i), &
                            sclr2(:,i), varnce_w_1, varnce_w_2, &
                            varnce_sclr1(:,i), varnce_sclr2(:,i), &
                            corr_w_sclr_1(:,i), corr_w_sclr_2(:,i), &
                            mixt_frac )

          wpsclrp2(:,i) &
          = calc_wpxp2_pdf( wm, sclrm(:,i), w_1, w_2, sclr1(:,i), &
                            sclr2(:,i), varnce_w_1, varnce_w_2, &
                            varnce_sclr1(:,i), varnce_sclr2(:,i), &
                            corr_w_sclr_1(:,i), corr_w_sclr_2(:,i), &
                            mixt_frac )

          wpsclrprtp(:,i) &
          = calc_wpxpyp_pdf( wm, sclrm(:,i), rtm, w_1, w_2, &
                             sclr1(:,i), sclr2(:,i), rt_1, rt_2, &
                             varnce_w_1, varnce_w_2, &
                             varnce_sclr1(:,i), varnce_sclr2(:,i), &
                             varnce_rt_1, varnce_rt_2, &
                             corr_w_sclr_1(:,i), corr_w_sclr_2(:,i), &
                             corr_w_rt_1, corr_w_rt_2, &
                             corr_sclr_rt_1(:,i), &
                             corr_sclr_rt_2(:,i), mixt_frac )

          wpsclrpthlp(:,i) &
          = calc_wpxpyp_pdf( wm, sclrm(:,i), thlm, w_1, w_2, &
                             sclr1(:,i), sclr2(:,i), thl_1, &
                             thl_2, varnce_w_1, varnce_w_2, &
                             varnce_sclr1(:,i), varnce_sclr2(:,i), &
                             varnce_thl_1, varnce_thl_2, &
                             corr_w_sclr_1(:,i), &
                             corr_w_sclr_2(:,i), corr_w_thl_1, &
                             corr_w_thl_2, corr_sclr_thl_1(:,i), &
                             corr_sclr_thl_2(:,i), mixt_frac )

       enddo

    endif


    ! Compute higher order moments that include theta_v.

    ! First compute some preliminary quantities.
    ! "1" denotes first Gaussian; "2" denotes 2nd Gaussian
    ! liq water temp (Sommeria & Deardorff 1977 (SD), eqn. 3)

    tl1  = thl_1*exner
    tl2  = thl_2*exner

#ifdef GFDL
    if ( sclr_dim > 0  .and.  (.not. do_liquid_only_in_clubb) ) then ! h1g, 2010-06-16 begin mod

       where ( tl1 > t1_combined )
          rsatl_1 = sat_mixrat_liq( p_in_Pa, tl1 )
       elsewhere ( tl1 > t2_combined )
          rsatl_1 = sat_mixrat_liq( p_in_Pa, tl1 ) &
                    * (tl1 - t2_combined)/(t1_combined - t2_combined) &
                    + sat_mixrat_ice( p_in_Pa, tl1 ) &
                      * (t1_combined - tl1)/(t1_combined - t2_combined)
       elsewhere ( tl1 > t3_combined )
          rsatl_1 = sat_mixrat_ice( p_in_Pa, tl1 ) &
                    + sat_mixrat_ice( p_in_Pa, tl1 ) * (RH_crit(1, 1) -one ) &
                      * ( t2_combined -tl1)/(t2_combined - t3_combined)
       elsewhere
          rsatl_1 = sat_mixrat_ice( p_in_Pa, tl1 ) * RH_crit(1, 1)
       endwhere

       where ( tl2 > t1_combined )
          rsatl_2 = sat_mixrat_liq( p_in_Pa, tl2 )
       elsewhere ( tl2 > t2_combined )
          rsatl_2 = sat_mixrat_liq( p_in_Pa, tl2 ) &
                    * (tl2 - t2_combined)/(t1_combined - t2_combined) &
                    + sat_mixrat_ice( p_in_Pa, tl2 ) &
                      * (t1_combined - tl2)/(t1_combined - t2_combined)
       elsewhere ( tl2 > t3_combined )
          rsatl_2 = sat_mixrat_ice( p_in_Pa, tl2 ) &
                    + sat_mixrat_ice( p_in_Pa, tl2 )* (RH_crit(1, 2) -one) &
                      * ( t2_combined -tl2)/(t2_combined - t3_combined)
       elsewhere
          rsatl_2 = sat_mixrat_ice( p_in_Pa, tl2 ) * RH_crit(1, 2)
       endwhere

    else ! sclr_dim <= 0  or  do_liquid_only_in_clubb = .T.

       rsatl_1 = sat_mixrat_liq( p_in_Pa, tl1 )
       rsatl_2 = sat_mixrat_liq( p_in_Pa, tl2 )

    endif !sclr_dim > 0
    
    ! Determine whether to compute ice_supersat_frac. We do not compute
    ! ice_supersat_frac for GFDL (unless do_liquid_only_in_clubb is true),
    ! because liquid and ice are both fed into rtm, ruining the calculation.
    if (do_liquid_only_in_clubb) then
      l_calc_ice_supersat_frac = .true.
    else
      l_calc_ice_supersat_frac = .false.
    end if

#else
    rsatl_1 = sat_mixrat_liq( p_in_Pa, tl1 )
    rsatl_2 = sat_mixrat_liq( p_in_Pa, tl2 ) ! h1g, 2010-06-16 end mod

    l_calc_ice_supersat_frac = .true.
#endif


    call transform_pdf_chi_eta_component( tl1, rsatl_1, rt_1, exner, &            ! Intent(in)
                                      varnce_thl_1, varnce_rt_1, corr_rt_thl_1, & ! Intent(in)
                                      chi_1, crt_1, cthl_1, &                     ! Intent(out)
                                      stdev_chi_1, stdev_eta_1, &                 ! Intent(out)
                                      covar_chi_eta_1, corr_chi_eta_1 )           ! Intent(out)
    
    ! Calculate cloud fraction component for pdf 1
    call calc_cloud_frac_component( chi_1, stdev_chi_1, &   ! Intent(in)
                                    chi_at_liq_sat, &       ! Intent(in)
                                    cloud_frac_1, rc_1 )    ! Intent(out)

    
    call transform_pdf_chi_eta_component( tl2, rsatl_2, rt_2, exner, &            ! Intent(in)
                                      varnce_thl_2, varnce_rt_2, corr_rt_thl_2, & ! Intent(in)
                                      chi_2, crt_2, cthl_2, &                     ! Intent(out)
                                      stdev_chi_2, stdev_eta_2, &                 ! Intent(out)
                                      covar_chi_eta_2, corr_chi_eta_2 )           ! Intent(out)

    
    ! Calculate cloud fraction component for pdf 2
    call calc_cloud_frac_component( chi_2, stdev_chi_2, &   ! Intent(in)
                                    chi_at_liq_sat, &       ! Intent(in)
                                    cloud_frac_2, rc_2 )    ! Intent(out)

    ! Calc ice_supersat_frac
    if ( l_calc_ice_supersat_frac ) then

        do i = 1, gr%nz

            if ( tl1(i) <= T_freeze_K ) then
    
                ! Temperature is freezing, we must compute chi_at_ice_sat and 
                ! calculate the a new cloud_frac_component
                chi_at_ice_sat1 = ( sat_mixrat_ice( p_in_Pa(i), tl1(i) ) - rsatl_1(i) ) * crt_1(i)

                call calc_cloud_frac_component( chi_1(i), stdev_chi_1(i), &
                                                chi_at_ice_sat1, &
                                                ice_supersat_frac_1(i), rc_1_ice(i) )
            else

                ! Temperature is warmer than freezing, the ice_supersat_frac calculation is 
                ! the same as cloud_frac
                ice_supersat_frac_1(i) = cloud_frac_1(i)
                rc_1_ice(i) = rc_1(i)

            end if

        end do


        do i = 1, gr%nz

            if ( tl2(i) <= T_freeze_K ) then

                ! Temperature is freezing, we must compute chi_at_ice_sat and 
                ! calculate the a new cloud_frac_component
                chi_at_ice_sat2 = ( sat_mixrat_ice( p_in_Pa(i), tl2(i) ) - rsatl_2(i) ) * crt_2(i)

                call calc_cloud_frac_component( chi_2(i), stdev_chi_2(i), &
                                                chi_at_ice_sat2, &
                                                ice_supersat_frac_2(i), rc_2_ice(i) )
            else

                ! Temperature is warmer than freezing, the ice_supersat_frac calculation is 
                ! the same as cloud_frac
                ice_supersat_frac_2(i) = cloud_frac_2(i)
                rc_2_ice(i) = rc_2(i)

            end if

        end do

        ! Compute ice cloud fraction, ice_supersat_frac
        ice_supersat_frac = mixt_frac * ice_supersat_frac_1 &
                            + ( one - mixt_frac ) * ice_supersat_frac_2

    else 

        ! ice_supersat_frac will be garbage if computed as above
        ice_supersat_frac = 0.0_core_rknd

        if (clubb_at_least_debug_level( 1 )) then
            write(fstderr,*) "Warning: ice_supersat_frac has garbage values if &
                            & do_liquid_only_in_clubb = .false."
        end if

    end if ! l_calc_ice_supersat_frac


    ! Compute cloud fraction and mean cloud water mixing ratio.
    ! Reference:
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_cloud_terms
    cloud_frac = mixt_frac * cloud_frac_1 + ( one - mixt_frac ) * cloud_frac_2
    rcm = mixt_frac * rc_1 + ( one - mixt_frac ) * rc_2
    rcm = max( zero_threshold, rcm )

    if ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 ) then

        ! corr_w_rt and corr_w_thl are zero for these pdf types so
        ! corr_w_chi and corr_w_eta are zero as well
        corr_w_chi_1 = zero
        corr_w_chi_2 = zero
        corr_w_eta_1 = zero
        corr_w_eta_2 = zero
        corr_u_w_1   = zero
        corr_u_w_2   = zero
        corr_v_w_1   = zero
        corr_v_w_2   = zero

    else 
        
        ! Correlation of w and chi for each component.
        corr_w_chi_1 &
        = calc_corr_chi_x( crt_1, cthl_1, sqrt(varnce_rt_1), sqrt(varnce_thl_1), &
                           stdev_chi_1, corr_w_rt_1, corr_w_thl_1 )

        corr_w_chi_2 &
        = calc_corr_chi_x( crt_2, cthl_2, sqrt(varnce_rt_2), sqrt(varnce_thl_2), &
                           stdev_chi_2, corr_w_rt_2, corr_w_thl_2 )

        ! Correlation of w and eta for each component.
        corr_w_eta_1 &
        = calc_corr_eta_x( crt_1, cthl_1, sqrt(varnce_rt_1), sqrt(varnce_thl_1), &
                           stdev_eta_1, corr_w_rt_1, corr_w_thl_1 )

        corr_w_eta_2 &
        = calc_corr_eta_x( crt_2, cthl_2, sqrt(varnce_rt_2), sqrt(varnce_thl_2),  &
                           stdev_eta_2, corr_w_rt_2, corr_w_thl_2 )

    end if


    ! Compute moments that depend on theta_v
    ! 
    ! The moments that depend on th_v' are calculated based on an approximated
    ! and linearized form of the theta_v equation:
    ! 
    ! theta_v = theta_l + { (R_v/R_d) - 1 } * thv_ds * r_t
    !                   + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c;
    ! 
    ! and therefore:
    ! 
    ! th_v' = th_l' + { (R_v/R_d) - 1 } * thv_ds * r_t'
    !               + [ {L_v/(C_p*exner)} - (R_v/R_d) * thv_ds ] * r_c';
    ! 
    ! where thv_ds is used as a reference value to approximate theta_l.
    ! 
    ! Reference:
    ! https://arxiv.org/pdf/1711.03675v1.pdf#nameddest=url:anl_int_buoy_terms
    
    ! Calculate the contributions to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
    ! <thl'rc'> from the 1st PDF component.
    call calc_xprcp_component( wm, rtm, thlm, um, vm, rcm,        & ! In
                               w_1, rt_1, thl_1, u_1, v_1,        & ! In
                               varnce_w_1,                        & ! In
                               chi_1, stdev_chi_1, stdev_eta_1,   & ! In
                               corr_w_chi_1, corr_chi_eta_1,      & ! In
                               corr_u_w_1, corr_v_w_1,            & ! In
                               crt_1, cthl_1, rc_1, cloud_frac_1, & ! In
                               wprcp_contrib_comp_1,              & ! Out
                               wp2rcp_contrib_comp_1,             & ! Out
                               rtprcp_contrib_comp_1,             & ! Out
                               thlprcp_contrib_comp_1,            & ! Out
                               uprcp_contrib_comp_1,              & ! Out
                               vprcp_contrib_comp_1               ) ! Out

    call calc_xprcp_component( wm, rtm, thlm, um, vm, rcm,        & ! In
                               w_2, rt_2, thl_2, u_2, v_2,        & ! In
                               varnce_w_2,                        & ! In
                               chi_2, stdev_chi_2, stdev_eta_2,   & ! In
                               corr_w_chi_2, corr_chi_eta_2,      & ! In
                               corr_u_w_2, corr_v_w_2,            & ! In
                               crt_2, cthl_2, rc_2, cloud_frac_2, & ! In
                               wprcp_contrib_comp_2,              & ! Out
                               wp2rcp_contrib_comp_2,             & ! Out
                               rtprcp_contrib_comp_2,             & ! Out
                               thlprcp_contrib_comp_2,            & ! Out
                               uprcp_contrib_comp_2,              & ! Out
                               vprcp_contrib_comp_2               ) ! Out

    
    ! Calculate rc_coef, which is the coefficient on <x'rc'> in the <x'thv'> equation.
    rc_coef = Lv / ( exner * Cp ) - ep2 * thv_ds


    ! Calculate <w'rc'>, <w'^2 rc'>, <rt'rc'>, and <thl'rc'>.
    wprcp = mixt_frac * wprcp_contrib_comp_1 &
            + ( one - mixt_frac ) * wprcp_contrib_comp_2

    wp2rcp = mixt_frac * wp2rcp_contrib_comp_1 &
             + ( one - mixt_frac ) * wp2rcp_contrib_comp_2

    rtprcp = mixt_frac * rtprcp_contrib_comp_1 & 
             + ( one - mixt_frac ) * rtprcp_contrib_comp_2

    thlprcp = mixt_frac * thlprcp_contrib_comp_1 &
              + ( one - mixt_frac ) * thlprcp_contrib_comp_2

    uprcp = mixt_frac * uprcp_contrib_comp_1 &
            + ( one - mixt_frac ) * uprcp_contrib_comp_2

    vprcp = mixt_frac * vprcp_contrib_comp_1 &
            + ( one - mixt_frac ) * vprcp_contrib_comp_2

    ! Calculate <w'thv'>, <w'^2 thv'>, <rt'thv'>, and <thl'thv'>.
    wpthvp = wpthlp + ep1 * thv_ds * wprtp + rc_coef * wprcp

    wp2thvp = wp2thlp + ep1 * thv_ds * wp2rtp + rc_coef * wp2rcp

    rtpthvp = rtpthlp + ep1 * thv_ds * rtp2 + rc_coef * rtprcp

    thlpthvp = thlp2 + ep1 * thv_ds * rtpthlp + rc_coef * thlprcp


    ! Add the precipitation loading term in the <x'thv'> equation.
    if ( l_liq_ice_loading_test ) then

       do hm_idx = 1, hydromet_dim, 1

          if ( l_mix_rat_hm(hm_idx) ) then

             wp2thvp  = wp2thvp  - thv_ds * wp2hmp(:,hm_idx)
             wpthvp   = wpthvp   - thv_ds * wphydrometp(:,hm_idx)
             thlpthvp = thlpthvp - thv_ds * thlphmp(:,hm_idx)
             rtpthvp  = rtpthvp  - thv_ds * rtphmp(:,hm_idx)

          endif

       enddo

    endif


    ! Account for subplume correlation of scalar, theta_v.
    ! See Eqs. A13, A8 from Larson et al. (2002) ``Small-scale...''
    !  where the ``scalar'' in this paper is w.
    if ( l_scalar_calc ) then
       do i = 1, sclr_dim
          sclrprcp(:,i) &
          = mixt_frac * ( ( sclr1(:,i) - sclrm(:,i) ) * rc_1 ) &
            + ( one - mixt_frac ) * ( ( sclr2(:,i) - sclrm(:,i) ) * rc_2 ) &
            + mixt_frac * corr_sclr_rt_1(:,i) * crt_1 &
              * sqrt( varnce_sclr1(:,i) * varnce_rt_1 ) * cloud_frac_1 & 
            + ( one - mixt_frac ) * corr_sclr_rt_2(:,i) * crt_2 &
              * sqrt( varnce_sclr2(:,i) * varnce_rt_2 ) * cloud_frac_2 & 
            - mixt_frac * corr_sclr_thl_1(:,i) * cthl_1 &
              * sqrt( varnce_sclr1(:,i) * varnce_thl_1 ) * cloud_frac_1 & 
            - ( one - mixt_frac ) * corr_sclr_thl_2(:,i) * cthl_2 &
              * sqrt( varnce_sclr2(:,i) * varnce_thl_2 ) * cloud_frac_2

          sclrpthvp(:,i) = sclrpthlp(:,i) + ep1*thv_ds*sclrprtp(:,i) &
                           + rc_coef*sclrprcp(:,i)
       enddo ! i=1, sclr_dim
    endif ! l_scalar_calc

    
    ! Compute variance of liquid water mixing ratio.
    ! This is not needed for closure.  Statistical Analysis only.

#ifndef CLUBB_CAM
    !  if CLUBB is used in CAM we want this variable computed no matter what
    if ( ircp2 > 0 ) then
#endif

      rcp2 = mixt_frac * ( chi_1*rc_1 + cloud_frac_1*stdev_chi_1**2 ) &
             + ( one-mixt_frac ) * ( chi_2*rc_2 + cloud_frac_2*stdev_chi_2**2 ) - rcm**2
      rcp2 = max( zero_threshold, rcp2 )

#ifndef CLUBB_CAM
    !  if CLUBB is used in CAM we want this variable computed no matter what
    end if
#endif


    ! Save PDF parameters
    pdf_params%w_1             = w_1
    pdf_params%w_2             = w_2
    pdf_params%varnce_w_1      = varnce_w_1
    pdf_params%varnce_w_2      = varnce_w_2
    pdf_params%rt_1            = rt_1
    pdf_params%rt_2            = rt_2
    pdf_params%varnce_rt_1     = varnce_rt_1
    pdf_params%varnce_rt_2     = varnce_rt_2
    pdf_params%thl_1           = thl_1
    pdf_params%thl_2           = thl_2
    pdf_params%varnce_thl_1    = varnce_thl_1
    pdf_params%varnce_thl_2    = varnce_thl_2
    pdf_params%corr_w_rt_1     = corr_w_rt_1
    pdf_params%corr_w_rt_2     = corr_w_rt_2
    pdf_params%corr_w_thl_1    = corr_w_thl_1
    pdf_params%corr_w_thl_2    = corr_w_thl_2
    pdf_params%corr_rt_thl_1   = corr_rt_thl_1
    pdf_params%corr_rt_thl_2   = corr_rt_thl_2
    pdf_params%alpha_thl       = alpha_thl
    pdf_params%alpha_rt        = alpha_rt
    pdf_params%crt_1           = crt_1
    pdf_params%crt_2           = crt_2
    pdf_params%cthl_1          = cthl_1
    pdf_params%cthl_2          = cthl_2
    pdf_params%chi_1           = chi_1
    pdf_params%chi_2           = chi_2
    pdf_params%stdev_chi_1     = stdev_chi_1
    pdf_params%stdev_chi_2     = stdev_chi_2
    pdf_params%stdev_eta_1     = stdev_eta_1
    pdf_params%stdev_eta_2     = stdev_eta_2
    pdf_params%covar_chi_eta_1 = covar_chi_eta_1
    pdf_params%covar_chi_eta_2 = covar_chi_eta_2
    pdf_params%corr_w_chi_1    = corr_w_chi_1
    pdf_params%corr_w_chi_2    = corr_w_chi_2
    pdf_params%corr_w_eta_1    = corr_w_eta_1
    pdf_params%corr_w_eta_2    = corr_w_eta_2
    pdf_params%corr_chi_eta_1  = corr_chi_eta_1
    pdf_params%corr_chi_eta_2  = corr_chi_eta_2
    pdf_params%rsatl_1         = rsatl_1
    pdf_params%rsatl_2         = rsatl_2
    pdf_params%rc_1            = rc_1
    pdf_params%rc_2            = rc_2
    pdf_params%cloud_frac_1    = cloud_frac_1
    pdf_params%cloud_frac_2    = cloud_frac_2
    pdf_params%mixt_frac       = mixt_frac
    pdf_params%ice_supersat_frac_1 = ice_supersat_frac_1
    pdf_params%ice_supersat_frac_2 = ice_supersat_frac_2


    if ( clubb_at_least_debug_level( 2 ) ) then

      call pdf_closure_check & 
           ( wp4, wprtp2, wp2rtp, wpthlp2, & 
             wp2thlp, cloud_frac, rcm, wpthvp, wp2thvp, & 
             rtpthvp, thlpthvp, wprcp, wp2rcp, & 
             rtprcp, thlprcp, rcp2, wprtpthlp, & 
             crt_1, crt_2, cthl_1, cthl_2, pdf_params, &
             sclrpthvp, sclrprcp, wpsclrp2, & 
             wpsclrprtp, wpsclrpthlp, wp2sclrp )

      ! Error Reporting
      ! Joshua Fasching February 2008

      if ( err_code == clubb_fatal_error ) then

        write(fstderr,*) "Error in pdf_closure_new"

        write(fstderr,*) "Intent(in)"

        write(fstderr,*) "p_in_Pa = ", p_in_Pa
        write(fstderr,*) "exner = ", exner
        write(fstderr,*) "thv_ds = ", thv_ds
        write(fstderr,*) "wm = ", wm
        write(fstderr,*) "wp2 = ", wp2
        write(fstderr,*) "wp3 = ", wp3
        write(fstderr,*) "sigma_sqd_w = ", sigma_sqd_w
        write(fstderr,*) "rtm = ", rtm
        write(fstderr,*) "rtp2 = ", rtp2
        write(fstderr,*) "wprtp = ", wprtp
        write(fstderr,*) "thlm = ", thlm
        write(fstderr,*) "thlp2 = ", thlp2
        write(fstderr,*) "wpthlp = ", wpthlp
        write(fstderr,*) "rtpthlp = ", rtpthlp

        if ( sclr_dim > 0 ) then
          write(fstderr,*) "sclrm = ", sclrm
          write(fstderr,*) "wpsclrp = ", wpsclrp
          write(fstderr,*) "sclrp2 = ", sclrp2
          write(fstderr,*) "sclrprtp = ", sclrprtp
          write(fstderr,*) "sclrpthlp = ", sclrpthlp
        end if

        write(fstderr,*) "Intent(out)"

        write(fstderr,*) "wp4 = ", wp4
        write(fstderr,*) "wprtp2 = ", wprtp2
        write(fstderr,*) "wp2rtp = ", wp2rtp
        write(fstderr,*) "wpthlp2 = ", wpthlp2
        write(fstderr,*) "cloud_frac = ", cloud_frac
        write(fstderr,*) "ice_supersat_frac = ", ice_supersat_frac
        write(fstderr,*) "rcm = ", rcm
        write(fstderr,*) "wpthvp = ", wpthvp
        write(fstderr,*) "wp2thvp = ", wp2thvp
        write(fstderr,*) "rtpthvp = ", rtpthvp
        write(fstderr,*) "thlpthvp = ", thlpthvp
        write(fstderr,*) "wprcp = ", wprcp
        write(fstderr,*) "wp2rcp = ", wp2rcp
        write(fstderr,*) "rtprcp = ", rtprcp
        write(fstderr,*) "thlprcp = ", thlprcp
        write(fstderr,*) "rcp2 = ", rcp2
        write(fstderr,*) "wprtpthlp = ", wprtpthlp
        write(fstderr,*) "pdf_params%w_1 = ", pdf_params%w_1
        write(fstderr,*) "pdf_params%w_2 = ", pdf_params%w_2
        write(fstderr,*) "pdf_params%varnce_w_1 = ", pdf_params%varnce_w_1
        write(fstderr,*) "pdf_params%varnce_w_2 = ", pdf_params%varnce_w_2
        write(fstderr,*) "pdf_params%rt_1 = ", pdf_params%rt_1
        write(fstderr,*) "pdf_params%rt_2 = ", pdf_params%rt_2
        write(fstderr,*) "pdf_params%varnce_rt_1 = ", pdf_params%varnce_rt_1
        write(fstderr,*) "pdf_params%varnce_rt_2 = ", pdf_params%varnce_rt_2
        write(fstderr,*) "pdf_params%thl_1 = ", pdf_params%thl_1
        write(fstderr,*) "pdf_params%thl_2 = ", pdf_params%thl_2
        write(fstderr,*) "pdf_params%varnce_thl_1 = ", pdf_params%varnce_thl_1
        write(fstderr,*) "pdf_params%varnce_thl_2 = ", pdf_params%varnce_thl_2
        write(fstderr,*) "pdf_params%corr_w_rt_1 = ", pdf_params%corr_w_rt_1
        write(fstderr,*) "pdf_params%corr_w_rt_2 = ", pdf_params%corr_w_rt_2
        write(fstderr,*) "pdf_params%corr_w_thl_1 = ", pdf_params%corr_w_thl_1
        write(fstderr,*) "pdf_params%corr_w_thl_2 = ", pdf_params%corr_w_thl_2
        write(fstderr,*) "pdf_params%corr_rt_thl_1 = ", pdf_params%corr_rt_thl_1
        write(fstderr,*) "pdf_params%corr_rt_thl_2 = ", pdf_params%corr_rt_thl_2
        write(fstderr,*) "pdf_params%alpha_thl = ", pdf_params%alpha_thl
        write(fstderr,*) "pdf_params%alpha_rt = ", pdf_params%alpha_rt
        write(fstderr,*) "pdf_params%crt_1 = ", pdf_params%crt_1
        write(fstderr,*) "pdf_params%crt_2 = ", pdf_params%crt_2
        write(fstderr,*) "pdf_params%cthl_1 = ", pdf_params%cthl_1
        write(fstderr,*) "pdf_params%cthl_2 = ", pdf_params%cthl_2
        write(fstderr,*) "pdf_params%chi_1 = ", pdf_params%chi_1
        write(fstderr,*) "pdf_params%chi_2 = ", pdf_params%chi_2
        write(fstderr,*) "pdf_params%stdev_chi_1 = ", pdf_params%stdev_chi_1
        write(fstderr,*) "pdf_params%stdev_chi_2 = ", pdf_params%stdev_chi_2
        write(fstderr,*) "pdf_params%stdev_eta_1 = ", pdf_params%stdev_eta_1
        write(fstderr,*) "pdf_params%stdev_eta_2 = ", pdf_params%stdev_eta_2
        write(fstderr,*) "pdf_params%covar_chi_eta_1 = ", &
                         pdf_params%covar_chi_eta_1
        write(fstderr,*) "pdf_params%covar_chi_eta_2 = ", &
                         pdf_params%covar_chi_eta_2
        write(fstderr,*) "pdf_params%corr_w_chi_1 = ", pdf_params%corr_w_chi_1
        write(fstderr,*) "pdf_params%corr_w_chi_2 = ", pdf_params%corr_w_chi_2
        write(fstderr,*) "pdf_params%corr_w_eta_1 = ", pdf_params%corr_w_eta_1
        write(fstderr,*) "pdf_params%corr_w_eta_2 = ", pdf_params%corr_w_eta_2
        write(fstderr,*) "pdf_params%corr_chi_eta_1 = ", &
                         pdf_params%corr_chi_eta_1
        write(fstderr,*) "pdf_params%corr_chi_eta_2 = ", &
                         pdf_params%corr_chi_eta_2
        write(fstderr,*) "pdf_params%rsatl_1 = ", pdf_params%rsatl_1
        write(fstderr,*) "pdf_params%rsatl_2 = ", pdf_params%rsatl_2
        write(fstderr,*) "pdf_params%rc_1 = ", pdf_params%rc_1
        write(fstderr,*) "pdf_params%rc_2 = ", pdf_params%rc_2
        write(fstderr,*) "pdf_params%cloud_frac_1 = ", pdf_params%cloud_frac_1
        write(fstderr,*) "pdf_params%cloud_frac_2 = ", pdf_params%cloud_frac_2
        write(fstderr,*) "pdf_params%mixt_frac = ", pdf_params%mixt_frac
        write(fstderr,*) "pdf_params%ice_supersat_frac_1 = ", &
                         pdf_params%ice_supersat_frac_1
        write(fstderr,*) "pdf_params%ice_supersat_frac_2 = ", &
                         pdf_params%ice_supersat_frac_2

        if ( sclr_dim > 0 )then
          write(fstderr,*) "sclrpthvp = ", sclrpthvp
          write(fstderr,*) "sclrprcp = ", sclrprcp
          write(fstderr,*) "wpsclrp2 = ", wpsclrp2
          write(fstderr,*) "wpsclrprtp = ", wpsclrprtp
          write(fstderr,*) "wpsclrpthlp = ", wpsclrpthlp
          write(fstderr,*) "wp2sclrp = ", wp2sclrp
        end if

      end if ! Fatal error

      ! Error check pdf parameters and moments to ensure consistency
      if ( iiPDF_type == iiPDF_3D_Luhar ) then

        ! Means
        wm_clubb_pdf = mixt_frac * w_1 + ( one - mixt_frac ) * w_2

        do k = 1, gr%nz, 1
           if ( abs( ( wm_clubb_pdf(k) - wm(k) ) &
                     / max( wm(k), eps ) ) > .05_core_rknd ) then
              write(fstderr,*) "wm error at thlm = ", thlm(k), &
                               ( ( wm_clubb_pdf(k) - wm(k) ) &
                                 / max( wm(k), eps ) )
           endif
        enddo ! k = 1, gr%nz, 1

        rtm_clubb_pdf = mixt_frac * rt_1 + ( one - mixt_frac ) * rt_2

        do k = 1, gr%nz, 1
           if ( abs( ( rtm_clubb_pdf(k) - rtm(k) ) &
                     / max( rtm(k), eps ) ) > .05_core_rknd ) then
              write(fstderr,*) "rtm error at thlm = ", thlm(k), &
                               ( ( rtm_clubb_pdf(k) - rtm(k) ) &
                                 / max( rtm(k), eps ) )
           endif
        enddo ! k = 1, gr%nz, 1

        thlm_clubb_pdf = mixt_frac * thl_1 + ( one - mixt_frac ) * thl_2

        do k = 1, gr%nz, 1
           if ( abs( ( thlm_clubb_pdf(k) - thlm(k) ) / thlm(k) ) &
                > .05_core_rknd ) then
              write(fstderr,*) "thlm error at thlm = ", thlm(k), &
                               ( ( thlm_clubb_pdf(k) - thlm(k) ) / thlm(k) )
           endif
        enddo ! k = 1, gr%nz, 1

        ! Variances
        wp2_clubb_pdf = mixt_frac * ( ( w_1 - wm )**2 + varnce_w_1 ) &
                        + ( one - mixt_frac ) * ( ( w_2 - wm )**2 + varnce_w_2 )

        do k = 1, gr%nz, 1
           if ( wp2(k) > w_tol**2 ) then
              if ( abs( ( wp2_clubb_pdf(k) - wp2(k) ) / wp2(k) ) &
                   > .05_core_rknd ) then
                 write(fstderr,*) "wp2 error at thlm = ", thlm(k), &
                                  ( ( wp2_clubb_pdf(k) - wp2(k) ) / wp2(k) )
              endif
           endif
        enddo ! k = 1, gr%nz, 1

        rtp2_clubb_pdf &
        = mixt_frac * ( ( rt_1 - rtm )**2 + varnce_rt_1 ) &
          + ( one - mixt_frac ) * ( ( rt_2 - rtm )**2 + varnce_rt_2 )

        do k = 1, gr%nz, 1
           if ( rtp2(k) > rt_tol**2 ) then
              if ( abs( ( rtp2_clubb_pdf(k) - rtp2(k) ) / rtp2(k) ) &
                   > .05_core_rknd ) then
                 write(fstderr,*) "rtp2 error at thlm = ", thlm(k), &
                 "Error = ", ( ( rtp2_clubb_pdf(k) - rtp2(k) ) / rtp2(k) )
              endif
           endif
        enddo ! k = 1, gr%nz, 1

        thlp2_clubb_pdf &
        = mixt_frac * ( ( thl_1 - thlm )**2 + varnce_thl_1 ) &
          + ( one - mixt_frac ) * ( ( thl_2 - thlm )**2 + varnce_thl_2 )

        do k = 1, gr%nz, 1
           if( thlp2(k) > thl_tol**2 ) then
              if ( abs( ( thlp2_clubb_pdf(k) - thlp2(k) ) / thlp2(k) ) &
                   > .05_core_rknd ) then
                 write(fstderr,*) "thlp2 error at thlm = ", thlm(k), &
                 "Error = ", ( ( thlp2_clubb_pdf(k) - thlp2(k) ) / thlp2(k) )
              endif
           endif
        enddo ! k = 1, gr%nz, 1

        ! Third order moments
        wp3_clubb_pdf &
        = mixt_frac * ( w_1 - wm ) &
                    * ( ( w_1 - wm )**2 + three * varnce_w_1 ) &
          + ( one - mixt_frac ) * ( w_2 - wm ) &
                                * ( ( w_2 - wm )**2 + three * varnce_w_2 )

        rtp3_clubb_pdf &
        = mixt_frac * ( rt_1 - rtm ) &
                    * ( ( rt_1 - rtm )**2 + three * varnce_rt_1 ) &
          + ( one - mixt_frac ) * ( rt_2 - rtm ) &
                                * ( ( rt_2 - rtm )**2 + three * varnce_rt_2 )

        thlp3_clubb_pdf &
        = mixt_frac * ( thl_1 - thlm ) &
                    * ( ( thl_1 - thlm )**2 + three * varnce_thl_1 ) &
          + ( one - mixt_frac ) * ( thl_2 - thlm ) &
                                * ( ( thl_2 - thlm )**2 + three * varnce_thl_2 )

        ! Skewness
        Skw_clubb_pdf &
        = wp3_clubb_pdf &
          / ( wp2_clubb_pdf + Skw_denom_coef * w_tol**2 )**1.5_core_rknd

        do k = 1, gr%nz, 1
           if ( Skw(k) > .05_core_rknd ) then
              if( abs( ( Skw_clubb_pdf(k) - Skw(k) ) / Skw(k) ) &
                  > .25_core_rknd ) then
                 write(fstderr,*) "Skw error at thlm = ", thlm(k), &
                 "Error = ", ( ( Skw_clubb_pdf(k) - Skw(k) ) / Skw(k) ), &
                 Skw_clubb_pdf(k), Skw(k)
              endif
           endif
        enddo ! k = 1, gr%nz, 1

        Skrt_clubb_pdf &
        = rtp3_clubb_pdf &
          / ( rtp2_clubb_pdf + Skw_denom_coef * rt_tol**2 )**1.5_core_rknd

        do k = 1, gr%nz, 1
           if ( Skrt(k) > .05_core_rknd ) then
              if( abs( ( Skrt_clubb_pdf(k) - Skrt(k) ) / Skrt(k) ) &
                  > .25_core_rknd ) then
                 write(fstderr,*) "Skrt error at thlm = ", thlm(k), &
                 "Error = ", ( ( Skrt_clubb_pdf(k) - Skrt(k) ) / Skrt(k) ), &
                 Skrt_clubb_pdf(k), Skrt(k)
              endif
           endif
        enddo ! k = 1, gr%nz, 1

        Skthl_clubb_pdf &
        = thlp3_clubb_pdf &
          / ( thlp2_clubb_pdf + Skw_denom_coef * thl_tol**2 )**1.5_core_rknd

        do k = 1, gr%nz, 1
           if ( Skthl(k) > .05_core_rknd ) then
              if ( abs( ( Skthl_clubb_pdf(k) - Skthl(k) ) / Skthl(k) ) &
                   > .25_core_rknd ) then
                 write(fstderr,*) "Skthl error at thlm = ", thlm(k), &
                 "Error = ", ( ( Skthl_clubb_pdf(k) - Skthl(k) ) / Skthl(k) ), &
                 Skthl_clubb_pdf(k), Skthl(k)
              endif
           endif
        enddo ! k = 1, gr%nz, 1

      endif ! iiPDF_type == iiPDF_3D_Luhar

    endif ! clubb_at_least_debug_level


    return

  end subroutine pdf_closure

  !===============================================================================================
  subroutine transform_pdf_chi_eta_component( tl, rsatl, rt, exner, &               ! intent(in)
                                              varnce_thl, varnce_rt, corr_rt_thl, & ! intent(in)
                                              chi, crt, cthl, &                     ! intent(out)
                                              stdev_chi, stdev_eta, &               ! intent(out)
                                              covar_chi_eta, corr_chi_eta )         ! intent(out)
    use grid_class, only: &
        gr    ! Variable type(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use constants_clubb, only: &
        zero, one, two, &
        ep, Lv, Rd, Cp, &
        chi_tol, &
        eta_tol, &
        max_mag_correlation

    implicit none

    ! ----------- Input Variables -----------
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      tl, &
      rsatl, &
      rt, &
      varnce_thl, &
      varnce_rt, &
      corr_rt_thl, &
      exner
    
    ! ----------- Output Variables -----------
    real( kind = core_rknd ), dimension(gr%nz), intent(out) ::  &
      chi, &            ! s from Lewellen and Yoh 1993 (LY) eqn. 1
      crt, &            ! Coefficients for s'
      cthl, &           ! Coefficients for s'
      stdev_chi, &      ! Standard deviation of chi for each component.
      stdev_eta, &      ! Standard deviation of eta for each component.
      covar_chi_eta, &  ! Covariance of chi and eta for each component.
      corr_chi_eta      ! Correlation of chi and eta for each component.

    ! ----------- Local Variables -----------
    real( kind = core_rknd ), dimension(gr%nz) ::  &
      beta

    real( kind = core_rknd ) :: &
      varnce_rt_term, &
      corr_rt_thl_term, &
      varnce_thl_term

    ! Loop variable
    integer :: i

    ! ----------- Begin Code -----------

    ! SD's beta (eqn. 8)
    beta = ep * ( Lv/(Rd*tl) ) * ( Lv/(Cp*tl) )

    ! s from Lewellen and Yoh 1993 (LY) eqn. 1
    chi = ( rt - rsatl ) / ( one + beta * rsatl )

    ! For each normal distribution in the sum of two normal distributions,
    ! s' = crt * rt'  +  cthl * thl';
    ! therefore, x's' = crt * x'rt'  +  cthl * x'thl'.
    ! Larson et al. May, 2001.
    crt  = one / ( one + beta * rsatl )
    cthl = (one + beta * rt) / ( one + beta * rsatl )**2 * ( Cp/Lv ) * beta * rsatl * exner

    ! Calculate covariance, correlation, and standard deviation of 
    ! chi and eta for each component
    ! Include subplume correlation of qt, thl
    do i = 1, gr%nz
       
        varnce_rt_term = crt(i)**2 * varnce_rt(i)
        varnce_thl_term = cthl(i)**2 * varnce_thl(i)

        covar_chi_eta(i) = varnce_rt_term - varnce_thl_term

        corr_rt_thl_term = two * corr_rt_thl(i) * crt(i) * cthl(i) &
                           * sqrt( varnce_rt(i) * varnce_thl(i) )

        stdev_chi(i) = sqrt( varnce_rt_term - corr_rt_thl_term + varnce_thl_term )
        stdev_eta(i) = sqrt( varnce_rt_term + corr_rt_thl_term + varnce_thl_term )
        
    end do

    ! We need to introduce a threshold value for the variance of chi and eta
    where ( stdev_chi < chi_tol .or. stdev_eta < eta_tol )

        where ( stdev_chi < chi_tol ) stdev_chi = zero  ! Treat chi as a delta function
        where ( stdev_eta < eta_tol ) stdev_eta = zero  ! Treat eta as a delta function

        corr_chi_eta = zero

    elsewhere

        corr_chi_eta = covar_chi_eta / ( stdev_chi * stdev_eta )
        corr_chi_eta = min( max_mag_correlation, max( -max_mag_correlation, corr_chi_eta ) )

    end where


  end subroutine transform_pdf_chi_eta_component
  
  !=============================================================================
  function calc_wp4_pdf( wm, w_1, w_2, varnce_w_1, varnce_w_2, mixt_frac ) &
  result( wp4 )

    ! Description:
    ! Calculates <w'^4> by integrating over the PDF of w.  The integral is:
    !
    ! <w'^4> = INT(-inf:inf) ( w - <w> )^4 P(w) dw;
    !
    ! where <w> is the overall mean of w and P(w) is a two-component normal
    ! distribution of w.  The integrated equation is:
    !
    ! <w'^4> = mixt_frac * ( 3 * sigma_w_1^4
    !                        + 6 * ( mu_w_1 - <w> )^2 * sigma_w_1^2
    !                        + ( mu_w_1 - <w> )^4 )
    !          + ( 1 - mixt_frac ) * ( 3 * sigma_w_2^4
    !                                  + 6 * ( mu_w_2 - <w> )^2 * sigma_w_2^2
    !                                  + ( mu_w_2 - <w> )^4 );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, sigma_w_1 is the standard deviation of w in
    ! the 1st PDF component, sigma_w_2 is the standard deviation of w in the 2nd
    ! PDF component, and mixt_frac is the mixture fraction, which is the weight
    ! of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        six,   & ! Variable(s)
        three, &
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      wm,         & ! Mean of w (overall)                           [m/s]
      w_1,        & ! Mean of w (1st PDF component)                 [m/s]
      w_2,        & ! Mean of w (2nd PDF component)                 [m/s]
      varnce_w_1, & ! Variance of w (1st PDF component)             [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)             [m^2/s^2]
      mixt_frac     ! Mixture fraction                              [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      wp4    ! <w'^4>                   [m^4/s^4]


    ! Calculate <w'^4> by integrating over the PDF.
    wp4 = mixt_frac * ( three * varnce_w_1**2 &
                        + six * ( ( w_1 - wm )**2 ) * varnce_w_1 &
                        + ( w_1 - wm )**4 ) & 
          + ( one - mixt_frac ) * ( three * varnce_w_2**2 &
                                    + six * ( (w_2 - wm )**2 )*varnce_w_2 &
                                    + ( w_2 - wm )**4 )


    return

  end function calc_wp4_pdf

  !=============================================================================
  function calc_wp2xp_pdf( wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, &
                           varnce_w_2, varnce_x_1, varnce_x_2, &
                           corr_w_x_1, corr_w_x_2, mixt_frac ) &
  result( wp2xp )

    ! Description:
    ! Calculates <w'^2 x'> by integrating over the PDF of w and x.  The integral
    ! is:
    !
    ! <w'^2 x'>
    ! = INT(-inf:inf) INT(-inf:inf) ( w - <w> )^2 ( x - <x> ) P(w,x) dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, and
    ! P(w,x) is a two-component bivariate normal distribution of w and x.  The
    ! integrated equation is:
    !
    ! <w'^2 x'>
    ! = mixt_frac * ( ( mu_x_1 - <x> ) * ( ( mu_w_1 - <w> )^2 + sigma_w_1^2 )
    !                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
    !                   * ( mu_w_1 - <w> ) )
    !   + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )
    !                           * ( ( mu_w_2 - <w> )^2 + sigma_w_2^2 )
    !                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
    !                             * ( mu_w_2 - <w> ) );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    ! the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    ! standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    ! standard deviation of x in the 1st PDF component, sigma_x_2 is the
    ! standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    ! correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    ! correlation of w and x in the 2nd PDF component, and mixt_frac is the
    ! mixture fraction, which is the weight of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        two,   & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      wm,         & ! Mean of w (overall)                       [m/s]
      xm,         & ! Mean of x (overall)                       [units vary]
      w_1,        & ! Mean of w (1st PDF component)             [m/s]
      w_2,        & ! Mean of w (2nd PDF component)             [m/s]
      x_1,        & ! Mean of x (1st PDF component)             [units vary]
      x_2,        & ! Mean of x (2nd PDF component)             [units vary]
      varnce_w_1, & ! Variance of w (1st PDF component)         [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)         [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF comp.)    [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF comp.)    [-]
      mixt_frac     ! Mixture fraction                          [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      wp2xp    ! <w'^2 x'>                   [m^2/s^2 (units vary)]


    ! Calculate <w'^2 x'> by integrating over the PDF.
    wp2xp  = mixt_frac &
             * ( ( ( w_1 - wm )**2 + varnce_w_1 ) * ( x_1 - xm ) &
                 + two * corr_w_x_1 * sqrt( varnce_w_1 * varnce_x_1 ) &
                   * ( w_1 - wm ) ) &
             + ( one - mixt_frac ) &
               * ( ( ( w_2 - wm )**2 + varnce_w_2 ) * ( x_2 - xm ) &
                   + two * corr_w_x_2 * sqrt( varnce_w_2 * varnce_x_2 ) &
                     * ( w_2 - wm ) )


    return

  end function calc_wp2xp_pdf

  !=============================================================================
  function calc_wpxp2_pdf( wm, xm, w_1, w_2, x_1, x_2, varnce_w_1, &
                           varnce_w_2, varnce_x_1, varnce_x_2, &
                           corr_w_x_1, corr_w_x_2, mixt_frac ) &
  result( wpxp2 )

    ! Description:
    ! Calculates <w'x'^2> by integrating over the PDF of w and x.  The integral
    ! is:
    !
    ! <w'x'^2>
    ! = INT(-inf:inf) INT(-inf:inf) ( w - <w> ) ( x - <x> )^2 P(w,x) dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, and
    ! P(w,x) is a two-component bivariate normal distribution of w and x.  The
    ! integrated equation is:
    !
    ! <w'x'^2>
    ! = mixt_frac * ( ( mu_w_1 - <w> ) * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
    !                 + 2 * corr_w_x_1 * sigma_w_1 * sigma_x_1
    !                   * ( mu_x_1 - <x> ) )
    !   + ( 1 - mixt_frac ) * ( ( mu_w_2 - <w> )
    !                           * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 )
    !                           + 2 * corr_w_x_2 * sigma_w_2 * sigma_x_2
    !                             * ( mu_x_2 - <x> ) );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, sigma_w_1 is
    ! the standard deviation of w in the 1st PDF component, sigma_w_2 is the
    ! standard deviation of w in the 2nd PDF component, sigma_x_1 is the
    ! standard deviation of x in the 1st PDF component, sigma_x_2 is the
    ! standard deviation of x in the 2nd PDF component, corr_w_x_1 is the
    ! correlation of w and x in the 1st PDF component, corr_w_x_2 is the
    ! correlation of w and x in the 2nd PDF component, and mixt_frac is the
    ! mixture fraction, which is the weight of the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        two,   & ! Variable(s)
        one

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      wm,         & ! Mean of w (overall)                       [m/s]
      xm,         & ! Mean of x (overall)                       [units vary]
      w_1,        & ! Mean of w (1st PDF component)             [m/s]
      w_2,        & ! Mean of w (2nd PDF component)             [m/s]
      x_1,        & ! Mean of x (1st PDF component)             [units vary]
      x_2,        & ! Mean of x (2nd PDF component)             [units vary]
      varnce_w_1, & ! Variance of w (1st PDF component)         [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)         [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)         [(units vary)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)         [(units vary)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF comp.)    [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF comp.)    [-]
      mixt_frac     ! Mixture fraction                          [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      wpxp2    ! <w'x'^2>                   [m/s (units vary)^2]


    ! Calculate <w'x'^2> by integrating over the PDF.
    wpxp2 = mixt_frac &
            * ( ( w_1 - wm ) * ( ( x_1 - xm )**2 + varnce_x_1 ) &
                + two * corr_w_x_1 * sqrt( varnce_w_1 * varnce_x_1 ) &
                  * ( x_1 - xm ) ) &
            + ( one - mixt_frac ) &
              * ( ( w_2 - wm ) * ( ( x_2 - xm )**2 + varnce_x_2 ) &
                  + two * corr_w_x_2 * sqrt( varnce_w_2 * varnce_x_2 ) &
                    * ( x_2 - xm ) )


    return

  end function calc_wpxp2_pdf

  !=============================================================================
  function calc_wpxpyp_pdf( wm, xm, ym, w_1, w_2, x_1, x_2, &
                            y_1, y_2, varnce_w_1, varnce_w_2, &
                            varnce_x_1, varnce_x_2, varnce_y_1, &
                            varnce_y_2, corr_w_x_1, corr_w_x_2, &
                            corr_w_y_1, corr_w_y_2, corr_x_y_1, &
                            corr_x_y_2, mixt_frac ) &
  result( wpxpyp )

    ! Description:
    ! Calculates <w'x'y'> by integrating over the PDF of w, x, and y.  The
    ! integral is:
    !
    ! <w'x'y'>
    ! = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> ) ( x - <x> ) ( y - <y> ) P(w,x,y) dy dx dw;
    !
    ! where <w> is the overall mean of w, <x> is the overall mean of x, <y> is
    ! the overall mean of y, and P(w,x,y) is a two-component trivariate normal
    ! distribution of w, x, and y.  The integrated equation is:
    !
    ! <w'x'y'>
    ! = mixt_frac 
    !   * ( ( mu_w_1 - <w> ) * ( mu_x_1 - <x> ) * ( mu_y_1 - <y> )
    !       + corr_x_y_1 * sigma_x_1 * sigma_y_1 * ( mu_w_1 - <w> )
    !       + corr_w_y_1 * sigma_w_1 * sigma_y_1 * ( mu_x_1 - <x> )
    !       + corr_w_x_1 * sigma_w_1 * sigma_x_1 * ( mu_y_1 - <y> ) )
    !   + ( 1 - mixt_frac )
    !     * ( ( mu_w_2 - <w> ) * ( mu_x_2 - <x> ) * ( mu_y_2 - <y> )
    !         + corr_x_y_2 * sigma_x_2 * sigma_y_2 * ( mu_w_2 - <w> )
    !         + corr_w_y_2 * sigma_w_2 * sigma_y_2 * ( mu_x_2 - <x> )
    !         + corr_w_x_2 * sigma_w_2 * sigma_x_2 * ( mu_y_2 - <y> ) );
    !
    ! where mu_w_1 is the mean of w in the 1st PDF component, mu_w_2 is the mean
    ! of w in the 2nd PDF component, mu_x_1 is the mean of x in the 1st PDF
    ! component, mu_x_2 is the mean of x in the 2nd PDF component, mu_y_1 is the
    ! mean of y in the 1st PDF component, mu_y_2 is the mean of y in the 2nd PDF
    ! component, sigma_w_1 is the standard deviation of w in the 1st PDF
    ! component, sigma_w_2 is the standard deviation of w in the 2nd PDF
    ! component, sigma_x_1 is the standard deviation of x in the 1st PDF
    ! component, sigma_x_2 is the standard deviation of x in the 2nd PDF
    ! component, sigma_y_1 is the standard deviation of y in the 1st PDF
    ! component, sigma_y_2 is the standard deviation of y in the 2nd PDF
    ! component, corr_w_x_1 is the correlation of w and x in the 1st PDF
    ! component, corr_w_x_2 is the correlation of w and x in the 2nd PDF
    ! component, corr_w_y_1 is the correlation of w and y in the 1st PDF
    ! component, corr_w_y_2 is the correlation of w and y in the 2nd PDF
    ! component, corr_x_y_1 is the correlation of x and y in the 1st PDF
    ! component, corr_x_y_2 is the correlation of x and y in the 2nd PDF
    ! component, and mixt_frac is the mixture fraction, which is the weight of
    ! the 1st PDF component.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        one    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) ::  &
      wm,         & ! Mean of w (overall)                          [m/s]
      xm,         & ! Mean of x (overall)                          [x units]
      ym,         & ! Mean of y (overall)                          [y units]
      w_1,        & ! Mean of w (1st PDF component)                [m/s]
      w_2,        & ! Mean of w (2nd PDF component)                [m/s]
      x_1,        & ! Mean of x (1st PDF component)                [x units]
      x_2,        & ! Mean of x (2nd PDF component)                [x units]
      y_1,        & ! Mean of y (1st PDF component)                [y units]
      y_2,        & ! Mean of y (2nd PDF component)                [y units]
      varnce_w_1, & ! Variance of w (1st PDF component)            [m^2/s^2]
      varnce_w_2, & ! Variance of w (2nd PDF component)            [m^2/s^2]
      varnce_x_1, & ! Variance of x (1st PDF component)            [(x units)^2]
      varnce_x_2, & ! Variance of x (2nd PDF component)            [(x units)^2]
      varnce_y_1, & ! Variance of y (1st PDF component)            [(y units)^2]
      varnce_y_2, & ! Variance of y (2nd PDF component)            [(y units)^2]
      corr_w_x_1, & ! Correlation of w and x (1st PDF component)   [-]
      corr_w_x_2, & ! Correlation of w and x (2nd PDF component)   [-]
      corr_w_y_1, & ! Correlation of w and y (1st PDF component)   [-]
      corr_w_y_2, & ! Correlation of w and y (2nd PDF component)   [-]
      corr_x_y_1, & ! Correlation of x and y (1st PDF component)   [-]
      corr_x_y_2, & ! Correlation of x and y (2nd PDF component)   [-]
      mixt_frac     ! Mixture fraction                             [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) ::  & 
      wpxpyp    ! <w'x'y'>                   [m/s (units vary)]


    ! Calculate <w'x'y'> by integrating over the PDF.
    wpxpyp &
    = mixt_frac &
      * ( ( w_1 - wm ) * ( x_1 - xm ) * ( y_1 - ym ) &
          + corr_x_y_1 * sqrt( varnce_x_1 * varnce_y_1 ) * ( w_1 - wm ) &
          + corr_w_y_1 * sqrt( varnce_w_1 * varnce_y_1 ) * ( x_1 - xm ) &
          + corr_w_x_1 * sqrt( varnce_w_1 * varnce_x_1 ) * ( y_1 - ym ) ) &
      + ( one - mixt_frac ) &
        * ( ( w_2 - wm ) * ( x_2 - xm ) * ( y_2 - ym ) &
            + corr_x_y_2 * sqrt( varnce_x_2 * varnce_y_2 ) * ( w_2 - wm ) &
            + corr_w_y_2 * sqrt( varnce_w_2 * varnce_y_2 ) * ( x_2 - xm ) &
            + corr_w_x_2 * sqrt( varnce_w_2 * varnce_x_2 ) * ( y_2 - ym ) )


    return

  end function calc_wpxpyp_pdf

  !=============================================================================
  elemental subroutine calc_cloud_frac_component( mean_chi_i, stdev_chi_i, &
                                                  chi_at_sat, &
                                                  cloud_frac_i, rc_i )
    ! Description:
    ! Calculates the PDF component cloud water mixing ratio, rc_i, and cloud
    ! fraction, cloud_frac_i, for the ith PDF component.
    !
    ! The equation for cloud water mixing ratio, rc, at any point is:
    !
    ! rc = chi * H(chi);
    !
    ! and the equation for cloud fraction at a point, fc, is:
    !
    ! fc = H(chi);
    !
    ! where where extended liquid water mixing ratio, chi, is equal to cloud
    ! water mixing ratio, rc, when positive.  When the atmosphere is saturated
    ! at this point, cloud water is found, and rc = chi, while fc = 1.
    ! Otherwise, clear air is found at this point, and rc = fc = 0.
    !
    ! The mean of rc and fc is calculated by integrating over the PDF, such
    ! that:
    !
    ! <rc> = INT(-inf:inf) chi * H(chi) * P(chi) dchi; and
    !
    ! cloud_frac = <fc> = INT(-inf:inf) H(chi) * P(chi) dchi.
    !
    ! This can be rewritten as:
    !
    ! <rc> = INT(0:inf) chi * P(chi) dchi; and
    !
    ! cloud_frac = <fc> = INT(0:inf) P(chi) dchi;
    !
    ! and further rewritten as:
    !
    ! <rc> = SUM(i=1,N) mixt_frac_i INT(0:inf) chi * P_i(chi) dchi; and
    !
    ! cloud_frac = SUM(i=1,N) mixt_frac_i INT(0:inf) P_i(chi) dchi;
    !
    ! where N is the number of PDF components.  The equation for mean rc in the
    ! ith PDF component is:
    !
    ! rc_i = INT(0:inf) chi * P_i(chi) dchi;
    !
    ! and the equation for cloud fraction in the ith PDF component is:
    ! 
    ! cloud_frac_i = INT(0:inf) P_i(chi) dchi.
    !
    ! The component values are related to the overall values by:
    !
    ! <rc> = SUM(i=1,N) mixt_frac_i * rc_i; and
    !
    ! cloud_frac = SUM(i=1,N) mixt_frac_i * cloud_frac_i.

    ! References:
    !----------------------------------------------------------------------

    use constants_clubb, only: &
        chi_tol,        & ! Tolerance for pdf parameter chi       [kg/kg]
        sqrt_2pi,       & ! sqrt(2*pi)
        sqrt_2,         & ! sqrt(2)
        one,            & ! 1
        one_half,       & ! 1/2
        zero,           & ! 0
        max_num_stdevs, &
        eps

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    !----------- Input Variables -----------
    real( kind = core_rknd ), intent(in) :: &
      mean_chi_i,  & ! Mean of chi (old s) (ith PDF component)           [kg/kg]
      stdev_chi_i, & ! Standard deviation of chi (ith PDF component)     [kg/kg]
      chi_at_sat     ! Value of chi at saturation (0--liquid; neg.--ice) [kg/kg]

    !----------- Output Variables -----------
    real( kind = core_rknd ), intent(out) :: &
      cloud_frac_i, & ! Cloud fraction (ith PDF component)               [-]
      rc_i            ! Mean cloud water mixing ratio (ith PDF comp.)    [kg/kg]

    !----------- Local Variables -----------
    real( kind = core_rknd) :: zeta_i

    !----------- Begin Code -----------

    if ( ( abs( mean_chi_i - chi_at_sat ) <= eps .and. stdev_chi_i <= chi_tol ) &
           .or. ( mean_chi_i - chi_at_sat < - max_num_stdevs * stdev_chi_i ) ) then

        ! The mean of chi is at saturation and does not vary in the ith PDF component
        cloud_frac_i = zero
        rc_i         = zero

    elseif ( mean_chi_i - chi_at_sat > max_num_stdevs * stdev_chi_i ) then

        ! The mean of chi is multiple standard deviations above the saturation point.
        ! Thus, all cloud in the ith PDF component.
        cloud_frac_i = one
        rc_i         = mean_chi_i - chi_at_sat
    
    else

        ! The mean of chi is within max_num_stdevs of the saturation point.
        ! Thus, layer is partly cloudy, requires calculation.
        zeta_i = ( mean_chi_i - chi_at_sat ) / stdev_chi_i

        cloud_frac_i = one_half * ( one + erf( zeta_i / sqrt_2 )  )

        rc_i = ( mean_chi_i - chi_at_sat ) * cloud_frac_i &
               + stdev_chi_i * exp( - one_half * zeta_i**2 ) / ( sqrt_2pi )

    end if

    return
    
  end subroutine calc_cloud_frac_component

  !=============================================================================
  function calc_cloud_frac( cloud_frac_1, cloud_frac_2, mixt_frac ) &
  result( cloud_frac )

  ! Description:
  !   Given the the two pdf components of a cloud fraction, and the weight
  !   of the first component, this fuction calculates the cloud fraction,
  !   cloud_frac
  !
  ! References:
  !-----------------------------------------------------------------------
 
    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: & ! Constant(s)
        one,            & ! 1
        fstderr,        & ! Standard error output
        zero_threshold    ! A physical quantity equal to zero
    
    use clubb_precision, only: &
        core_rknd        ! Precision

    use error_code, only: &
        clubb_at_least_debug_level  ! Procedure
      
    implicit none
    
    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      cloud_frac_1, & ! First PDF component of cloud_frac               [-]
      cloud_frac_2, & ! Second PDF component of cloud_frac              [-]
      mixt_frac       ! Weight of 1st PDF component (Sk_w dependent)    [-]
    
    ! Output Variables
    real( kind = core_rknd), dimension(gr%nz) :: &
      cloud_frac    ! Cloud fraction    [-]

    ! Local Variable
    integer :: k    ! Vertical level loop index

    !----- Begin Code -----

    cloud_frac = mixt_frac * cloud_frac_1 + ( one - mixt_frac ) * cloud_frac_2
    
    ! Note: Brian added the following lines to ensure that there
    ! are never any negative liquid water values (or any negative
    ! cloud fraction values, for that matter).  According to
    ! Vince Larson, the analytic formula should not produce any
    ! negative results, but such computer-induced errors such as
    ! round-off error may produce such a value.  This has been
    ! corrected because Brian found a small negative value of
    ! rcm in the first timestep of the FIRE case.

    cloud_frac = max( zero_threshold, cloud_frac )

    if ( clubb_at_least_debug_level( 2 ) ) then
       do k = 1, gr%nz, 1
          if ( cloud_frac(k) > one ) then
             write(fstderr,*) "Cloud fraction > 1 at k = ", k
          endif
       enddo ! k = 1, gr%nz, 1
    endif

    cloud_frac = min( one, cloud_frac )


    return
    
  end function calc_cloud_frac

  !=============================================================================
  subroutine calc_xprcp_component( wm, rtm, thlm, um, vm, rcm,        & ! In
                                   w_i, rt_i, thl_i, u_i, v_i,        & ! In
                                   varnce_w_i,                        & ! In
                                   chi_i, stdev_chi_i, stdev_eta_i,   & ! In
                                   corr_w_chi_i, corr_chi_eta_i,      & ! In
                                   corr_u_w_i, corr_v_w_i,            & ! In
                                   crt_i, cthl_i, rc_i, cloud_frac_i, & ! In
                                   wprcp_contrib_comp_i,              & ! Out
                                   wp2rcp_contrib_comp_i,             & ! Out
                                   rtprcp_contrib_comp_i,             & ! Out
                                   thlprcp_contrib_comp_i,            & ! Out
                                   uprcp_contrib_comp_i,              & ! Out
                                   vprcp_contrib_comp_i               ) ! Out

    ! Description:
    ! Calculates the contribution to <w'rc'>, <w'^2 rc'>, <rt'rc'>, and
    ! <thl'rc'> from the ith PDF component.
    !
    !
    ! <w'rc'>
    ! -------
    !
    ! The value of <w'rc'> is calculated by integrating over the PDF:
    !
    ! <w'rc'>
    ! = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> ) ( rc - <rc> ) P(w,rt,thl) dthl drt dw;
    !
    ! where <w> is the overall mean of w, <rc> is the overall mean of rc, and
    ! P(w,rt,thl) is a two-component trivariate normal distribution of w, rt,
    ! and thl.  This equation is rewritten as:
    !
    ! <w'rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !     ( w - <w> ) ( rc - <rc> ) P_1(w,rt,thl) dthl drt dw
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !       ( w - <w> ) ( rc - <rc> ) P_2(w,rt,thl) dthl drt dw;
    !
    ! where mixt_frac is the mixture fraction, which is the weight of the 1st
    ! PDF component, and where P_1(w,rt,thl) and P_2(w,rt,thl) are the equations
    ! for the trivariate normal PDF of w, rt, and thl in the 1st and 2nd PDF
    ! components, respectively.  The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> ) ( rc - <rc> ) P_i(w,rt,thl) dthl drt dw;
    !
    ! where P_i(w,rt,thl) is the trivariate normal PDF of w, rt, and thl in the
    ! ith PDF component.  The PDF undergoes a PDF transformation in each PDF
    ! component, which is a change of variables and a translation, stretching,
    ! and rotation of the axes.  The PDF becomes a trivariate normal PDF that is
    ! written in terms of w, chi, and eta coordinates.  Cloud water mixing
    ! ratio, rc, is written in terms of extended liquid water mixing ratio, chi,
    ! such that:
    !
    ! rc = chi H(chi);
    !
    ! where H(chi) is the Heaviside step function.  The contribution from the
    ! ith PDF component to <w'rc'> can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> ) ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw;
    !
    ! where P_i(w,chi) is the bivariate normal PDF of w and chi in the ith PDF
    ! component.  The solved equation for the <w'rc'> contribution from the ith
    ! PDF component (wprcp_contrib_comp_i) is:
    !
    ! wprcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> ) ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw
    ! = ( mu_w_i - <w> )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + corr_w_chi_i * sigma_w_i * sigma_chi_i
    !     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
    !
    ! where mu_w_i is the mean of w in the ith PDF component, mu_chi_i is the
    ! mean of chi in the ith PDF component, sigma_w_i is the standard deviation
    ! of w in the ith PDF component, sigma_chi_i is the standard deviation of
    ! chi in the ith PDF component, and corr_w_chi_i is the correlation of w and
    ! chi in the ith PDF component.
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! wprcp_contrib_comp_i
    ! = | ( mu_w_i - <w> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
    !   | ( mu_w_i - <w> ) * ( -<rc> ); when mu_chi_i <= 0.
    !
    !
    ! <w'^2 rc'>
    ! ----------
    !
    ! The value of <w'^2 rc'> is calculated by integrating over the PDF:
    !
    ! <w'^2 rc'>
    ! = INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> )^2 ( rc - <rc> ) P(w,rt,thl) dthl drt dw.
    !
    ! This equation is rewritten as:
    !
    ! <w'^2 rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !     ( w - <w> )^2 ( rc - <rc> ) P_1(w,rt,thl) dthl drt dw
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    !       ( w - <w> )^2 ( rc - <rc> ) P_2(w,rt,thl) dthl drt dw.
    !
    ! The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> )^2 ( rc - <rc> ) P_i(w,rt,thl) dthl drt dw.
    !
    ! The PDF undergoes a PDF transformation in each PDF component, and becomes
    ! a trivariate normal PDF that is written in terms of w, chi, and eta
    ! coordinates.  The contribution from the ith PDF component to <w'^2 rc'>
    ! can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( w - <w> )^2 ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw.
    !
    ! The solved equation for the <w'^2 rc'> contribution from the ith PDF
    ! component (wp2rcp_contrib_comp_i) is:
    !
    ! wp2rcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( w - <w> )^2 ( chi H(chi) - <rc> ) P_i(w,chi) dchi dw
    ! = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + ( mu_w_i - <w> ) * corr_w_chi_i * sigma_w_i * sigma_chi_i
    !     * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !   + 1/sqrt(2*pi) * corr_w_chi_i^2 * sigma_w_i^2 * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }.
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! wp2rcp_contrib_comp_i
    ! = | ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( mu_chi_i - <rc> );
    !   |     when mu_chi_i > 0;
    !   | ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( -<rc> );
    !   |     when mu_chi_i <= 0.
    !
    !
    ! <rt'rc'>
    ! --------
    !
    ! The value of <rt'rc'> is calculated by integrating over the PDF:
    !
    ! <rt'rc'>
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( rt - <rt> ) ( rc - <rc> ) P(rt,thl) dthl drt;
    !
    ! where <rt> is the overall mean of rt, and where P(rt,thl) is a
    ! two-component bivariate normal distribution of rt and thl.  This equation
    ! is rewritten as:
    !
    ! <rt'rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf)
    !     ( rt - <rt> ) ( rc - <rc> ) P_1(rt,thl) dthl drt
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf)
    !       ( rt - <rt> ) ( rc - <rc> ) P_2(rt,thl) dthl drt;
    !
    ! where P_1(rt,thl) and P_2(rt,thl) are the equations for the bivariate
    ! normal PDF of rt and thl in the 1st and 2nd PDF components, respectively.
    ! The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( rt - <rt> ) ( rc - <rc> ) P_i(rt,thl) dthl drt;
    !
    ! where P_i(rt,thl) is the bivariate normal PDF of rt and thl in the ith PDF
    ! component.  The PDF undergoes a PDF transformation in each PDF component,
    ! and becomes a bivariate normal PDF that is written in terms of chi and
    ! eta coordinates.  Total water mixing ratio, rt, is rewritten in terms of
    ! chi and eta by:
    !
    ! rt = mu_rt_i
    !      + ( ( eta - mu_eta_i ) + ( chi - mu_chi_i ) ) / ( 2 * crt_i );
    !
    ! where mu_rt_i is the mean of rt in the ith PDF component, mu_eta_i is the
    ! mean of eta in the ith PDF component, and crt_i is a coefficient on rt in
    ! the chi/eta transformation equations.  The contribution from the ith PDF
    ! component to <rt'rc'> can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( mu_rt_i - <rt> + ( eta - mu_eta_i ) / ( 2 * crt_i )
    !   + ( chi - mu_chi_i ) / ( 2 * crt_i ) )
    ! * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi;
    !
    ! where P_i(chi,eta) is the bivariate normal PDF of chi and eta in the ith
    ! PDF component.  The solved equation for the <rt'rc'> contribution from the
    ! ith PDF component (rtprcp_contrib_comp_i) is:
    !
    ! rtprcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( mu_rt_i - <rt> + ( eta - mu_eta_i ) / ( 2 * crt_i )
    !     + ( chi - mu_chi_i ) / ( 2 * crt_i ) )
    !   * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi
    ! = ( mu_rt_i - <rt> )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i + sigma_chi_i ) / ( 2 * crt_i )
    !     * sigma_chi_i
    !     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) );
    !
    ! where sigma_eta_i is the standard deviation of eta in the ith PDF
    ! component and corr_chi_eta_i is the correlation of chi and eta in the ith
    ! PDF component.
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! rtprcp_contrib_comp_i
    ! = | ( mu_rt_i - <rt> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
    !   | ( mu_rt_i - <rt> ) * ( -<rc> ); when mu_chi_i <= 0.
    !
    !
    ! <thl'rc'>
    ! ---------
    !
    ! The value of <thl'rc'> is calculated by integrating over the PDF:
    !
    ! <thl'rc'>
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( thl - <thl> ) ( rc - <rc> ) P(rt,thl) dthl drt;
    !
    ! where <thl> is the overall mean of thl.  This equation is rewritten as:
    !
    ! <thl'rc'>
    ! = mixt_frac 
    !   * INT(-inf:inf) INT(-inf:inf)
    !     ( thl - <thl> ) ( rc - <rc> ) P_1(rt,thl) dthl drt
    !   + ( 1 - mixt_frac )
    !     * INT(-inf:inf) INT(-inf:inf)
    !       ( thl - <thl> ) ( rc - <rc> ) P_2(rt,thl) dthl drt.
    !
    ! The contribution from the ith PDF component is:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( thl - <thl> ) ( rc - <rc> ) P_i(rt,thl) dthl drt.
    !
    ! The PDF undergoes a PDF transformation in each PDF component, and becomes
    ! a bivariate normal PDF that is written in terms of chi and eta
    ! coordinates.  Liquid water potential temperature, thl, is rewritten in
    ! terms of chi and eta by:
    !
    ! thl = mu_thl_i
    !       + ( ( eta - mu_eta_i ) - ( chi - mu_chi_i ) ) / ( 2 * cthl_i );
    !
    ! where mu_thl_i is the mean of thl in the ith PDF component and cthl_i is a
    ! coefficient on thl in the chi/eta transformation equations.  The
    ! contribution from the ith PDF component to <thl'rc'> can be written as:
    !
    ! INT(-inf:inf) INT(-inf:inf)
    ! ( mu_thl_i - <thl> + ( eta - mu_eta_i ) / ( 2 * cthl_i )
    !   - ( chi - mu_chi_i ) / ( 2 * cthl_i ) )
    ! * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi.
    !
    ! The solved equation for the <thl'rc'> contribution from the ith PDF
    ! component (thlprcp_contrib_comp_i) is:
    !
    ! thlprcp_contrib_comp_i
    ! = INT(-inf:inf) INT(-inf:inf)
    !   ( mu_thl_i - <thl> + ( eta - mu_eta_i ) / ( 2 * cthl_i )
    !     - ( chi - mu_chi_i ) / ( 2 * cthl_i ) )
    !   * ( chi H(chi) - <rc> ) P_i(chi,eta) deta dchi
    ! = ( mu_thl_i - <thl> )
    !   * ( mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !       + 1/sqrt(2*pi) * sigma_chi_i
    !         * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) } - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i - sigma_chi_i ) / ( 2 * cthl_i )
    !     * sigma_chi_i
    !     * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
    !
    ! Special case:  sigma_chi_i = 0.
    !
    ! In the special case that sigma_chi_i = 0, chi, as well as rc, are constant
    ! in the ith PDF component.  The equation becomes:
    !
    ! thlprcp_contrib_comp_i
    ! = | ( mu_thl_i - <thl> ) * ( mu_chi_i - <rc> ); when mu_chi_i > 0;
    !   | ( mu_thl_i - <thl> ) * ( -<rc> ); when mu_chi_i <= 0.
    !
    !
    ! Use equations for PDF component cloud fraction cloud water mixing ratio
    ! -----------------------------------------------------------------------
    !
    ! The equation for cloud fraction in the ith PDF component, fc_i, is:
    !
    ! fc_i = 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) ).
    !
    ! In the special case that sigma_chi_i = 0, the equation becomes:
    !
    ! fc_i = | 1; when mu_chi_i > 0;
    !        | 0; when mu_chi_i <= 0.
    !
    ! The equation for mean cloud water mixing ratio in the ith PDF component,
    ! rc_i, is:
    !
    ! rc_i
    ! = mu_chi_i * 1/2 * ( 1 + erf( mu_chi_i / ( sqrt(2) * sigma_chi_i ) ) )
    !   + 1/sqrt(2*pi) * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }
    ! = mu_chi_i * fc_i
    !   + 1/sqrt(2*pi) * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) }.
    !
    ! In the special case that sigma_chi_i = 0, the equation becomes:
    !
    ! rc_i = | mu_chi_i; when mu_chi_i > 0;
    !        | 0; when mu_chi_i <= 0.
    !
    ! The above equations can be substituted into the equations for
    ! wprcp_contrib_comp_i, wp2rcp_contrib_comp_i, rtprcp_contrib_comp_i, and
    ! thlprcp_contrib_comp_i.  The new equations are:
    !
    ! wprcp_contrib_comp_i
    ! = ( mu_w_i - <w> ) * ( rc_i - <rc> )
    !   + corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i;
    !
    ! wp2rcp_contrib_comp_i
    ! = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( rc_i - <rc> )
    !   + 2 * ( mu_w_i - <w> ) * corr_w_chi_i * sigma_w_i * sigma_chi_i * fc_i
    !   + 1/sqrt(2*pi) * corr_w_chi_i^2 * sigma_w_i^2 * sigma_chi_i
    !     * exp{ - mu_chi_i^2 / ( 2 * sigma_chi_i^2 ) };
    !
    ! rtprcp_contrib_comp_i
    ! = ( mu_rt_i - <rt> ) * ( rc_i - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i + sigma_chi_i ) / ( 2 * crt_i )
    !     * sigma_chi_i * fc_i; and
    !
    ! thlprcp_contrib_comp_i
    ! = ( mu_thl_i - <thl> ) * ( rc_i - <rc> )
    !   + ( corr_chi_eta_i * sigma_eta_i - sigma_chi_i ) / ( 2 * cthl_i )
    !     * sigma_chi_i * fc_i.
    !
    ! While the above equations reduce to their listed versions in the special
    ! case that sigma_chi_i = 0, those versions are faster to calculate.  When
    ! mu_chi_i > 0, they are:
    !
    ! wprcp_contrib_comp_i = ( mu_w_i - <w> ) * ( mu_chi_i - <rc> );
    ! wp2rcp_contrib_comp_i
    ! = ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * ( mu_chi_i - <rc> );
    ! rtprcp_contrib_comp_i = ( mu_rt_i - <rt> ) * ( mu_chi_i - <rc> ); and
    ! thlprcp_contrib_comp_i = ( mu_thl_i - <thl> ) * ( mu_chi_i - <rc> );
    !
    ! and when mu_chi_i <= 0, they are:
    !
    ! wprcp_contrib_comp_i = - ( mu_w_i - <w> ) * <rc>;
    ! wp2rcp_contrib_comp_i = - ( ( mu_w_i - <w> )^2 + sigma_w_i^2 ) * <rc>;
    ! rtprcp_contrib_comp_i = - ( mu_rt_i - <rt> ) * <rc>; and
    ! thlprcp_contrib_comp_i = - ( mu_thl_i - <thl> ) * <rc>.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        sqrt_2pi,       & ! Variable(s)
        two,            &
        zero,           &
        chi_tol,        &
        cloud_frac_min

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wm,             & ! Mean of w (overall)                          [m/s]
      rtm,            & ! Mean of rt (overall)                         [kg/kg]
      thlm,           & ! Mean of thl (overall)                        [K]
      um,             & ! Mean of eastward wind (overall)              [m/s]
      vm,             & ! Mean of northward wind (overall)             [m/s]
      rcm,            & ! Mean of rc (overall)                         [kg/kg]
      w_i,            & ! Mean of w (ith PDF component)                [m/s]
      rt_i,           & ! Mean of rt (ith PDF component)               [kg/kg]
      thl_i,          & ! Mean of thl (ith PDF component)              [K]
      u_i,            & ! Mean of eastward wind (ith PDF component)    [m/s]
      v_i,            & ! Mean of northward wind (ith PDF component)   [m/s]
      varnce_w_i,     & ! Variance of w (ith PDF component)            [m^2/s^2]
      chi_i,          & ! Mean of chi (ith PDF component)              [kg/kg]
      stdev_chi_i,    & ! Standard deviation of chi (ith PDF comp.)    [kg/kg]
      stdev_eta_i,    & ! Standard deviation of eta (ith PDF comp.)    [kg/kg]
      corr_w_chi_i,   & ! Correlation of w and chi (ith PDF component) [-]
      corr_chi_eta_i, & ! Correlation of chi and eta (ith PDF comp.)   [-]
      corr_u_w_i,     & ! Correlation of u and w (ith PDF component)   [-]
      corr_v_w_i,     & ! Correlation of v and w (ith PDF component)   [-]
      crt_i,          & ! Coef. on rt in chi/eta eqns. (ith PDF comp.) [-]
      cthl_i,         & ! Coef. on thl: chi/eta eqns. (ith PDF comp.)  [kg/kg/K]
      rc_i,           & ! Mean of rc (ith PDF component)               [kg/kg]
      cloud_frac_i      ! Cloud fraction (ith PDF component)           [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      wprcp_contrib_comp_i,   & ! <w'rc'> contrib. (ith PDF comp.)  [m/s(kg/kg)]
      wp2rcp_contrib_comp_i,  & ! <w'^2rc'> contrib. (ith comp) [m^2/s^2(kg/kg)]
      rtprcp_contrib_comp_i,  & ! <rt'rc'> contrib. (ith PDF comp.)  [kg^2/kg^2]
      thlprcp_contrib_comp_i, & ! <thl'rc'> contrib. (ith PDF comp.)  [K(kg/kg)]
      uprcp_contrib_comp_i,   & ! <u'rc'> contrib. (ith PDF comp.)  [m/s(kg/kg)]
      vprcp_contrib_comp_i      ! <v'rc'> contrib. (ith PDF comp.)  [m/s(kg/kg)]

    ! ---------------------- Begin Code ------------------
    
    ! Changing these conditionals may result in inconsistencies with the conditional
    ! statements located in calc_cloud_frac_component

    wprcp_contrib_comp_i = ( w_i - wm ) * ( rc_i - rcm )

    wp2rcp_contrib_comp_i = ( ( w_i - wm )**2 + varnce_w_i ) * ( rc_i - rcm )

    rtprcp_contrib_comp_i = ( rt_i - rtm ) * ( rc_i - rcm ) &
                            + ( corr_chi_eta_i * stdev_eta_i + stdev_chi_i ) &
                              / ( two * crt_i ) * stdev_chi_i * cloud_frac_i

    thlprcp_contrib_comp_i = ( thl_i - thlm ) * ( rc_i - rcm ) &
                             + ( corr_chi_eta_i * stdev_eta_i - stdev_chi_i ) &
                               / ( two * cthl_i ) * stdev_chi_i * cloud_frac_i

    uprcp_contrib_comp_i = ( u_i - um ) * ( rc_i - rcm )

    vprcp_contrib_comp_i = ( v_i - vm ) * ( rc_i - rcm )

    ! If iiPDF_type isn't iiPDF_ADG1 or iiPDF_ADG2, so corr_w_chi_i /= 0
    !   (and perhaps corr_u_w_i /= 0).
    if ( .not. ( iiPDF_type == iiPDF_ADG1 .or. iiPDF_type == iiPDF_ADG2 ) ) then

        ! Chi varies significantly in the ith PDF component (stdev_chi > chi_tol)
        ! and there is some cloud (0 < cloud_frac <= 1)
        where ( stdev_chi_i > chi_tol .and. cloud_frac_i > zero )

            wprcp_contrib_comp_i = wprcp_contrib_comp_i &
                                   + corr_w_chi_i * sqrt( varnce_w_i ) * stdev_chi_i * cloud_frac_i

            wp2rcp_contrib_comp_i = wp2rcp_contrib_comp_i &
                                    + two * ( w_i - wm ) * corr_w_chi_i &
                                      * sqrt( varnce_w_i ) * stdev_chi_i * cloud_frac_i &
                                    + corr_w_chi_i**2 * varnce_w_i * stdev_chi_i &
                                      * exp( - chi_i**2 / ( two * stdev_chi_i**2 ) ) / sqrt_2pi

            ! In principle, uprcp_contrib_comp_i might depend on corr_u_w_i here.

        end where

    end if 


    return

  end subroutine calc_xprcp_component

  !=============================================================================
  subroutine calc_vert_avg_cf_component &
                  ( nz, k, z_vals, chi, stdev_chi, chi_at_sat, &
                    cloud_frac_i, rc_i )
  ! Description:
  !   This subroutine is similar to calc_cloud_frac_component, but
  !   resolves cloud_frac and rc at an arbitrary number of vertical levels
  !   in the vicinity of the desired level. This may give a better
  !   parameterization of sub-grid atmospheric conditions.
  !
  ! References:
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd

    implicit none

    intrinsic :: sum

    ! Local Constants
    integer, parameter :: &
      n_points = 9       ! Number of vertical levels to use in averaging
                         ! (arbitrary, but must be odd)

    ! Input Variables
    integer, intent(in) :: &
      nz, &       ! Number of vertical levels                         [count]
      k           ! Level at which cloud_frac is to be computed       [count]

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals,    & ! Height at each vertical level                   [m]
      chi,       & ! Value of chi (old s)                            [kg/kg]
      stdev_chi, & ! Standard deviation of chi                       [kg/kg]
      chi_at_sat   ! Value of chi at saturation with respect to ice  [kg/kg]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      cloud_frac_i, & ! Vertically averaged cloud fraction               [-]
      rc_i            ! Vertically averaged cloud water mixing ratio     [kg/kg]

    ! Local Variables
    real( kind = core_rknd ), dimension(n_points) :: &
      chi_ref,           &   ! chi (old s) evaluated on refined grid     [kg/kg]
      stdev_chi_ref,     &   ! stdev_chi evaluated on refined grid       [kg/kg]
      cloud_frac_ref,    &   ! cloud_frac evaluated on refined grid      [-]
      rc_ref                 ! r_c evaluated on refined grid             [kg/kg]
      
  !-----------------------------------------------------------------------

    !----- Begin Code -----
    chi_ref = interp_var_array( n_points, nz, k, z_vals, chi )
    stdev_chi_ref = interp_var_array( n_points, nz, k, z_vals, stdev_chi )
    ! We could optionally compute chi_at_sat in an analogous manner. For now,
    ! use chi_at_sat(k) as an approximation.

    ! Compute cloud_frac and r_c at each refined grid level
    call calc_cloud_frac_component( chi_ref(:), stdev_chi_ref(:), chi_at_sat(k), & ! Intent(in)
                                    cloud_frac_ref(:), rc_ref(:) )                  ! Intent(out)

    cloud_frac_i = sum( cloud_frac_ref(:) ) / real( n_points, kind=core_rknd )
    rc_i = sum( rc_ref(:) ) / real( n_points, kind=core_rknd )

    return
  end subroutine calc_vert_avg_cf_component

  !=============================================================================
  function interp_var_array( n_points, nz, k, z_vals, var )

  ! Description:
  !   Interpolates a variable to an array of values about a given level

  ! References
  !-----------------------------------------------------------------------

  use clubb_precision, only: &
    core_rknd           ! Constant

  implicit none

  ! Input Variables
  integer, intent(in) :: &
    n_points, & ! Number of points to interpolate to (must be odd and >= 3)
    nz,       & ! Total number of vertical levels
    k           ! Center of interpolation array

  real( kind = core_rknd ), dimension(nz), intent(in) :: &
    z_vals, &         ! Height at each vertical level           [m]
    var               ! Variable values on grid                 [units vary]

  ! Output Variables
  real( kind = core_rknd ), dimension(n_points) :: &
    interp_var_array  ! Interpolated values of variable         [units vary]

  ! Local Variables
  real( kind = core_rknd ) :: &
    dz    ! Distance between vertical levels

  real( kind = core_rknd ) :: &
    z_val             ! Height at some sub-grid level

  integer :: &
    i, &                      ! Loop iterator

    subgrid_lev_count         ! Number of refined grid points located between
                              ! two defined grid levels

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    ! Place a point at each of k-1, k, and k+1.
    interp_var_array(1) = var_value_integer_height( nz, k-1, z_vals, var )
    interp_var_array((n_points+1)/2) = var_value_integer_height( nz, k, z_vals, var )
    interp_var_array(n_points) = var_value_integer_height( nz, k+1, z_vals, var )

    subgrid_lev_count = (n_points - 3) / 2

    ! Lower half
    if ( k == 1 ) then
      dz = (z_vals(2) - z_vals(1)) / real( subgrid_lev_count+1, kind=core_rknd )
    else
      dz = (z_vals(k) - z_vals(k-1)) / real( subgrid_lev_count+1, kind=core_rknd )
    end if
    do i=1, subgrid_lev_count
      z_val = z_vals(k) - real( i, kind=core_rknd ) * dz
      interp_var_array(1+i) &
      = var_subgrid_interp( nz, k, z_vals, var, z_val, l_below=.true. )
    end do

    ! Upper half
    if ( k == nz ) then
      dz = ( z_vals(nz) - z_vals(nz-1) ) / real( subgrid_lev_count+1, kind=core_rknd )
    else
      dz = ( z_vals(k+1) - z_vals(k) ) / real( subgrid_lev_count+1, kind=core_rknd )
    end if
    do i=1, (n_points-3)/2
      z_val = z_vals(k) + real( i, kind=core_rknd ) * dz
      interp_var_array((n_points+1)/2+i) &
      = var_subgrid_interp( nz, k, z_vals, var, z_val, l_below=.false. )
    end do

    return
  end function interp_var_array

  !=============================================================================
  function var_value_integer_height( nz, k, z_vals, var_grid_value ) result( var_value )

  ! Description
  !   Returns the value of a variable at an integer height between 0 and
  !   nz+1 inclusive, using extrapolation when k==0 or k==nz+1

  ! References
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd       ! Constant

    use interpolation, only: &
      mono_cubic_interp  ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz,    & ! Total number of vertical levels
      k        ! Level to resolve variable value

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals,            & ! Height at each vertical level                  [m]
      var_grid_value       ! Value of variable at each grid level           [units vary]

    ! Output Variables
    real( kind = core_rknd ) :: &
      var_value            ! Value of variable at height level              [units vary]

    ! Local Variables
    integer :: km1, k00, kp1, kp2
  !-----------------------------------------------------------------------

    !----- Begin Code -----

    if ( k >= 1 .and. k <= nz ) then
      ! This is the simple case. No extrapolation necessary.
      var_value = var_grid_value(k)
    else if ( k == 0 ) then
      ! Extrapolate below the lower boundary
      km1 = nz
      k00 = 1
      kp1 = 2
      kp2 = 3
      var_value = mono_cubic_interp( z_vals(1)-(z_vals(2)-z_vals(1)), &
                                     km1, k00, kp1, kp2, &
                                     z_vals(km1), z_vals(k00), z_vals(kp1), z_vals(kp2), &
                                     var_grid_value(km1), var_grid_value(k00), &
                                     var_grid_value(kp1), var_grid_value(kp2) )
    else if ( k == nz+1 ) then
      ! Extrapolate above the upper boundary
      km1 = nz
      k00 = nz-1
      kp1 = nz
      kp2 = nz
      var_value = mono_cubic_interp( z_vals(nz)+(z_vals(nz)-z_vals(nz-1)), &
                                     km1, k00, kp1, kp2, &
                                     z_vals(km1), z_vals(k00), z_vals(kp1), z_vals(kp2), &
                                     var_grid_value(km1), var_grid_value(k00), &
                                     var_grid_value(kp1), var_grid_value(kp2) )
    else
      ! Invalid height requested
      var_value = -999._core_rknd
    end if ! k > 1 .and. k < nz
    return
  end function var_value_integer_height

  !=============================================================================
  function var_subgrid_interp( nz, k, z_vals, var, z_interp, l_below ) result( var_value )

  ! Description
  !   Interpolates (or extrapolates) a variable to a value between grid
  !   levels

  ! References
  !-----------------------------------------------------------------------

    use clubb_precision, only: &
      core_rknd       ! Constant

    use interpolation, only: &
      mono_cubic_interp   ! Procedure

    implicit none

    ! Input Variables
    integer, intent(in) :: &
      nz, &         ! Number of vertical levels
      k             ! Grid level near interpolation target

    real( kind = core_rknd ), dimension(nz), intent(in) :: &
      z_vals, &     ! Height at each grid level          [m]
      var           ! Variable values at grid levels     [units vary]

    real( kind = core_rknd ), intent(in) :: &
      z_interp      ! Interpolation target height        [m]

    logical, intent(in) :: &
      l_below       ! True if z_interp < z_vals(k), false otherwise

    ! Output Variable
    real( kind = core_rknd ) :: &
      var_value     ! Interpolated value of variable     [units vary]

    ! Local Variables
    integer :: km1, k00, kp1, kp2 ! Parameters for call to mono_cubic_interp
  !----------------------------------------------------------------------

    !----- Begin Code -----
    if ( l_below ) then

      if ( k == 1 ) then ! Extrapolation
        km1 = nz
        k00 = 1
        kp1 = 2
        kp2 = 3
      else if ( k == 2 ) then
        km1 = 1
        k00 = 1
        kp1 = 2
        kp2 = 3
      else if ( k == nz ) then
        km1 = nz-2
        k00 = nz-1
        kp1 = nz
        kp2 = nz
      else
        km1 = k-2
        k00 = k-1
        kp1 = k
        kp2 = k+1
      end if ! k == 1

    else ! .not. l_below

      if ( k == 1 ) then
        km1 = 1
        k00 = 1
        kp1 = 2
        kp2 = 3
      else if ( k == nz-1 ) then
        km1 = nz-2
        k00 = nz-1
        kp1 = nz
        kp2 = nz
      else if ( k == nz ) then ! Extrapolation
        km1 = nz
        k00 = nz-1
        kp1 = nz
        kp2 = nz
      else
        km1 = k-1
        k00 = k
        kp1 = k+1
        kp2 = k+2
      end if ! k == 1

    end if ! l_below

    ! Now perform the interpolation
    var_value = mono_cubic_interp( z_interp, km1, k00, kp1, kp2, &
                                   z_vals(km1), z_vals(k00), z_vals(kp1), z_vals(kp2), &
                                   var(km1), var(k00), var(kp1), var(kp2) )

    return

  end function var_subgrid_interp

  !=============================================================================

end module pdf_closure_module
