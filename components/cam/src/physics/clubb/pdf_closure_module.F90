!---------------------------------------------------------------------------
! $Id: pdf_closure_module.F90 8106 2016-05-17 23:29:01Z raut@uwm.edu $
!===============================================================================
module pdf_closure_module

  implicit none

  public :: pdf_closure, calc_vert_avg_cf_component

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
  subroutine pdf_closure( hydromet_dim, p_in_Pa, exner, thv_ds, wm, &
                          wp2, wp3, sigma_sqd_w,                    &
                          Skw, Skthl, Skrt, rtm, rtp2,              &
                          wprtp, thlm, thlp2,                       &
                          wpthlp, rtpthlp, sclrm,                   &
                          wpsclrp, sclrp2, sclrprtp,                &
                          sclrpthlp, level,                         &
#ifdef GFDL
                          RH_crit,  do_liquid_only_in_clubb,        & ! h1g, 2010-06-15
#endif
                          wphydrometp, wp2hmp,                      &
                          rtphmp, thlphmp,                          &
                          wp4, wprtp2, wp2rtp,                      &
                          wpthlp2, wp2thlp, wprtpthlp,              &
                          cloud_frac, ice_supersat_frac,            &
                          rcm, wpthvp, wp2thvp, rtpthvp,            &
                          thlpthvp, wprcp, wp2rcp, rtprcp,          &
                          thlprcp, rcp2, pdf_params,                &
                          err_code,                                 &
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
    !   Eqn. 29, 30, 31, 32 & 33  on p. 3547 of
    !   ``A PDF-Based Model for Boundary Layer Clouds. Part I:
    !   Method and Model Description'' Golaz, et al. (2002)
    !   JAS, Vol. 59, pp. 3540--3551.
    !----------------------------------------------------------------------

    use constants_clubb, only: &  ! Constants
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
        w_tol_sqd,      & ! Tolerance for w'^2                  [m^2/s^2]
        rt_tol,         & ! Tolerance for r_t                   [kg/kg]
        thl_tol,        & ! Tolerance for th_l                  [K]
        T_freeze_K,     & ! Freezing point of water             [K]
        fstderr,        &
        zero_threshold, &
        chi_tol, &
        eps, &
        w_tol


    use parameters_model, only: &
        sclr_tol,          & ! Array of passive scalar tolerances  [units vary]
        sclr_dim,          & ! Number of passive scalar variables
        mixt_frac_max_mag    ! Maximum values for PDF parameter 'mixt_frac'

    use parameters_tunable, only: & 
        beta, &  ! Variable(s)
        Skw_denom_coef

    use pdf_parameter_module, only:  &
        pdf_parameter  ! type

    use array_index, only: &
        l_mix_rat_hm  ! Variable(s)

    use anl_erf, only:  & 
        erf ! Procedure(s)
    ! The error function

    use numerical_check, only:  & 
        pdf_closure_check ! Procedure(s)

    use saturation, only:  & 
        sat_mixrat_liq, & ! Procedure(s)
        sat_mixrat_ice

    use error_code, only: & 
        clubb_no_error ! Constant(s)

    use error_code, only:  & 
        clubb_at_least_debug_level, & ! Procedure(s)
        fatal_error

    use stats_variables, only: &
        iwp4,       & ! Variables
        ircp2,      &
        iwprtp2,    &
        iwprtpthlp, &
        iwpthlp2

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    use model_flags, only:&
        l_use_ADG2, &
        l_use_3D_closure

    implicit none

    intrinsic :: sqrt, exp, min, max, abs, present

    ! Input Variables
    integer, intent(in) :: &
      hydromet_dim   ! Number of hydrometeor species              [#]

    real( kind = core_rknd ), intent(in) ::  & 
      p_in_Pa,     & ! Pressure                                   [Pa]
      exner,       & ! Exner function                             [-]
      thv_ds,      & ! Dry, base-state theta_v (ref. th_l here)   [K]
      wm,          & ! mean w-wind component (vertical velocity)  [m/s] 
      wp2,         & ! w'^2                                       [m^2/s^2] 
      wp3,         & ! w'^3                                       [m^3/s^3]
      Skw,         & ! Skewness of w                              [-]
      Skthl,       & ! Skewness of thl                            [-]
      Skrt,        & ! Skewness of rt                             [-]
      rtm,         & ! Mean total water mixing ratio              [kg/kg]
      rtp2,        & ! r_t'^2                                     [(kg/kg)^2]
      wprtp,       & ! w'r_t'                                     [(kg/kg)(m/s)]
      thlm,        & ! Mean liquid water potential temperature    [K]
      thlp2,       & ! th_l'^2                                    [K^2]
      wpthlp,      & ! w'th_l'                                    [K(m/s)]
      rtpthlp        ! r_t'th_l'                                  [K(kg/kg)]

    real( kind = core_rknd ), dimension(sclr_dim), intent(in) ::  & 
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

    integer, intent(in) ::  &
      level  ! Thermodynamic level for which calculations are taking place.

    real( kind = core_rknd ), dimension(hydromet_dim), intent(in) :: &
      wphydrometp, & ! Covariance of w and a hydrometeor    [(m/s) <hm units>]
      wp2hmp,      & ! Third-order moment:  < w'^2 hm' >    [(m/s)^2 <hm units>]
      rtphmp,      & ! Covariance of rt and a hydrometeor   [(kg/kg) <hm units>]
      thlphmp        ! Covariance of thl and a hydrometeor  [K <hm units>]

    real( kind = core_rknd ), intent(inout) :: &
      ! If l_use_ADG2, this gets overwritten. Therefore, intent(inout).
      ! otherwise it should be intent(in)
      sigma_sqd_w   ! Width of individual w plumes               [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) ::  & 
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

    type(pdf_parameter), intent(out) :: & 
      pdf_params     ! pdf paramters         [units vary]

    integer, intent(out) :: & 
      err_code       ! Are the outputs usable numbers?

    ! Output (passive scalar variables)

    real( kind = core_rknd ), intent(out), dimension(sclr_dim) ::  & 
      sclrpthvp, & 
      sclrprcp, & 
      wpsclrp2, & 
      wpsclrprtp, & 
      wpsclrpthlp, & 
      wp2sclrp

    ! Local Variables

    real( kind = core_rknd ) ::  & 
      w_1_n, w_2_n,              &
      thl_1_n, thl_2_n,          &
      rt_1_n, rt_2_n

    ! Variables that are stored in derived data type pdf_params.
    real( kind = core_rknd ) ::  &
      w_1,          & ! Mean of w (1st PDF component)                       [m/s]
      w_2,          & ! Mean of w (2nd PDF component)                       [m/s]
      varnce_w_1,   & ! Variance of w (1st PDF component)               [m^2/s^2]
      varnce_w_2,   & ! Variance of w (2nd PDF component)               [m^2/s^2]
      rt_1,         & ! Mean of r_t (1st PDF component)                   [kg/kg]
      rt_2,         & ! Mean of r_t (2nd PDF component)                   [kg/kg]
      varnce_rt_1,  & ! Variance of r_t (1st PDF component)           [kg^2/kg^2]
      varnce_rt_2,  & ! Variance of r_t (2nd PDF component)           [kg^2/kg^2]
      thl_1,        & ! Mean of th_l (1st PDF component)                      [K]
      thl_2,        & ! Mean of th_l (2nd PDF component)                      [K]
      varnce_thl_1, & ! Variance of th_l (1st PDF component)                [K^2]
      varnce_thl_2, & ! Variance of th_l (2nd PDF component)                [K^2]
      rrtthl,       & ! Correlation of r_t and th_l (both components)         [-]
      alpha_thl,    & ! Factor relating to normalized variance for th_l       [-]
      alpha_rt,     & ! Factor relating to normalized variance for r_t        [-]
      crt_1,        & ! Coef. on r_t in s/t eqns. (1st PDF comp.)             [-]
      crt_2,        & ! Coef. on r_t in s/t eqns. (2nd PDF comp.)             [-]
      cthl_1,       & ! Coef. on th_l in s/t eqns. (1st PDF comp.)    [(kg/kg)/K]
      cthl_2          ! Coef. on th_l in s/t eqns. (2nd PDF comp.)    [(kg/kg)/K]

    real( kind = core_rknd ) :: &
      chi_1,           & ! Mean of chi (old s) (1st PDF component)       [kg/kg]
      chi_2,           & ! Mean of chi (old s) (2nd PDF component)       [kg/kg]
      stdev_chi_1,     & ! Standard deviation of chi (1st PDF component) [kg/kg]
      stdev_chi_2,     & ! Standard deviation of chi (2nd PDF component) [kg/kg]
      stdev_eta_1,     & ! Standard dev. of eta (old t) (1st PDF comp.)  [kg/kg]
      stdev_eta_2,     & ! Standard dev. of eta (old t) (2nd PDF comp.)  [kg/kg]
      covar_chi_eta_1, & ! Covariance of chi and eta (1st PDF comp.) [kg^2/kg^2]
      covar_chi_eta_2, & ! Covariance of chi and eta (2nd PDF comp.) [kg^2/kg^2]
      corr_chi_eta_1,  & ! Correlation of chi and eta (1st PDF component)    [-]
      corr_chi_eta_2,  & ! Correlation of chi and eta (2nd PDF component)    [-]
      rsatl_1,         & ! Mean of r_sl (1st PDF component)              [kg/kg]
      rsatl_2,         & ! Mean of r_sl (2nd PDF component)              [kg/kg]
      rc_1,            & ! Mean of r_c (1st PDF component)               [kg/kg]
      rc_2,            & ! Mean of r_c (2nd PDF component)               [kg/kg]
      cloud_frac_1,    & ! Cloud fraction (1st PDF component)                [-]
      cloud_frac_2,    & ! Cloud fraction (2nd PDF component)                [-]
      mixt_frac          ! Weight of 1st PDF component (Sk_w dependent)      [-]

    real( kind = core_rknd ) :: & ! If l_use_ADG2 == .true., or l_use_3D_closure
      sigma_sqd_w_1,       & !
      sigma_sqd_w_2,       & !
      sigma_sqd_thl_1,     & !
      sigma_sqd_thl_2,     & !
      sigma_sqd_rt_1,      & !
      sigma_sqd_rt_2

    ! Note:  alpha coefficients = 0.5 * ( 1 - correlations^2 ).
    !        These are used to calculate the scalar widths
    !        varnce_thl_1, varnce_thl_2, varnce_rt_1, and varnce_rt_2 as in Eq. (34)
    !        of Larson and Golaz (2005)

    ! Passive scalar local variables

    real( kind = core_rknd ), dimension(sclr_dim) ::  & 
      sclr1, sclr2,  &
      varnce_sclr1, varnce_sclr2, & 
      alpha_sclr,  & 
      rsclrthl, rsclrrt
!     sclr1_n, sclr2_n,

    logical :: &
      l_scalar_calc, &  ! True if sclr_dim > 0
      l_calc_ice_supersat_frac ! True if we should calculate ice_supersat_frac

    ! Quantities needed to predict higher order moments
    real( kind = core_rknd ) ::  & 
      tl1, tl2,  & 
      beta1, beta2

    real( kind = core_rknd ) :: sqrt_wp2

    ! Thermodynamic quantity

    real( kind = core_rknd ), intent(out) :: rc_coef

    real( kind = core_rknd ) :: &
      wp2rxp,  & ! Sum total < w'^2 r_x' > for all hm species x [(m/s)^2(kg/kg)]
      wprxp,   & ! Sum total < w'r_x' > for all hm species x      [(m/s)(kg/kg)]
      thlprxp, & ! Sum total < th_l'r_x' > for all hm species x       [K(kg/kg)]
      rtprxp     ! Sum total < r_t'r_x' > for all hm species x       [(kg/kg)^2]

    ! variables for a generalization of Chris Golaz' closure
    ! varies width of plumes in theta_l, rt
    real( kind = core_rknd ) :: width_factor_1, width_factor_2
    
    ! variables for the ADG2 and 3D-Luhar closure
    real( kind = core_rknd ) :: big_m_w, small_m_w, &
                                big_m_thl, small_m_thl, &
                                big_m_rt, small_m_rt

    ! variables for computing ice cloud fraction
    real( kind = core_rknd) :: &
      ice_supersat_frac_1, & ! Ice supersaturation fraction (1st PDF comp.)  [-]
      ice_supersat_frac_2, & ! Ice supersaturation fraction (2nd PDF comp.)  [-]
      rt_at_ice_sat1, rt_at_ice_sat2, &
      chi_at_ice_sat1, chi_at_ice_sat2, rc_1_ice, rc_2_ice
    
    ! To test pdf parameters
    real( kind = core_rknd ) ::  &
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

    real( kind = core_rknd ), parameter :: &
      chi_at_liq_sat  = 0.0_core_rknd    ! Always zero

    logical, parameter :: &
      l_liq_ice_loading_test = .false. ! Temp. flag liq./ice water loading test

    integer :: i, hm_idx   ! Indices

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

    err_code = clubb_no_error ! Initialize to the value for no errors

    ! If there is no variance in vertical velocity, then treat rt and theta-l as
    ! constant, as well.  Otherwise width parameters (e.g. varnce_w_1,
    ! varnce_w_2, etc.) are non-zero.
    if ( (wp2 <= w_tol_sqd) .and. (.not. l_use_3D_closure) )  then

      mixt_frac    = one_half
      w_1          = wm
      w_2          = wm
      varnce_w_1   = 0._core_rknd
      varnce_w_2   = 0._core_rknd
      rt_1         = rtm
      rt_2         = rtm
      alpha_rt     = one_half
      varnce_rt_1  = 0._core_rknd
      varnce_rt_2  = 0._core_rknd
      thl_1        = thlm
      thl_2        = thlm
      alpha_thl    = one_half
      varnce_thl_1 = 0._core_rknd
      varnce_thl_2 = 0._core_rknd
      rrtthl       = 0._core_rknd

      if ( l_scalar_calc ) then
        do i = 1, sclr_dim, 1
          sclr1(i)        = sclrm(i)
          sclr2(i)        = sclrm(i)
          varnce_sclr1(i) = 0.0_core_rknd
          varnce_sclr2(i) = 0.0_core_rknd
          alpha_sclr(i)   = one_half
          rsclrrt(i)      = 0.0_core_rknd
          rsclrthl(i)     = 0.0_core_rknd
        end do ! 1..sclr_dim
      end if

    else ! Width (standard deviation) parameters are non-zero

      ! To avoid recomputing
      sqrt_wp2 = sqrt( wp2 )

      if( (.not. l_use_ADG2) .and. (.not. l_use_3D_closure) ) then ! use ADG1
        call  ADG1_w_closure(Skw, wm, wp2, sigma_sqd_w, sqrt_wp2, mixt_frac_max_mag,& 
                             mixt_frac, varnce_w_1, varnce_w_2, w_1_n, w_2_n, w_1, w_2 )

      elseif( l_use_ADG2 ) then ! use ADG2

        ! Reproduce ADG2_w_closure using separate functions
        call calc_Luhar_params( Skw, Skw, &                     ! intent(in)
                                mixt_frac, big_m_w, small_m_w ) ! intent(out)

        call close_Luhar_pdf( wm, wp2, mixt_frac,           & ! intent(in)
                              small_m_w, Skw, Skw,          & ! intent(in)
                              sigma_sqd_w_1, sigma_sqd_w_2, & ! intent(out)
                              varnce_w_1, varnce_w_2,       & ! intent(out)
                              w_1_n, w_2_n, w_1, w_2 )        ! intent(out)

        ! Overwrite sigma_sqd_w for consistency with ADG1
        sigma_sqd_w = min( one / ( one + small_m_w**2 ), 0.99_core_rknd )

      endif ! l_use_ADG2

      if( .not. l_use_3D_closure ) then ! proceed as usual
        ! The normalized variance for thl, rt, and sclr for "plume" 1 is:
        !
        ! { 1 - [1/(1-sigma_sqd_w)]*[ (w'x')^2 / (w'^2 * x'^2) ] / mixt_frac }
        ! * { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) };
        !
        ! where "x" stands for thl, rt, or sclr; "mixt_frac" is the weight of Gaussian
        ! "plume" 1, and 0 <= beta <= 3.
        !
        ! The factor { (1/3)*beta + mixt_frac*( 1 - (2/3)*beta ) } does not depend on
        ! which varable "x" stands for.  The factor is multiplied by 2 and defined
        ! as width_factor_1.
        !
        ! The factor { 1 - [1/(1-sigma_sqd_w)]*[ (w'x')^2 / (w'^2 * x'^2) ] / mixt_frac }
        ! depends on which variable "x" stands for.  It is multiplied by one_half and
        ! defined as alpha_x, where "x" stands for thl, rt, or sclr.

        ! Vince Larson added a dimensionless factor so that the
        ! width of plumes in theta_l, rt can vary.
        ! beta is a constant defined in module parameters_tunable
        !   Set 0<beta<3.
        ! beta=1.5_core_rknd recovers Chris Golaz' simplified formula.
        ! 3 Nov 2003

        width_factor_1 = ( 2.0_core_rknd/3.0_core_rknd )*beta + 2.0_core_rknd&
             *mixt_frac*( one - ( 2.0_core_rknd/3.0_core_rknd )*beta )
        width_factor_2 = 2.0_core_rknd - width_factor_1

        if ( thlp2 <= thl_tol**2 ) then
          thl_1        = thlm
          thl_2        = thlm
          varnce_thl_1 = 0.0_core_rknd
          varnce_thl_2 = 0.0_core_rknd
          alpha_thl    = one_half
        else
!         thl_1_n = - (wpthlp/(sqrt( wp2 )*sqrt( thlp2 )))/w_2_n
!         thl_2_n = - (wpthlp/(sqrt( wp2 )*sqrt( thlp2 )))/w_1_n

          thl_1 = thlm - ( wpthlp/sqrt_wp2 )/w_2_n
          thl_2 = thlm - ( wpthlp/sqrt_wp2 )/w_1_n

          alpha_thl = one_half * ( one - wpthlp*wpthlp / &
             ((one-sigma_sqd_w)*wp2*thlp2) )

          alpha_thl = max( min( alpha_thl, one ), zero_threshold )

          ! Vince Larson multiplied original expressions by width_factor_1,2
          !   to generalize scalar skewnesses.  05 Nov 03
          varnce_thl_1 = ( alpha_thl / mixt_frac * thlp2 ) * width_factor_1
          varnce_thl_2 = ( alpha_thl / (one-mixt_frac) * thlp2 ) * width_factor_2

        end if ! thlp2 <= thl_tol**2

        if ( rtp2 <= rt_tol**2 ) then
          rt_1        = rtm
          rt_2        = rtm
          varnce_rt_1 = 0.0_core_rknd
          varnce_rt_2 = 0.0_core_rknd
          alpha_rt    = one_half
        else
!         rt_1_n = -( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) ) ) / w_2_n
!         rt_2_n = -( wprtp / ( sqrt( wp2 )*sqrt( rtp2 ) ) ) / w_1_n

          rt_1 = rtm - ( wprtp / sqrt_wp2 ) / w_2_n
          rt_2 = rtm - ( wprtp / sqrt_wp2 ) / w_1_n

          alpha_rt = one_half * ( one - wprtp*wprtp / &
             ((one-sigma_sqd_w)*wp2*rtp2) )

          alpha_rt = max( min( alpha_rt, one ), zero_threshold )

          ! Vince Larson multiplied original expressions by width_factor_1,2
          !   to generalize scalar skewnesses.  05 Nov 03
          varnce_rt_1 = ( alpha_rt / mixt_frac * rtp2 ) * width_factor_1
          varnce_rt_2 = ( alpha_rt / (one-mixt_frac) * rtp2 ) * width_factor_2

        end if ! rtp2 <= rt_tol**2

      else ! use 3D_Luhar closure
        
        if ( ( abs(Skw) >= abs(Skthl)) .and. ( abs(Skw) >= abs(Skrt) ) ) then

          ! w has the greatest magnitude of skewness.

          ! Solve for the w PDF
          call calc_Luhar_params( Skw, Skw, &                     ! intent(in)
                                  mixt_frac, big_m_w, small_m_w ) ! intent(out)

          call close_Luhar_pdf( wm, wp2, mixt_frac,           & ! intent(in)
                                small_m_w, Skw, Skw,          & ! intent(in)
                                sigma_sqd_w_1, sigma_sqd_w_2, & ! intent(out)
                                varnce_w_1, varnce_w_2,       & ! intent(out)
                                w_1_n, w_2_n, w_1, w_2 )        ! intent(out)

          ! Solve for the thl PDF
          call backsolve_Luhar_params( Skw, Skthl,         &    ! intent(in)
                                       big_m_w, mixt_frac, &    ! intent(in)
                                       big_m_thl, small_m_thl ) ! intent(out)

          call close_Luhar_pdf( thlm, thlp2, mixt_frac,           &! intent(in)
                                small_m_thl, Skthl, Skw,          &! intent(in)
                                sigma_sqd_thl_1, sigma_sqd_thl_2, &! intent(out)
                                varnce_thl_1, varnce_thl_2,       &! intent(out)
                                thl_1_n, thl_2_n, thl_1, thl_2 )   ! intent(out)

          ! Solve for the rt PDF
          call backsolve_Luhar_params( Skw, Skrt,          &  ! intent(in)
                                       big_m_w, mixt_frac, &  ! intent(in)
                                       big_m_rt, small_m_rt ) ! intent(out)

          call close_Luhar_pdf( rtm, rtp2, mixt_frac,           & ! intent(in)
                                small_m_rt, Skrt, Skw,          & ! intent(in)
                                sigma_sqd_rt_1, sigma_sqd_rt_2, & ! intent(out)
                                varnce_rt_1, varnce_rt_2,       & ! intent(out)
                                rt_1_n, rt_2_n, rt_1, rt_2 )      ! intent(out)

        elseif ( ( abs(Skthl) > abs(Skw) ) &
                   .and. ( abs(Skthl) >= abs(Skrt) )  ) then

          ! theta-l has the greatest magnitude of skewness.

          ! Solve for the thl PDF
          call calc_Luhar_params( Skthl, Skw, &                      !intent(in)
                                  mixt_frac, big_m_thl, small_m_thl )!intent(out)

          ! Solve for the thl PDF
          call close_Luhar_pdf( thlm, thlp2, mixt_frac,           &! intent(in)
                                small_m_thl, Skthl, Skw,          &! intent(in)
                                sigma_sqd_thl_1, sigma_sqd_thl_2, &! intent(out)
                                varnce_thl_1, varnce_thl_2,       &! intent(out)
                                thl_1_n, thl_2_n, thl_1, thl_2 )   ! intent(out)

          ! Solve for the w PDF
          call backsolve_Luhar_params( Skthl, Skw,           & ! intent(in)
                                       big_m_thl, mixt_frac, & ! intent(in)
                                       big_m_w, small_m_w )    ! intent(out)

          call close_Luhar_pdf( wm, wp2, mixt_frac,           & ! intent(in)
                                small_m_w, Skw, Skw,          & ! intent(in)
                                sigma_sqd_w_1, sigma_sqd_w_2, & ! intent(out)
                                varnce_w_1, varnce_w_2,       & ! intent(out)
                                w_1_n, w_2_n, w_1, w_2 )        ! intent(out)

          ! Solve for the rt PDF
          call backsolve_Luhar_params( Skthl, Skrt,           & ! intent(in)
                                       big_m_thl, mixt_frac,  & ! intent(in)
                                       big_m_rt, small_m_rt )   ! intent(out)

          call close_Luhar_pdf( rtm, rtp2, mixt_frac,           & ! intent(in)
                                small_m_rt, Skrt, Skw,          & ! intent(in)
                                sigma_sqd_rt_1, sigma_sqd_rt_2, & ! intent(out)
                                varnce_rt_1, varnce_rt_2,       & ! intent(out)
                                rt_1_n, rt_2_n, rt_1, rt_2 )      ! intent(out)

        else

          ! rt has the greatest magnitude of skewness.

          ! Solve for the rt PDF
          call calc_Luhar_params( Skrt, Skw, &                      ! intent(in)
                                  mixt_frac, big_m_rt, small_m_rt ) ! intent(out)

          ! Solve for the rt PDF
          call close_Luhar_pdf( rtm, rtp2, mixt_frac,           & ! intent(in)
                                small_m_rt, Skrt, Skw,          & ! intent(in)
                                sigma_sqd_rt_1, sigma_sqd_rt_2, & ! intent(out)
                                varnce_rt_1, varnce_rt_2,       & ! intent(out)
                                rt_1_n, rt_2_n, rt_1, rt_2 )      ! intent(out)

          ! Solve for the w PDF
          call backsolve_Luhar_params( Skrt, Skw,           & ! intent(in)
                                       big_m_rt, mixt_frac, & ! intent(in)
                                       big_m_w, small_m_w )   ! intent(out)

          call close_Luhar_pdf( wm, wp2, mixt_frac,           & ! intent(in)
                                small_m_w, Skw, Skw,          & ! intent(in)
                                sigma_sqd_w_1, sigma_sqd_w_2, & ! intent(out)
                                varnce_w_1, varnce_w_2,       & ! intent(out)
                                w_1_n, w_2_n, w_1, w_2 )        ! intent(out)

          ! Solve for the thl PDF
          call backsolve_Luhar_params( Skrt, Skthl,           & ! intent(in)
                                       big_m_rt, mixt_frac,   & ! intent(in)
                                       big_m_thl, small_m_thl ) ! intent(out)

          call close_Luhar_pdf( thlm, thlp2, mixt_frac,           &! intent(in)
                                small_m_thl, Skthl, Skw,          &! intent(in)
                                sigma_sqd_thl_1, sigma_sqd_thl_2, &! intent(out)
                                varnce_thl_1, varnce_thl_2,       &! intent(out)
                                thl_1_n, thl_2_n, thl_1, thl_2 )   ! intent(out)

        endif

        ! CLUBB still uses ADG1 elsewhere in the code. This makes things a
        ! little more consistent.
        sigma_sqd_w = min( one / ( one + small_m_w**2 ), 0.99_core_rknd )

        ! Set to default values when using the 3D_Luhar closure
        alpha_thl = one_half
        alpha_rt = one_half

    endif ! if( .not. l_use_3D_closure )

      ! Compute pdf parameters for passive scalars
      if ( l_scalar_calc ) then
        do i = 1, sclr_dim
          if ( sclrp2(i) <= sclr_tol(i)**2 ) then
            ! Set plume sclr for plume 1,2 to the mean
            sclr1(i)= sclrm(i)
            sclr2(i)= sclrm(i)
            ! Set the variance to zero
            varnce_sclr1(i) = 0.0_core_rknd
            varnce_sclr2(i) = 0.0_core_rknd

            alpha_sclr(i) = one_half
          else
!           sclr1_n(i) = - ( wpsclrp(i) / (sqrt( wp2 ) &
!                        * sqrt( sclrp2(i) )) )/w_2_n
!           sclr2_n(i) = - ( wpsclrp(i) / (sqrt( wp2 ) &
!                        * sqrt( sclrp2(i) )) )/w_1_n

            sclr1(i) = sclrm(i)  & 
                     - ( wpsclrp(i) / sqrt_wp2 ) / w_2_n
            sclr2(i) = sclrm(i)  & 
                     - ( wpsclrp(i) / sqrt_wp2 ) / w_1_n

            alpha_sclr(i) = one_half * ( one - wpsclrp(i)*wpsclrp(i) & 
                    / ((one-sigma_sqd_w)*wp2*sclrp2(i)) )

            alpha_sclr(i) = max( min( alpha_sclr(i), one ), zero_threshold )

            ! Vince Larson multiplied original expressions by width_factor_1,2
            !  to generalize scalar skewnesses.  05 Nov 03
            varnce_sclr1(i) = ( alpha_sclr(i) / mixt_frac * sclrp2(i) ) * width_factor_1
            varnce_sclr2(i) = ( alpha_sclr(i) / (one-mixt_frac) * &
                sclrp2(i) ) * width_factor_2
          end if ! sclrp2(i) <= sclr_tol(i)**2
        end do ! i=1, sclr_dim
      end if ! l_scalar_calc

      ! We include sub-plume correlation with coeff rrtthl.

      if ( varnce_rt_1*varnce_thl_1 > 0._core_rknd .and. &
             varnce_rt_2*varnce_thl_2 > 0._core_rknd ) then
        rrtthl = ( rtpthlp - mixt_frac * ( rt_1-rtm ) * ( thl_1-thlm ) & 
                   - (one-mixt_frac) * ( rt_2-rtm ) * ( thl_2-thlm ) ) & 
                / ( mixt_frac*sqrt( varnce_rt_1*varnce_thl_1 ) &
                   + (one-mixt_frac)*sqrt( varnce_rt_2*varnce_thl_2 ) )
        if ( rrtthl < -one ) then
          rrtthl = -one
        end if
        if ( rrtthl > one ) then
          rrtthl = one
        end if
      else
        rrtthl = 0.0_core_rknd
      end if ! varnce_rt_1*varnce_thl_1 > 0 .and. varnce_rt_2*varnce_thl_2 > 0

      ! Sub-plume correlation, rsclrthl, of passive scalar and theta_l.
      if ( l_scalar_calc ) then
        do i=1, sclr_dim
          if ( varnce_sclr1(i)*varnce_thl_1 > 0._core_rknd .and. &
               varnce_sclr2(i)*varnce_thl_2 > 0._core_rknd ) then
            rsclrthl(i) = ( sclrpthlp(i)  & 
            - mixt_frac * ( sclr1(i)-sclrm(i) ) * ( thl_1-thlm ) & 
            - (one-mixt_frac) * ( sclr2(i)-sclrm(i) ) * ( thl_2-thlm ) ) & 
                / ( mixt_frac*sqrt( varnce_sclr1(i)*varnce_thl_1 )  & 
                         + (one-mixt_frac)*sqrt( varnce_sclr2(i)*varnce_thl_2 ) )
            if ( rsclrthl(i) < -one ) then
              rsclrthl(i) = -one
            end if
            if ( rsclrthl(i) > one ) then
              rsclrthl(i) = one
            end if
          else
            rsclrthl(i) = 0.0_core_rknd
          end if

          ! Sub-plume correlation, rsclrrt, of passive scalar and total water.

          if ( varnce_sclr1(i)*varnce_rt_1 > 0._core_rknd .and. &
               varnce_sclr2(i)*varnce_rt_2 > 0._core_rknd ) then
            rsclrrt(i) = ( sclrprtp(i) - mixt_frac * ( sclr1(i)-sclrm(i) ) * ( rt_1-rtm )&
                         - (one-mixt_frac) * ( sclr2(i)-sclrm(i) ) * ( rt_2-rtm ) ) & 
                       / ( mixt_frac*sqrt( varnce_sclr1(i)*varnce_rt_1 ) &
                         + (one-mixt_frac)*sqrt( varnce_sclr2(i)*varnce_rt_2 ) )
            if ( rsclrrt(i) < -one ) then
              rsclrrt(i) = -one
            end if
            if ( rsclrrt(i) > one ) then
              rsclrrt(i) = one
            end if
          else
            rsclrrt(i) = 0.0_core_rknd
          end if
        end do ! i=1, sclr_dim
      end if ! l_scalar_calc

    end if  ! Widths non-zero

    ! Compute higher order moments (these are interactive)
    wp2rtp  = mixt_frac * ( (w_1-wm)**2+varnce_w_1 ) * ( rt_1-rtm ) & 
            + (one-mixt_frac) * ( (w_2-wm)**2+varnce_w_2 ) * ( rt_2-rtm )

    wp2thlp = mixt_frac * ( (w_1-wm)**2+varnce_w_1 ) * ( thl_1-thlm ) & 
            + (one-mixt_frac) * ( (w_2-wm)**2+varnce_w_2 ) * ( thl_2-thlm )

    ! Compute higher order moments (these are non-interactive diagnostics)
    if ( iwp4 > 0 ) then
      wp4 = mixt_frac * ( 3._core_rknd*varnce_w_1**2 + &
          6._core_rknd*((w_1-wm)**2)*varnce_w_1 + (w_1-wm)**4 ) & 
          + (one-mixt_frac) * ( 3._core_rknd*varnce_w_2**2 + &
          6._core_rknd*((w_2-wm)**2)*varnce_w_2 + (w_2-wm)**4 )
    end if

    if ( iwprtp2 > 0 ) then
      wprtp2  = mixt_frac * ( w_1-wm )*( (rt_1-rtm)**2 + varnce_rt_1 )  & 
              + (one-mixt_frac) * ( w_2-wm )*( (rt_2-rtm)**2 + varnce_rt_2)
    end if

    if ( iwpthlp2 > 0 ) then
      wpthlp2 = mixt_frac * ( w_1-wm )*( (thl_1-thlm)**2 + varnce_thl_1 )  & 
              + (one-mixt_frac) * ( w_2-wm )*( (thl_2-thlm)**2+varnce_thl_2 )
    end if

    if ( iwprtpthlp > 0 ) then
      wprtpthlp = mixt_frac * ( w_1-wm )*( (rt_1-rtm)*(thl_1-thlm)  & 
                + rrtthl*sqrt( varnce_rt_1*varnce_thl_1 ) ) & 
                + ( one-mixt_frac ) * ( w_2-wm )*( (rt_2-rtm)*(thl_2-thlm) & 
                + rrtthl*sqrt( varnce_rt_2*varnce_thl_2 ) )
    end if


    ! Scalar Addition to higher order moments
    if ( l_scalar_calc ) then
      do i=1, sclr_dim

        wp2sclrp(i)  = mixt_frac * ( (w_1-wm)**2+varnce_w_1 )*( sclr1(i)-sclrm(i) ) & 
                     + (one-mixt_frac) * ( (w_2-wm)**2+varnce_w_2 ) * ( sclr2(i)-sclrm(i) )

        wpsclrp2(i) = mixt_frac * ( w_1-wm ) * ( (sclr1(i)-sclrm(i))**2 + varnce_sclr1(i) )  & 
                    + (one-mixt_frac) * ( w_2-wm ) * &
                    ( (sclr2(i)-sclrm(i))**2 + varnce_sclr2(i) )

        wpsclrprtp(i) = mixt_frac * ( w_1-wm ) * ( ( rt_1-rtm )*( sclr1(i)-sclrm(i) )  & 
          + rsclrrt(i)*sqrt( varnce_rt_1*varnce_sclr1(i) ) ) &
          + ( one-mixt_frac )*( w_2-wm ) *  &
            ( ( rt_2-rtm )*( sclr2(i)-sclrm(i) ) + rsclrrt(i)*sqrt( varnce_rt_2*varnce_sclr2(i) ) )

        wpsclrpthlp(i) = mixt_frac * ( w_1-wm ) * ( ( sclr1(i)-sclrm(i) )*( thl_1-thlm )  & 
          + rsclrthl(i)*sqrt( varnce_sclr1(i)*varnce_thl_1 ) ) & 
          + ( one-mixt_frac ) * ( w_2-wm ) * &
            ( ( sclr2(i)-sclrm(i) )*( thl_2-thlm ) &
              + rsclrthl(i)*sqrt( varnce_sclr2(i)*varnce_thl_2 ) )

      end do ! i=1, sclr_dim
    end if ! l_scalar_calc

    ! Compute higher order moments that include theta_v.

    ! First compute some preliminary quantities.
    ! "1" denotes first Gaussian; "2" denotes 2nd Gaussian
    ! liq water temp (Sommeria & Deardorff 1977 (SD), eqn. 3)

    tl1  = thl_1*exner
    tl2  = thl_2*exner

#ifdef GFDL
    if( sclr_dim > 0  .and.  (.not. do_liquid_only_in_clubb) ) then ! h1g, 2010-06-16 begin mod

      if( tl1 > t1_combined ) then
        rsatl_1 = sat_mixrat_liq( p_in_Pa, tl1 )
      elseif( tl1 > t2_combined )  then
        rsatl_1 = sat_mixrat_liq( p_in_Pa, tl1 ) * (tl1 - t2_combined)/(t1_combined - t2_combined) &
             + sat_mixrat_ice( p_in_Pa, tl1 ) * (t1_combined - tl1)/(t1_combined - t2_combined)
      elseif( tl1 > t3_combined )  then
        rsatl_1 = sat_mixrat_ice( p_in_Pa, tl1 ) &
             + sat_mixrat_ice( p_in_Pa, tl1 ) * (RH_crit(1, 1) -one ) &
               * ( t2_combined -tl1)/(t2_combined - t3_combined)
      else
        rsatl_1 = sat_mixrat_ice( p_in_Pa, tl1 ) * RH_crit(1, 1)
      endif

      if( tl2 > t1_combined ) then
        rsatl_2 = sat_mixrat_liq( p_in_Pa, tl2 )
      elseif( tl2 > t2_combined )  then
        rsatl_2 = sat_mixrat_liq( p_in_Pa, tl2 ) * (tl2 - t2_combined)/(t1_combined - t2_combined) &
             + sat_mixrat_ice( p_in_Pa, tl2 ) * (t1_combined - tl2)/(t1_combined - t2_combined)
      elseif( tl2 > t3_combined )  then
        rsatl_2 = sat_mixrat_ice( p_in_Pa, tl2 ) &
             + sat_mixrat_ice( p_in_Pa, tl2 )* (RH_crit(1, 2) -one) &
               * ( t2_combined -tl2)/(t2_combined - t3_combined)
      else
        rsatl_2 = sat_mixrat_ice( p_in_Pa, tl2 ) * RH_crit(1, 2)
      endif

    else !sclr_dim <= 0  or  do_liquid_only_in_clubb = .T.
      rsatl_1 = sat_mixrat_liq( p_in_Pa, tl1 )
      rsatl_2 = sat_mixrat_liq( p_in_Pa, tl2 )

    endif !sclr_dim > 0
#else
    rsatl_1 = sat_mixrat_liq( p_in_Pa, tl1 )
    rsatl_2 = sat_mixrat_liq( p_in_Pa, tl2 ) ! h1g, 2010-06-16 end mod
#endif

    ! SD's beta (eqn. 8)
    beta1 = ep * ( Lv/(Rd*tl1) ) * ( Lv/(Cp*tl1) )
    beta2 = ep * ( Lv/(Rd*tl2) ) * ( Lv/(Cp*tl2) )

    ! s from Lewellen and Yoh 1993 (LY) eqn. 1
    chi_1 = ( rt_1 - rsatl_1 ) / ( one + beta1 * rsatl_1 )
    chi_2 = ( rt_2 - rsatl_2 ) / ( one + beta2 * rsatl_2 )

    ! Coefficients for s'
    ! For each normal distribution in the sum of two normal distributions,
    ! s' = crt * rt'  +  cthl * thl';
    ! therefore, x's' = crt * x'rt'  +  cthl * x'thl'.
    ! Larson et al. May, 2001.

    crt_1  = one/( one + beta1*rsatl_1)
    crt_2  = one/( one + beta2*rsatl_2)

    cthl_1 = ( (one + beta1 * rt_1) / ( one + beta1*rsatl_1)**2 ) & 
             * ( Cp/Lv ) * beta1 * rsatl_1 * exner
    cthl_2 = ( (one + beta2 * rt_2) / ( one + beta2*rsatl_2 )**2 ) & 
             * ( Cp/Lv ) * beta2 * rsatl_2 * exner

    ! Standard deviation of chi for each component.
    ! Include subplume correlation of qt, thl
    ! Because of round-off error,
    ! stdev_chi_1 (and probably stdev_chi_2) can become negative when rrtthl=1
    ! One could also write this as a squared term
    ! plus a postive correction; this might be a neater format
    stdev_chi_1 = sqrt( max( crt_1**2 * varnce_rt_1  &
                          - two * rrtthl * crt_1 * cthl_1  &
                                * sqrt( varnce_rt_1 * varnce_thl_1 )  &
                          + cthl_1**2 * varnce_thl_1,  &
                          zero_threshold )  )

    stdev_chi_2 = sqrt( max( crt_2**2 * varnce_rt_2  &
                          - two * rrtthl * crt_2 * cthl_2  &
                                * sqrt( varnce_rt_2 * varnce_thl_2 )  &
                          + cthl_2**2 * varnce_thl_2,  &
                          zero_threshold )  )

    ! We need to introduce a threshold value for the variance of chi
    if ( stdev_chi_1 <= chi_tol ) then
      ! Treat chi as a delta function in this component.
      stdev_chi_1 = zero
    end if

    if ( stdev_chi_2 <= chi_tol ) then
      ! Treat chi as a delta function in this component.
      stdev_chi_2 = zero
    end if

    ! Standard deviation of eta for each component.
    stdev_eta_1 = sqrt( max( crt_1**2 * varnce_rt_1  &
                          + two * rrtthl * crt_1 * cthl_1  &
                                * sqrt( varnce_rt_1 * varnce_thl_1 )  &
                          + cthl_1**2 * varnce_thl_1,  &
                          zero_threshold )  )

    stdev_eta_2 = sqrt( max( crt_2**2 * varnce_rt_2  &
                          + two * rrtthl * crt_2 * cthl_2  &
                                * sqrt( varnce_rt_2 * varnce_thl_2 )  &
                          + cthl_2**2 * varnce_thl_2,  &
                          zero_threshold )  )

    ! Covariance of chi and eta for each component.
    covar_chi_eta_1 = crt_1**2 * varnce_rt_1 - cthl_1**2 * varnce_thl_1

    covar_chi_eta_2 = crt_2**2 * varnce_rt_2 - cthl_2**2 * varnce_thl_2

    ! Correlation of chi and eta for each component.
    if ( stdev_chi_1 * stdev_eta_1 > zero ) then
      corr_chi_eta_1 = covar_chi_eta_1 / ( stdev_chi_1 * stdev_eta_1 )
    else
      corr_chi_eta_1 = zero
    endif

    if ( stdev_chi_2 * stdev_eta_2 > zero ) then
      corr_chi_eta_2 = covar_chi_eta_2 / ( stdev_chi_2 * stdev_eta_2 )
    else
      corr_chi_eta_2 = zero
    endif
    
    ! Determine whether to compute ice_supersat_frac. We do not compute
    ! ice_supersat_frac for GFDL (unless do_liquid_only_in_clubb is true),
    ! because liquid and ice are both fed into rtm, ruining the calculation.
#ifdef GFDL
    if (do_liquid_only_in_clubb) then
      l_calc_ice_supersat_frac = .true.
    else
      l_calc_ice_supersat_frac = .false.
    end if
#else
    l_calc_ice_supersat_frac = .true.
#endif

    ! Calculate cloud_frac_1 and rc_1
    call calc_cloud_frac_component(chi_1, stdev_chi_1, chi_at_liq_sat, cloud_frac_1, rc_1)

    ! Calculate cloud_frac_2 and rc_2
    call calc_cloud_frac_component(chi_2, stdev_chi_2, chi_at_liq_sat, cloud_frac_2, rc_2)

    if ( l_calc_ice_supersat_frac ) then
      ! We must compute chi_at_ice_sat1 and chi_at_ice_sat2
      if (tl1 <= T_freeze_K) then
        rt_at_ice_sat1 = sat_mixrat_ice( p_in_Pa, tl1 )
        chi_at_ice_sat1 = ( rt_at_ice_sat1 - rsatl_1 ) / ( one + beta1 * rsatl_1 )
      else
        ! If the temperature is warmer than freezing (> 0C) then ice_supersat_frac
        ! is not defined, so we use chi_at_liq_sat
        chi_at_ice_sat1 = chi_at_liq_sat
      end if

      if (tl2 <= T_freeze_K) then
        rt_at_ice_sat2 = sat_mixrat_ice( p_in_Pa, tl2 )
        chi_at_ice_sat2 = ( rt_at_ice_sat2 - rsatl_2 ) / ( one + beta2 * rsatl_2 )
      else
        ! If the temperature is warmer than freezing (> 0C) then ice_supersat_frac
        ! is not defined, so we use chi_at_liq_sat
        chi_at_ice_sat2 = chi_at_liq_sat
      end if

      ! Calculate ice supersaturation fraction in the 1st PDF component.
      call calc_cloud_frac_component( chi_1, stdev_chi_1, chi_at_ice_sat1, &
                                      ice_supersat_frac_1, rc_1_ice )
      
      ! Calculate ice supersaturation fraction in the 2nd PDF component.
      call calc_cloud_frac_component( chi_2, stdev_chi_2, chi_at_ice_sat2, &
                                      ice_supersat_frac_2, rc_2_ice )
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

    rc_coef = Lv / (exner*Cp) - ep2 * thv_ds

    wp2rxp  = zero
    wprxp   = zero
    thlprxp = zero
    rtprxp  = zero
    if ( l_liq_ice_loading_test ) then
       do hm_idx = 1, hydromet_dim, 1
          if ( l_mix_rat_hm(hm_idx) ) then
             wp2rxp  = wp2rxp + wp2hmp(hm_idx)
             wprxp   = wprxp + wphydrometp(hm_idx)
             thlprxp = thlprxp + thlphmp(hm_idx)
             rtprxp  = rtprxp + rtphmp(hm_idx)
          endif
       enddo ! hm_idx = 1, hydromet_dim, 1
    endif ! l_liq_ice_loading_test

    wp2rcp = mixt_frac * ((w_1-wm)**2 + varnce_w_1)*rc_1 &
               + (one-mixt_frac) * ((w_2-wm)**2 + varnce_w_2)*rc_2 & 
             - wp2 * (mixt_frac*rc_1+(one-mixt_frac)*rc_2)

    wp2thvp = wp2thlp + ep1*thv_ds*wp2rtp + rc_coef*wp2rcp - thv_ds * wp2rxp

    wprcp = mixt_frac * (w_1-wm)*rc_1 + (one-mixt_frac) * (w_2-wm)*rc_2

    wpthvp = wpthlp + ep1*thv_ds*wprtp + rc_coef*wprcp - thv_ds * wprxp

    ! Account for subplume correlation in qt-thl
    thlprcp  = mixt_frac * ( (thl_1-thlm)*rc_1 - (cthl_1*varnce_thl_1)*cloud_frac_1 ) & 
             + (one-mixt_frac) * ( (thl_2-thlm)*rc_2 - (cthl_2*varnce_thl_2)*cloud_frac_2 ) & 
             + mixt_frac*rrtthl*crt_1*sqrt( varnce_rt_1*varnce_thl_1 )*cloud_frac_1 & 
             + (one-mixt_frac)*rrtthl*crt_2*sqrt( varnce_rt_2*varnce_thl_2 )*cloud_frac_2
    thlpthvp = thlp2 + ep1*thv_ds*rtpthlp + rc_coef*thlprcp - thv_ds * thlprxp

    ! Account for subplume correlation in qt-thl
    rtprcp = mixt_frac * ( (rt_1-rtm)*rc_1 + (crt_1*varnce_rt_1)*cloud_frac_1 ) & 
           + (one-mixt_frac) * ( (rt_2-rtm)*rc_2 + (crt_2*varnce_rt_2)*cloud_frac_2 ) & 
           - mixt_frac*rrtthl*cthl_1*sqrt( varnce_rt_1*varnce_thl_1 )*cloud_frac_1 & 
           - (one-mixt_frac)*rrtthl*cthl_2*sqrt( varnce_rt_2*varnce_thl_2 )*cloud_frac_2

    rtpthvp  = rtpthlp + ep1*thv_ds*rtp2 + rc_coef*rtprcp - thv_ds * rtprxp

    ! Account for subplume correlation of scalar, theta_v.
    ! See Eqs. A13, A8 from Larson et al. (2002) ``Small-scale...''
    !  where the ``scalar'' in this paper is w.
    if ( l_scalar_calc ) then
      do i=1, sclr_dim
        sclrprcp(i) &
        = mixt_frac * ( ( sclr1(i)-sclrm(i) ) * rc_1 ) &
          + (one-mixt_frac) * ( ( sclr2(i)-sclrm(i) ) * rc_2 ) & 
          + mixt_frac*rsclrrt(i) * crt_1 &
            * sqrt( varnce_sclr1(i) * varnce_rt_1 ) * cloud_frac_1 & 
          + (one-mixt_frac) * rsclrrt(i) * crt_2 &
            * sqrt( varnce_sclr2(i) * varnce_rt_2 ) * cloud_frac_2 & 
          - mixt_frac * rsclrthl(i) * cthl_1 &
            * sqrt( varnce_sclr1(i) * varnce_thl_1 ) * cloud_frac_1 & 
          - (one-mixt_frac) * rsclrthl(i) * cthl_2 &
            * sqrt( varnce_sclr2(i) * varnce_thl_2 ) * cloud_frac_2

        sclrpthvp(i) = sclrpthlp(i) + ep1*thv_ds*sclrprtp(i) + rc_coef*sclrprcp(i)
      end do ! i=1, sclr_dim
    end if ! l_scalar_calc

    ! Compute mean cloud fraction and cloud water
    cloud_frac = calc_cloud_frac(cloud_frac_1, cloud_frac_2, mixt_frac)
    rcm        = mixt_frac * rc_1         + (one-mixt_frac) * rc_2
    
    rcm = max( zero_threshold, rcm )
    
    if (l_calc_ice_supersat_frac) then
      ! Compute ice cloud fraction, ice_supersat_frac
      ice_supersat_frac = calc_cloud_frac( ice_supersat_frac_1, &
                                           ice_supersat_frac_2, mixt_frac )
    else
      ! ice_supersat_frac will be garbage if computed as above
      ice_supersat_frac = 0.0_core_rknd
      if (clubb_at_least_debug_level( 1 )) then
         write(fstderr,*) "Warning: ice_supersat_frac has garbage values if &
                         & do_liquid_only_in_clubb = .false."
      end if
    end if
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
    pdf_params%rrtthl          = rrtthl
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
             wpsclrprtp, wpsclrpthlp, wp2sclrp, &
             err_code )

      ! Error Reporting
      ! Joshua Fasching February 2008

      if ( fatal_error( err_code ) ) then

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

        write(fstderr,*) "level = ", level

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
        write(fstderr,*) "pdf_params%rrtthl = ", pdf_params%rrtthl
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
      if(l_use_3D_closure) then

        ! Means
        wm_clubb_pdf = mixt_frac * w_1   + ( one - mixt_frac ) * w_2

        if( abs( (wm_clubb_pdf - wm) / max(wm,eps) ) > .05_core_rknd ) then
          write(fstderr,*) "wm error at thlm = ", thlm, ( (wm_clubb_pdf - wm) / max(wm,eps) )
        endif

        rtm_clubb_pdf = mixt_frac * rt_1  + ( one - mixt_frac ) * rt_2

        if( abs( (rtm_clubb_pdf - rtm) / max(rtm,eps) ) > .05_core_rknd ) then
          write(fstderr,*) "rtm error at thlm = ", thlm, ( (rtm_clubb_pdf - rtm) / max(rtm,eps) )
        endif

        thlm_clubb_pdf = mixt_frac * thl_1 + ( one - mixt_frac ) * thl_2

        if( abs( (thlm_clubb_pdf - thlm) / thlm ) > .05_core_rknd ) then
          write(fstderr,*) "thlm error at thlm = ", thlm, ( (thlm_clubb_pdf - thlm) / thlm )
        endif

        ! Variances
        if(wp2 > w_tol**2) then

          wp2_clubb_pdf &
          = mixt_frac * ( ( w_1 - wm )**2 + varnce_w_1 ) &
          + ( one - mixt_frac ) * ( ( w_2 - wm )**2 + varnce_w_2 )

          if( ( abs( (wp2_clubb_pdf - wp2) / wp2 ) > .05_core_rknd ) ) then
            write(fstderr,*) "wp2 error at thlm = ", thlm, ( (wp2_clubb_pdf - wp2) / wp2 )
          endif

        endif

        if(rtp2 > rt_tol**2) then

          rtp2_clubb_pdf &
          = mixt_frac * ( ( rt_1 - rtm )**2 + varnce_rt_1 ) &
          + ( one - mixt_frac ) * ( ( rt_2 - rtm )**2 + varnce_rt_2 )

          if( abs( (rtp2_clubb_pdf - rtp2) / rtp2 ) > .05_core_rknd ) then
            write(fstderr,*) "rtp2 error at thlm = ", thlm, &
            "Error = ", ( (rtp2_clubb_pdf - rtp2) / rtp2 )
          endif

        endif

        if(thlp2 > thl_tol**2) then

          thlp2_clubb_pdf &
          = mixt_frac * ( ( thl_1 - thlm )**2 + varnce_thl_1 ) &
          + ( one - mixt_frac ) * ( ( thl_2 - thlm )**2 + varnce_thl_2 )

          if( abs( (thlp2_clubb_pdf - thlp2) / thlp2 ) > .05_core_rknd ) then
            write(fstderr,*) "thlp2 error at thlm = ", thlm, &
            "Error = ", ( (thlp2_clubb_pdf - thlp2) / thlp2 )
          endif

        endif

        ! Third order moments
        wp3_clubb_pdf &
        = mixt_frac * ( w_1 - wm ) &
                    * ( ( w_1 - wm )**2 + 3.0_core_rknd*varnce_w_1 ) &
          + ( one - mixt_frac ) * ( w_2 - wm ) &
                          * ( ( w_2 - wm )**2 + 3.0_core_rknd*varnce_w_2 )

        rtp3_clubb_pdf &
        = mixt_frac * ( rt_1 - rtm ) &
                    * ( ( rt_1 - rtm )**2 + 3.0_core_rknd*varnce_rt_1 ) &
          + ( one - mixt_frac ) * ( rt_2 - rtm ) &
                                * ( ( rt_2 - rtm )**2 + 3.0_core_rknd*varnce_rt_2 )

        thlp3_clubb_pdf &
        = mixt_frac * ( thl_1 - thlm ) &
                    * ( ( thl_1 - thlm )**2 + 3.0_core_rknd*varnce_thl_1 ) &
          + ( one - mixt_frac ) * ( thl_2 - thlm ) &
                                * ( ( thl_2 - thlm )**2 + 3.0_core_rknd*varnce_thl_2 )

        ! Skewness
        Skw_clubb_pdf = wp3_clubb_pdf / &
        ( wp2_clubb_pdf + Skw_denom_coef * w_tol**2 )**1.5_core_rknd

        if(Skw > .05_core_rknd) then
          if( abs( (Skw_clubb_pdf - Skw) / Skw ) > .25_core_rknd ) then
            write(fstderr,*) "Skw error at thlm = ", thlm, &
            "Error = ",( (Skw_clubb_pdf - Skw) / Skw ), Skw_clubb_pdf, Skw
          endif
        endif

        Skrt_clubb_pdf = rtp3_clubb_pdf / &
        ( rtp2_clubb_pdf + Skw_denom_coef * rt_tol**2 )**1.5_core_rknd

        if(Skrt > .05_core_rknd) then
          if( abs( (Skrt_clubb_pdf - Skrt) / Skrt ) > .25_core_rknd ) then
            write(fstderr,*) "Skrt error at thlm = ", thlm, &
              "Error = ", ( (Skrt_clubb_pdf - Skrt) / Skrt ), Skrt_clubb_pdf, Skrt
          endif
        endif

        Skthl_clubb_pdf = thlp3_clubb_pdf / &
        ( thlp2_clubb_pdf + Skw_denom_coef * thl_tol**2 )**1.5_core_rknd

        if(Skthl > .05_core_rknd) then
          if( abs( (Skthl_clubb_pdf - Skthl) / Skthl ) > .25_core_rknd ) then
            write(fstderr,*) "Skthl error at thlm = ", thlm, &
              "Error = ", ( (Skthl_clubb_pdf - Skthl) / Skthl ), Skthl_clubb_pdf, Skthl
          endif
        endif

      end if !l_use_3D_closure

    end if ! clubb_at_least_debug_level

    return
  end subroutine pdf_closure
  
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
    !-----------------------------------------------------------------------
    
    use constants_clubb, only: &
        chi_tol,  & ! Tolerance for pdf parameter chi       [kg/kg]
        sqrt_2pi, & ! sqrt(2*pi)
        sqrt_2,   & ! sqrt(2)
        one,      & ! 1
        one_half, & ! 1/2
        zero        ! 0

    use anl_erf, only:  & 
        erf ! Procedure(s) -- The error function

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      mean_chi_i,  & ! Mean of chi (old s) (ith PDF component)           [kg/kg]
      stdev_chi_i, & ! Standard deviation of chi (ith PDF component)     [kg/kg]
      chi_at_sat     ! Value of chi at saturation (0--liquid; neg.--ice) [kg/kg]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      cloud_frac_i, & ! Cloud fraction (ith PDF component)               [-]
      rc_i            ! Mean cloud water mixing ratio (ith PDF comp.)    [kg/kg]

    ! Local Variables
    real( kind = core_rknd) :: zeta_i

    !----- Begin Code -----
    if ( stdev_chi_i > chi_tol ) then

       ! The value of chi varies in the ith PDF component.

       zeta_i = ( mean_chi_i - chi_at_sat ) / stdev_chi_i

       cloud_frac_i = one_half * ( one + erf( zeta_i / sqrt_2 )  )

       rc_i = ( mean_chi_i - chi_at_sat ) * cloud_frac_i &
              + stdev_chi_i * exp( - one_half * zeta_i**2 ) / ( sqrt_2pi )

    else ! stdev_chi_i <= chi_tol

       ! The value of chi does not vary in the ith PDF component.
       if ( ( mean_chi_i - chi_at_sat ) < zero ) then
          ! All clear air in the ith PDF component.
          cloud_frac_i = zero
          rc_i         = zero
       else ! mean_chi_i >= 0
          ! All cloud in the ith PDF component.
          cloud_frac_i = one
          rc_i         = mean_chi_i - chi_at_sat
       endif ! mean_chi_i < 0

    endif ! stdev_chi_i > chi_tol


    return
    
  end subroutine calc_cloud_frac_component

  !=============================================================================
  function calc_cloud_frac( cloud_frac_1, cloud_frac_2, mixt_frac )

  ! Description:
  !   Given the the two pdf components of a cloud fraction, and the weight
  !   of the first component, this fuction calculates the cloud fraction,
  !   cloud_frac
  !
  ! References:
  !-----------------------------------------------------------------------
    
    use constants_clubb, only: & ! Constant(s)
        one,            & ! 1
        fstderr,        & ! Standard error output
        zero_threshold    ! A physical quantity equal to zero
    
    use clubb_precision, only: &
        core_rknd        ! Precision
      
    use error_code, only: &
        clubb_at_least_debug_level ! Function to check whether clubb is in
                                   !  at least the specified debug level 
      
    implicit none
    
    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      cloud_frac_1, & ! First PDF component of cloud_frac
      cloud_frac_2, & ! Second PDF component of cloud_frac
      mixt_frac       ! Weight of 1st PDF component (Sk_w dependent)
    
    ! Output Variables
    real( kind = core_rknd) :: &
      calc_cloud_frac ! Cloud fraction
    
    ! Local Variables
    real( kind = core_rknd) :: &
      cloud_frac      ! Cloud fraction (used as a holding variable for
                      !                    output)

  !-----------------------------------------------------------------------
    !----- Begin Code -----
    cloud_frac = mixt_frac * cloud_frac_1 + (one-mixt_frac) * cloud_frac_2
    
    ! Note: Brian added the following lines to ensure that there
    ! are never any negative liquid water values (or any negative
    ! cloud fraction values, for that matter).  According to
    ! Vince Larson, the analytic formula should not produce any
    ! negative results, but such computer-induced errors such as
    ! round-off error may produce such a value.  This has been
    ! corrected because Brian found a small negative value of
    ! rcm in the first timestep of the FIRE case.

    cloud_frac  = max( zero_threshold, cloud_frac )
    if ( clubb_at_least_debug_level( 2 ) ) then
      if ( cloud_frac > one ) then
        write(fstderr,*) "Cloud fraction > 1"
      end if
    end if
    cloud_frac = min( one, cloud_frac )

    calc_cloud_frac = cloud_frac
    return
    
  end function calc_cloud_frac

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
  elemental subroutine ADG1_w_closure(Skw, wm, wp2, sigma_sqd_w, sqrt_wp2, mixt_frac_max_mag,&
                                      mixt_frac, varnce_w_1, varnce_w_2, w_1_n, w_2_n, &
                                      w_1, w_2 )               
  ! Description:
  !   The Analytic Double Gaussian 1 closure is used by default in CLUBB. It
  !   assumes the widths of both w Gaussians to be the same.
  !
  ! References:
  !   Golaz, J-C., V. E. Larson, and W. R. Cotton, 2002a: A PDF-based model for
  !   boundary layer clouds. Part I: Method and model description. J. Atmos.
  !   Sci., 59, 3540–3551.
  !
  !   Vincent E. Larson and Jean-Christophe Golaz, 2005: Using Probability
  !   Density Functions to Derive Consistent Closure Relationships among
  !   Higher-Order Moments. Mon. Wea. Rev., 133, 1023–1042.
  !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one, &
        one_half         

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skw,              & ! Skewness of w                         [-] 
      wm,               & ! Mean w                                [m / s] 
      wp2,              & ! w'^2                                  [m^2/s^2]
      sigma_sqd_w,      & ! Widths of each w Gaussian             [-]
      sqrt_wp2,         & ! w'                                    [m/s^1]
      mixt_frac_max_mag   !                                       [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mixt_frac,    & ! Mixture fraction                          [-]
      varnce_w_1,   & ! Variance of w (1st PDF component)         [m^2/s^2]
      varnce_w_2,   & ! Variance of w (2nd PDF component)         [m^2/s^2] 
      w_1_n,        & ! Normalized Mean of w (1st PDF component)  [-]  
      w_2_n,        & ! Normalized Mean of w (2nd PDF component)  [-]
      w_1,          & ! Mean of w (1st PDF component)             [m/s] 
      w_2             ! Mean of w (2nd PDF component)             [m/s]     

  !-----------------------------------------------------------------------
    !----- Begin Code -----

      ! The variable "mixt_frac" is the weight of the 1st PDF component.  The
      ! weight of the 2nd PDF component is "1-mixt_frac".  If there isn't any
      ! skewness of w (Sk_w = 0 because w'^3 = 0), mixt_frac = 0.5, and both
      ! PDF components are equally weighted.  If there is positive skewness of
      ! w (Sk_w > 0 because w'^3 > 0), 0 < mixt_frac < 0.5, and the 2nd PDF
      ! component has greater weight than does the 1st PDF component.  If there
      ! is negative skewness of w (Sk_w < 0 because w'^3 < 0),
      ! 0.5 < mixt_frac < 1, and the 1st PDF component has greater weight than
      ! does the 2nd PDF component.
       if ( abs( Skw ) <= 1e-5_core_rknd ) then
          mixt_frac = one_half
       else
          mixt_frac = one_half * ( one - Skw/ &
             sqrt( 4.0_core_rknd*( one - sigma_sqd_w )**3 + Skw**2 ) )
       endif

      ! Clip mixt_frac, 1-mixt_frac, to avoid dividing by zero
      ! Formula for mixt_frac_max_mag =
      ! 1 - ( 1/2 * ( 1 - Skw_max/sqrt( 4*( 1 - sigma_sqd_w )^3 + Skw_max^2 ) ) )
      ! Where sigma_sqd_w is fixed at 0.4.
      mixt_frac = min( max( mixt_frac, one-mixt_frac_max_mag ), mixt_frac_max_mag )

      ! The normalized mean of w for Gaussian "plume" 1 is w_1_n.  It's value
      ! will always be greater than 0.  As an example, a value of 1.0 would
      ! indicate that the actual mean of w for Gaussian "plume" 1 is found
      ! 1.0 standard deviation above the overall mean for w.
      w_1_n = sqrt( ( (one-mixt_frac)/mixt_frac )*(one-sigma_sqd_w) )
      ! The normalized mean of w for Gaussian "plume" 2 is w_2_n.  It's value
      ! will always be less than 0.  As an example, a value of -0.5 would
      ! indicate that the actual mean of w for Gaussian "plume" 2 is found
      ! 0.5 standard deviations below the overall mean for w.
      w_2_n = -sqrt( ( mixt_frac/(one-mixt_frac) )*(one-sigma_sqd_w) )
      ! The mean of w for Gaussian "plume" 1 is w_1.
      w_1 = wm + sqrt_wp2*w_1_n
      ! The mean of w for Gaussian "plume" 2 is w_2.
      w_2 = wm + sqrt_wp2*w_2_n

      ! The variance of w for Gaussian "plume" 1 for varnce_w_1.
      varnce_w_1  = sigma_sqd_w*wp2
      ! The variance of w for Gaussian "plume" 2 for varnce_w_2.
      ! The variance in both Gaussian "plumes" is defined to be the same.
      varnce_w_2  = sigma_sqd_w*wp2

  end subroutine ADG1_w_closure

  !=============================================================================
  elemental subroutine calc_Luhar_params( Skx, Skw, &
                                          mixt_frac, big_m, small_m )

    ! Description:
    ! For the Luhar closure, this subroutine takes Skx (and Skw) as input and
    ! outputs the mixture fraction, big_m, and small_m. This code was written
    ! using the equations and nomenclature of Larson et al. (2002) Appendix
    ! section e.
    !
    ! The relationship between skewness of x (Skx), mixture fraction (a), and
    ! Luhar's small m (m) is given by:
    !
    ! Skx^2 = ( m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3 )
    !         * ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).
    !
    ! Luhar's large M (M) is used to more easily express the factor involving
    ! the m's:
    !
    ! M = ( m^2 + 1 )^3 / ( m^2 * ( m^2 + 3 )^2 ).
    !
    ! The equation involving skewness of x becomes:
    !
    ! Skx^2 = ( 1 / M ) * ( 1 - 2*a )^2 / ( a * ( 1 - a ) );
    !
    ! or:
    !
    ! M * Skx^2 = ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).
    !
    ! This equation can be rewritten as:
    !
    ! ( a * ( 1 - a ) ) * M * Skx^2 = ( 1 - 2*a )^2;
    !
    ! as well as:
    !
    ! ( a - a^2 ) * M * Skx^2 = 1 - 4*a + 4*a^2;
    !
    ! and eventually as:
    !
    ! ( 4 + M * Skx^2 ) * a^2 - ( 4 + M * Skx^2 ) * a + 1 = 0.
    !
    ! Solving the quadratic equation for a:
    !
    ! a = (1/2) * ( 1 +- Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).
    !
    ! Since by definition, mu_w_1 >= mu_w_2, a < 0.5 when Skw > 0, the equation
    ! for mixture fraction is:
    !
    ! a = (1/2) * ( 1 - Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).
    !
    ! For 3-D Luhar, the variable (w, rt, or theta-l) with the greatest
    ! magnitude of skewness is used to calculate mixture fraction.  Since it is
    ! desirable to still have a < 0.5 when Skw > 0 and a > 0.5 when Skw < 0, the
    ! sign function is used.  The value of Skx is replaced by:
    !
    ! Skx|_adj = sign(Skw) * sign(Skx) * Skx;
    !
    ! where
    !
    ! sign(Skx) = | 1 when x >= 0
    !             | -1 when x < 0.
    !
    ! Since Skx|_adj^2 = ( sign(Skw) * sign(Skx) * Skx )^2
    ! = ( sign(Skw) * sign(Skx) )^2 * Skx^2 = Skx^2, the equation for mixture
    ! fraction is:
    !
    ! a = (1/2)
    !     * ( 1 - sign(Skw) * sign(Skx) * Skx * sqrt( 1 / ( 4/M + Skx^2 ) ) ).
    !
    ! When using the ADG2 closure or when using the 3-D Luhar closure when the
    ! variable with the greatest magnitude of skewness is w, Skw = Skx and
    ! sign(Skw) * sign(Skx) is always equal to 1, reducing the equation to its
    ! previous form.

    ! References:
    !    Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
    !    Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
    !    Probability Density Functions. J. Atmos. Sci., 59, 3519–3539.
    !
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        four,       & ! Constant(s)
        three,      &
        one,        &
        two_thirds, &
        one_half,   &
        one_third,  &
        zero

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skx, & ! Skewness of x ( <x'^3> / <x'2>^(3/2) )                     [-]
      Skw    ! Skewness of w ( <w'^3> / <w'2>^(3/2) )                     [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mixt_frac, & ! Mixture fraction                                     [-]
      big_m,     & ! Luhar's M                                            [-]
      small_m      ! Luhar's m                                            [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      small_m_sqd, & ! Luhar's m^2                                        [-]
      sign_Skw,    & ! Sign( Skw ); 1 when Skw >= 0 or -1 when Skw < 0    [-]
      sign_Skx       ! Sign( Skx ); 1 when Skx >= 0 or -1 when Skx < 0    [-]


    ! Calculate Luhar's m (small m).
    ! If Skx is very small, then small_m will tend to zero which risks
    ! divide-by-zero.  To ameliorate this problem, we enforce abs( x_1_n ) and
    ! abs( x_2_n ) > 0.05.
    ! Note:  Luhar's small_m (m) is the only tunable parameter in the Luhar
    !        closure, so this equation can be changed.  However, the value of m
    !        should go toward 0 as Skx goes toward 0 so that the double Gaussian
    !        reduces to a single Gaussian when the distribution is unskewed.
    small_m = max( two_thirds * abs( Skx )**one_third, 0.05_core_rknd )

    ! Calculate m^2.
    small_m_sqd = small_m**2

    ! Calculate Luhar's M (big M).
    big_m = ( one + small_m_sqd )**3 &
            / ( ( three + small_m_sqd )**2 * small_m_sqd )

    ! Calculate sign( Skw ).
    if ( Skw >= zero ) then
       sign_Skw = one
    else ! Skw < 0
       sign_Skw = -one
    endif ! Skw >= 0

    ! Calculate sign( Skx ).
    if ( Skx >= zero ) then
       sign_Skx = one
    else ! Skx < 0
       sign_Skx = -one
    endif ! Skx >= 0

    ! Calculate mixture fraction.
    mixt_frac = one_half &
                * ( one - sign_Skw * sign_Skx * Skx &
                          * sqrt( one / ( ( four / big_m ) + Skx**2 ) ) )


    return

  end subroutine calc_Luhar_params

  !=============================================================================
  elemental subroutine close_Luhar_pdf( xm, xp2, mixt_frac, &
                                        small_m, Skx, Skw, &
                                        sigma_sqd_x_1, sigma_sqd_x_2, &
                                        varnce_x_1, varnce_x_2, &
                                        x_1_n, x_2_n, x_1, x_2 )

    ! Description:
    ! For the Luhar closure, this subroutine takes Skx, xm, xp2, and mixt_frac,
    ! big_m, and small_m (calculated in calc_Luhar_params) as input and outputs
    ! the PDF component means and variances of a variable x in the joint-PDF
    ! according to Luhar et al. (1996).  This code was written using the
    ! equations and nomenclature of Larson et al. (2002) Appendix section e.

    ! References:
    !    Vincent E. Larson, Jean-Christophe Golaz, and William R. Cotton, 2002:
    !    Small-Scale and Mesoscale Variability in Cloudy Boundary Layers: Joint
    !    Probability Density Functions. J. Atmos. Sci., 59, 3519–3539.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,      & ! Constant(s)
        one_half, &
        zero,     &
        eps

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,        & ! Mean (overall) of x, <x>                      [(x units)]
      xp2,       & ! Variance (overall) of x, <x'^2>               [(x units)^2]
      mixt_frac, & ! Mixture fraction                              [-]
      small_m,   & ! Luhar's small m                               [-]
      Skx,       & ! Skewness of x ( <x'^3> / <x'2>^(3/2) )        [-]
      Skw          ! Skewness of w ( <w'^3> / <w'2>^(3/2) )        [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      sigma_sqd_x_1, & ! Normalized width parameter of x (1st PDF component) [-]
      sigma_sqd_x_2, & ! Normalized width parameter of x (1st PDF component) [-]
      varnce_x_1,    & ! Variance of x (1st PDF component)         [(x units)^2]
      varnce_x_2,    & ! Variance of x (2nd PDF component)         [(x units)^2]
      x_1_n,         & ! Normalized mean of x (1st PDF component)            [-]
      x_2_n,         & ! Normalized mean of x (2nd PDF component)            [-]
      x_1,           & ! Mean of x (1st PDF component)               [(x units)]
      x_2              ! Mean of x (2nd PDF component)               [(x units)]

    ! Local Variables
    real( kind = core_rknd) :: &
      sqrt_xp2, & ! Square root of the variance of x                 [(x units)]
      sign_Skw, & ! Sign( Skw ); 1 when Skw >= 0 or -1 when Skw < 0  [-]
      sign_Skx    ! Sign( Skx ); 1 when Skx >= 0 or -1 when Skx < 0  [-]


    ! Calculate sign( Skw ).
    if ( Skw >= zero ) then
       sign_Skw = one
    else ! Skw < 0
       sign_Skw = -one
    endif ! Skw >= 0

    ! Calculate sign( Skx ).
    if ( Skx >= zero ) then
       sign_Skx = one
    else ! Skx < 0
       sign_Skx = -one
    endif ! Skx >= 0

    ! Calculate the square root of the overall variance of x.
    sqrt_xp2 = sqrt( xp2 )

    ! Normalized width parameter of x in the 1st PDF component.
    sigma_sqd_x_1 = ( one - mixt_frac ) / ( mixt_frac * ( one + small_m**2 ) )

    ! The variance of x in the 1st PDF component.
    varnce_x_1 = sigma_sqd_x_1 * xp2

    ! Normalized width parameter of x in the 2nd PDF component.
    sigma_sqd_x_2 = mixt_frac / ( ( one - mixt_frac ) * ( one + small_m**2 ) )

    ! The variance of x in the 2nd PDF component.
    varnce_x_2 = sigma_sqd_x_2 * xp2

    ! Normalized mean of x in the 1st PDF component.
    x_1_n = sign_Skw * sign_Skx * small_m * sqrt( sigma_sqd_x_1 )

    ! Normalized mean of x in the 2nd PDF component.
    x_2_n = -sign_Skw * sign_Skx * small_m * sqrt( sigma_sqd_x_2 )

    ! The mean of x in the 1st PDF component.
    x_1 = xm + sqrt_xp2 * x_1_n

    ! The mean of x in the 2nd PDF component.
    x_2 = xm + sqrt_xp2 * x_2_n


    return

  end subroutine close_Luhar_pdf

  !=============================================================================
  elemental subroutine backsolve_Luhar_params( Sk_max, Skx, &
                                               big_m_max, mixt_frac, &
                                               big_m_x, small_m_x )

    ! Description:
    ! This subroutine calculates Luhar's big_m and small_m for the variate 'x'
    ! consistent with the mixture fraction of the variate with the largest
    ! skewness.
    !
    ! The relationship between skewness of x (Skx), mixture fraction (a), and
    ! Luhar's small m (m) is given by:
    !
    ! Skx^2 = ( m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3 )
    !         * ( 1 - 2*a )^2 / ( a * ( 1 - a ) ).
    !
    ! Moving the factor involving mixture fraction to the right-hand side:
    !
    ! ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2
    ! = m^2 * ( m^2 + 3 )^2 / ( m^2 + 1 )^3.
    !
    ! This can be rewritten as:
    !
    ! ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2
    ! = ( m^6 + 6*m^4 + 9*m^2 ) / ( m^6 + 3*m^4 + 3*m^2 + 1 ).
    !
    ! Setting alpha = ( ( a * ( 1 - a ) ) / ( 1 - 2*a )^2 ) * Skx^2, the
    ! equation can be rewritten as:
    !
    ! ( m^6 + 3*m^4 + 3*m^2 + 1 ) * alpha = m^6 + 6*m^4 + 9*m^2.
    !
    ! This can be rearranged and rewritten as:
    !
    ! ( alpha - 1 ) * m^6 + ( 3 * alpha - 6 ) * m^4
    ! + ( 3 * alpha - 9 ) * m^2 + alpha = 0.
    !
    ! This can be rewritten again as:
    !
    ! ( alpha - 1 ) * (m^2)^3 + ( 3 * alpha - 6 ) * (m^2)^2
    ! + ( 3 * alpha - 9 ) * (m^2) + alpha = 0.
    !
    ! The goal is to solve for m^2, and then take the square root of m^2 to
    ! solve for m.  This can be accomplished by using the cubic formula (with
    ! the l_use_cubic_backsolve option), or else by a quadratic approximation.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three,    &
        two,      &
        one,      &
        one_half, &
        zero,     &
        eps,      &
        fstderr

    use clubb_precision, only: &
        core_rknd     ! Precision

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Sk_max,    & ! Maximum skewness
      Skx,       & ! Skewness of the variate solving small_m and big_m for
      big_m_max, & ! Luhar's big_m of the variate with maximum skewness
      mixt_frac    ! Mixture fraction                                      [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      big_m_x,  & ! Luhar's big_m for the variate being solved for
      small_m_x   ! Luhar's small_m for the variate being solved for

    ! Local Variables
    real( kind = core_rknd ) :: &
      alpha,     & ! 1 / big_m_x
      a,         & ! For readability, quadratic equation
      b,         &
      c,         &
      alpha_upr, &
      alpha_low, &
      discrim

    ! Flag to backsolve for m^2 using cubic formula
    logical, parameter :: &
      l_use_cubic_backsolve = .true.


    if ( l_use_cubic_backsolve ) then

       if ( abs( mixt_frac - one_half ) < 0.001_core_rknd ) then

          ! When mixture fraction = 0.5 (based on the variable with the largest
          ! magnitude of skewness), all variables must have a skewness of 0.
          ! Set m to the minimum threshold of 0.05.
          small_m_x = 0.05_core_rknd

          ! Calculate the corresponding value of big_m_x.
          big_m_x = ( one + small_m_x**2 )**3 &
                    / ( ( three + small_m_x**2 )**2 * small_m_x**2 )

       elseif ( Skx == zero ) then

          ! Mixture fraction /= 0.5 because the variable with the largest
          ! magnitude of skewness has a skewness /= 0.  However, variable x has
          ! a skewness of 0.  In order to reproduce the correct skewness for
          ! variable x, set m to 0 (regardless of minimum thresholds used in
          ! other parts of the code).
          small_m_x = zero

          ! The value of big_m_x should be inf.  Set it to huge.  This is not
          ! used in any calculation, anyway.
          big_m_x = huge( big_m_x )

       else  ! mixt_frac /= 0.5 and Skx /= 0

          ! Backsolve for m, given mixt_frac and Skx.

          ! alpha = 1/M is given by:
          ! [ mixt_frac * ( 1 - mixt_frac ) / ( 1 - 2 * mixt_frac )^2 ] * Skx^2.
          alpha = ( mixt_frac * ( one - mixt_frac ) &
                    / ( one - two * mixt_frac )**2 ) * Skx**2

          ! Calculate big_m_x.
          big_m_x = one / alpha

          ! Solve the cubic equation for m^2:
          ! ( alpha - 1 ) * (m^2)^3 + ( 3 * alpha - 6 ) * (m^2)^2
          ! + ( 3 * alpha - 9 ) * (m^2) + alpha = 0.
          ! The largest root is preferred.
          small_m_x &
          = sqrt( max( max_cubic_root( alpha - one, three * alpha - 6.0_core_rknd, &
                                       three * alpha - 9.0_core_rknd, alpha ), &
                       0.05_core_rknd**2 ) )

       endif

    else ! original formualation

       alpha = ( Skx**2 / (max(Sk_max**2, eps) * big_m_max) )  ! 1 / big_m_x

       ! This limit keeps the discriminant >= 0
       alpha_upr = 2.0_core_rknd*sqrt( 13.0_core_rknd ) - 5.0_core_rknd

       alpha_low = eps

       ! For this approximation, alpha must be less than 2*sqrt(13) - 5 to get a real ans.
       alpha = min(alpha, alpha_upr)

       ! For testing, eliminate possibility of divide by zero
       alpha = max(alpha,alpha_low)

       ! Use a piece-wise approximation
       if(alpha < 1.0_core_rknd) then
         a = max(3.0_core_rknd * alpha - 6.0_core_rknd, eps) ! Prevent divide by zero
         b = 3.0_core_rknd * alpha - 9.0_core_rknd
         c = alpha

         discrim = b**2 - 4.0_core_rknd * a * c
         small_m_x = sqrt( (-b - sqrt(discrim)) / (2.0_core_rknd * a) )
       else
         ! For this approximation, alpha must be less than 2*sqrt(13) - 5 to get a real ans.
         alpha = min(alpha, 2.0_core_rknd)

         a = max(6.0_core_rknd * alpha - 9.0_core_rknd, eps) ! Prevent divide by zero
         b = -6.0_core_rknd
         c = 2.0_core_rknd * alpha - 1.0_core_rknd

         discrim = b**2 - 4.0_core_rknd * a * c
         small_m_x = sqrt( (-b - sqrt(discrim)) / (2.0_core_rknd * a) )
       endif

       ! Clip consistently with subroutine calc_Luhar_params
       small_m_x = max( 5e-2_core_rknd, small_m_x)

       big_m_x = 1.0_core_rknd / alpha

    endif ! l_use_cubic_backsolve


  end subroutine backsolve_Luhar_params

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
  pure function max_cubic_root( a_coef, b_coef, c_coef, d_coef ) &
  result( max_root )

    ! Description:
    ! Calculates the largest root that results from solving a cubic equation of
    ! the form a*x^3 + b*x^2 + c*x + d = 0.
    !
    ! This is done to backsolve for m^2 for the 3-D Luhar closure, given the
    ! values of mixt_frac and Skx.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        zero    ! Constant(s)

    use calc_roots, only: &
        cubic_solve,     & ! Procedure(s)
        quadratic_solve

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      a_coef, & ! Coefficient a (of x^3) in a*x^3 + b*x^2 + c^x + d = 0    [-]
      b_coef, & ! Coefficient b (of x^2) in a*x^3 + b*x^2 + c^x + d = 0    [-]
      c_coef, & ! Coefficient c (of x) in a*x^3 + b*x^2 + c^x + d = 0      [-]
      d_coef    ! Coefficient d in a*x^3 + b*x^2 + c^x + d = 0             [-]

    ! Return Variable
    real( kind = core_rknd ) :: &
      max_root    ! Maximum root that solves the cubic equation            [-]

    ! Local Variables
    complex( kind = core_rknd ), dimension(3) :: &
      cubic_roots    ! Roots of x that satisfy a*x^3 + b*x^2 + c*x + d = 0 [-]

    complex( kind = core_rknd ), dimension(2) :: &
      quadratic_roots    ! Roots of x that satisfy b*x^2 + c*x + d = 0     [-]

    real( kind = core_rknd ) :: &
      a_coef_thresh, & ! Minimum threshold of |a| to use cubic solver      [-]
      b_coef_thresh    ! Minimum threshold of |b| to use quadratic solver  [-]


    ! Calculate a minimum threshold for |a| to call this a cubic equation.
    a_coef_thresh = 0.001_core_rknd &
                    * max( abs(b_coef), abs(c_coef), abs(d_coef) )

    ! Calculate a minimum threshold for |b| to call this a quadratic equation.
    ! This only matters when |a| <= a_coef_thresh.
    b_coef_thresh = 0.001_core_rknd * max( abs(c_coef), abs(d_coef) )

    if ( abs( a_coef ) > a_coef_thresh ) then

       ! The equation is a cubic equation.
       cubic_roots = cubic_solve( a_coef, b_coef, c_coef, d_coef )

       if ( aimag( cubic_roots(2) ) == zero  &
            .and. aimag( cubic_roots(3) ) == zero ) then

          ! Find the maximum root of the three roots.
          max_root = max( real( cubic_roots(1), kind = core_rknd ), &
                          real( cubic_roots(2), kind = core_rknd ), &
                          real( cubic_roots(3), kind = core_rknd ) )

       else  ! cubic_roots(2) and cubic_roots(3) are complex.

          max_root = real( cubic_roots(1), kind = core_rknd )

       endif

    elseif ( abs( b_coef ) > b_coef_thresh ) then

       ! The equation is a quadratic equation, since a = 0, but b /= 0.
       ! This should very rarely occur for 3-D Luhar.  When it does, the result
       ! will always be two real-valued roots.
       quadratic_roots = quadratic_solve( b_coef, c_coef, d_coef )

       ! Find the maximum root of the two roots.
       max_root = max( real( quadratic_roots(1), kind = core_rknd ), &
                       real( quadratic_roots(2), kind = core_rknd ) )

    else ! |a| = 0 and |b| = 0

       ! The equation is a linear equation.
       ! This won't happen for 3-D Luhar.
       max_root = - d_coef / c_coef

    endif ! |a| > 0


    return

  end function max_cubic_root

!===============================================================================

end module pdf_closure_module
