! $Id$
!===============================================================================
module new_pdf_main

  ! Description:
  ! The portion of CLUBB's multivariate, two-component PDF that is the
  ! trivariate, two-component normal PDF of vertical velocity (w), total water
  ! mixing ratio (rt), and liquid water potential temperature (thl).

  ! References:
  !-------------------------------------------------------------------------

  implicit none

  public :: new_pdf_driver    ! Procedure(s)

  private :: calc_responder_var,     & ! Procedure(s)
             calc_F_x_zeta_x_setter, &
             calc_F_x_responder

  private

  contains

  !=============================================================================
  subroutine new_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2, Skw,    & ! In
                             wprtp, wpthlp, rtpthlp,                  & ! In
                             Skrt, Skthl,                             & ! In/Out
                             mu_w_1, mu_w_2, mu_rt_1, mu_rt_2,        & ! Out
                             mu_thl_1, mu_thl_2, sigma_w_1_sqd,       & ! Out
                             sigma_w_2_sqd, sigma_rt_1_sqd,           & ! Out
                             sigma_rt_2_sqd, sigma_thl_1_sqd,         & ! Out
                             sigma_thl_2_sqd, mixt_frac,              & ! Out
                             pdf_implicit_coefs_terms,                & ! Out
                             F_w, F_rt, F_thl, min_F_w, max_F_w,      & ! Out
                             min_F_rt, max_F_rt, min_F_thl, max_F_thl ) ! Out
                             

    ! Description:
    ! Selects which variable is used to set the mixture fraction for the PDF
    ! ("the setter") and which variables are handled after that mixture fraction
    ! has been set ("the responders").  Traditionally, w has been used to set
    ! the PDF.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        four,                & ! Variable(s)
        two,                 &
        one,                 &
        zero,                &
        rt_tol,              &
        thl_tol,             &
        max_mag_correlation

    use new_pdf, only: &
        calc_setter_var_params,     & ! Procedure(s)
        calc_coef_wp4_implicit,     &
        calc_coef_wpxp2_implicit,   &
        calc_coefs_wp2xp_semiimpl,  &
        calc_coefs_wpxpyp_semiimpl

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable Type

    use parameters_tunable, only: &
        slope_coef_spread_DG_means_w, & ! Variable(s)
        pdf_component_stdev_factor_w, &
        coef_spread_DG_means_rt, &
        coef_spread_DG_means_thl

    use model_flags, only: &
        l_explicit_turbulent_adv_wp3,  & ! Variable(s)
        l_explicit_turbulent_adv_wpxp, &
        l_explicit_turbulent_adv_xpyp

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wm,      & ! Mean of w (overall)                 [m/s]
      rtm,     & ! Mean of rt (overall)                [kg/kg]
      thlm,    & ! Mean of thl (overall)               [K]
      wp2,     & ! Variance of w (overall)             [m^2/s^2]
      rtp2,    & ! Variance of rt (overall)            [kg^2/kg^2]
      thlp2,   & ! Variance of thl (overall)           [K^2]
      Skw,     & ! Skewness of w (overall)             [-]
      wprtp,   & ! Covariance of w and rt (overall)    [(m/s)kg/kg]
      wpthlp,  & ! Covariance of w and thl (overall)   [(m/s)K]
      rtpthlp    ! Covariance of rt and thl (overall)  [(kg/kg)K]

    ! Input/Output Variables
    ! These variables are input/output because their values may be clipped.
    ! Otherwise, as long as it is not necessary to clip them, their values
    ! will stay the same.
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      Skrt,  & ! Skewness of rt (overall)            [-]
      Skthl    ! Skewness of thl (overall)           [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      mu_w_1,          & ! Mean of w (1st PDF component)        [m/s]
      mu_w_2,          & ! Mean of w (2nd PDF component)        [m/s]
      mu_rt_1,         & ! Mean of rt (1st PDF component)       [kg/kg]
      mu_rt_2,         & ! Mean of rt (2nd PDF component)       [kg/kg]
      mu_thl_1,        & ! Mean of thl (1st PDF component)      [K]
      mu_thl_2,        & ! Mean of thl (2nd PDF component)      [K]
      sigma_w_1_sqd,   & ! Variance of w (1st PDF component)    [m^2/s^2]
      sigma_w_2_sqd,   & ! Variance of w (2nd PDF component)    [m^2/s^2]
      sigma_rt_1_sqd,  & ! Variance of rt (1st PDF component)   [kg^2/kg^2]
      sigma_rt_2_sqd,  & ! Variance of rt (2nd PDF component)   [kg^2/kg^2]
      sigma_thl_1_sqd, & ! Variance of thl (1st PDF component)  [K^2]
      sigma_thl_2_sqd, & ! Variance of thl (2nd PDF component)  [K^2]
      mixt_frac          ! Mixture fraction                     [-]

    type(implicit_coefs_terms), dimension(gr%nz), intent(out) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Output only for recording statistics.
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

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      sigma_w_1,   & ! Standard deviation of w (1st PDF component)      [m/s]
      sigma_w_2,   & ! Standard deviation of w (2nd PDF component)      [m/s]
      sgn_wprtp,   & ! Sign of the covariance of w and rt (overall)     [-]
      sgn_wpthlp,  & ! Sign of the covariance of w and thl (overall)    [-]
      sgn_wp2        ! Sign of the variance of w (overall); always pos. [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      coef_sigma_w_1_sqd,   & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>    [-]
      coef_sigma_w_2_sqd,   & ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>    [-]
      coef_sigma_rt_1_sqd,  & ! sigma_rt_1^2 = coef_sigma_rt_1_sqd * <rt'^2> [-]
      coef_sigma_rt_2_sqd,  & ! sigma_rt_2^2 = coef_sigma_rt_2_sqd * <rt'^2> [-]
      coef_sigma_thl_1_sqd, & ! sigma_thl_1^2=coef_sigma_thl_1_sqd*<thl'^2>  [-]
      coef_sigma_thl_2_sqd    ! sigma_thl_2^2=coef_sigma_thl_2_sqd*<thl'^2>  [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      max_Skx2_pos_Skx_sgn_wpxp, & ! Maximum Skx^2 when Skx*sgn(<w'x'>) >= 0 [-]
      max_Skx2_neg_Skx_sgn_wpxp    ! Maximum Skx^2 when Skx*sgn(<w'x'>) < 0  [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      zeta_w    ! Parameter for the PDF component variances of w           [-]

    real ( kind = core_rknd ) :: &
      lambda_w    ! Param. that increases or decreases Skw dependence  [-]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      exp_factor_rt,   & ! Factor of the form 1 - exp{} that reduces F_rt   [-]
      exp_factor_thl,  & ! Don't reduce F_thl by exp_factor_thl        [-]
      adj_corr_rt_thl    ! Adjusted (overall) correlation of rt and theta-l [-]

    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp4_implicit,     & ! <w'^4> = coef_wp4_implicit * <w'^2>^2       [-]
      coef_wprtp2_implicit,  & ! <w'rt'^2> = coef_wprtp2_implicit*<rt'^2>  [m/s]
      coef_wpthlp2_implicit    ! <w'thl'^2>=coef_wpthlp2_implicit*<thl'^2> [m/s]

    ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2rtp_implicit, & ! Coefficient that is multiplied by <w'rt'>  [m/s]
      term_wp2rtp_explicit    ! Term that is on the RHS          [m^2/s^2 kg/kg]

    ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wp2thlp_implicit, & ! Coef. that is multiplied by <w'thl'>      [m/s]
      term_wp2thlp_explicit    ! Term that is on the RHS             [m^2/s^2 K]

    ! <w'rt'thl'> = coef_wprtpthlp_implicit*<rt'thl'> + term_wprtpthlp_explicit
    real ( kind = core_rknd ), dimension(gr%nz) :: &
      coef_wprtpthlp_implicit, & ! Coef. that is multiplied by <rt'thl'>   [m/s]
      term_wprtpthlp_explicit    ! Term that is on the RHS         [m/s(kg/kg)K]


    ! Calculate sgn( <w'rt'> ).
    where ( wprtp >= zero )
       sgn_wprtp = one
    elsewhere ! wprtp < 0
       sgn_wprtp = -one
    endwhere ! wprtp >= 0

    ! Calculate sgn( <w'thl'> ).
    where ( wpthlp >= zero )
       sgn_wpthlp = one
    elsewhere ! wpthlp < 0
       sgn_wpthlp = -one
    endwhere ! wpthlp >= 0

    ! Sign of the variance of w (overall), which is always positive.
    sgn_wp2 = one

    lambda_w = 0.5_core_rknd

    ! Calculate the adjusted (overall) correlation of rt and theta-l, and the
    ! value of exp_factor_rt.
    where ( rtp2 >= rt_tol**2 .and. thlp2 >= thl_tol**2 )
       adj_corr_rt_thl = rtpthlp / sqrt( rtp2 * thlp2 ) * sgn_wprtp * sgn_wpthlp
       adj_corr_rt_thl = min( max( adj_corr_rt_thl, -max_mag_correlation ), &
                              max_mag_correlation )
       exp_factor_rt = one &
                       - exp( -0.2_core_rknd * ( adj_corr_rt_thl + one )**5 )
    elsewhere ! <rt'^2> < rt_tol^2 or <thl'^2> < thl_tol^2
       adj_corr_rt_thl = zero  ! adj_corr_rt_thl is undefined in this scenario.
       exp_factor_rt = one     ! Set exp_factor_rt to 1.
    endwhere ! <rt'^2> >= rt_tol^2 and <thl'^2> >= thl_tol^2

    ! The value of F_thl is not reduced by exp_factor_thl.
    exp_factor_thl = one


    ! Vertical velocity, w, will always be the setter variable.
    call calc_F_x_zeta_x_setter( Skw,                          & ! In
                                 slope_coef_spread_DG_means_w, & ! In
                                 pdf_component_stdev_factor_w, & ! In
                                 lambda_w,                     & ! In
                                 F_w, zeta_w,                  & ! Out
                                 min_F_w, max_F_w              ) ! Out

    ! Calculate the PDF parameters, including mixture fraction, for the
    ! setter variable, w.
    call calc_setter_var_params( wm, wp2, Skw, sgn_wp2,     & ! In
                                 F_w, zeta_w,               & ! In
                                 mu_w_1, mu_w_2, sigma_w_1, & ! Out
                                 sigma_w_2, mixt_frac,      & ! Out
                                 coef_sigma_w_1_sqd,        & ! Out
                                 coef_sigma_w_2_sqd         ) ! Out

    sigma_w_1_sqd = sigma_w_1**2
    sigma_w_2_sqd = sigma_w_2**2

    ! Calculate the upper limit on the magnitude of skewness for responding
    ! variables.
    max_Skx2_pos_Skx_sgn_wpxp = four * ( one - mixt_frac )**2 &
                                / ( mixt_frac * ( two - mixt_frac ) )

    max_Skx2_neg_Skx_sgn_wpxp = four * mixt_frac**2 / ( one - mixt_frac**2 )

    ! Calculate the PDF parameters for responder variable rt.
    call calc_responder_var( rtm, rtp2, sgn_wprtp, mixt_frac, & ! In
                             coef_spread_DG_means_rt,         & ! In
                             exp_factor_rt,                   & ! In
                             max_Skx2_pos_Skx_sgn_wpxp,       & ! In
                             max_Skx2_neg_Skx_sgn_wpxp,       & ! In
                             Skrt,                            & ! In/Out
                             mu_rt_1, mu_rt_2,                & ! Out
                             sigma_rt_1_sqd, sigma_rt_2_sqd,  & ! Out
                             coef_sigma_rt_1_sqd,             & ! Out
                             coef_sigma_rt_2_sqd,             & ! Out
                             F_rt, min_F_rt, max_F_rt         ) ! Out

    ! Calculate the PDF parameters for responder variable thl.
    call calc_responder_var( thlm, thlp2, sgn_wpthlp, mixt_frac, & ! In
                             coef_spread_DG_means_thl,           & ! In
                             exp_factor_thl,                     & ! In
                             max_Skx2_pos_Skx_sgn_wpxp,          & ! In
                             max_Skx2_neg_Skx_sgn_wpxp,          & ! In
                             Skthl,                              & ! In/Out
                             mu_thl_1, mu_thl_2,                 & ! Out
                             sigma_thl_1_sqd, sigma_thl_2_sqd,   & ! Out
                             coef_sigma_thl_1_sqd,               & ! Out
                             coef_sigma_thl_2_sqd,               & ! Out
                             F_thl, min_F_thl, max_F_thl         ) ! Out


    if ( .not. l_explicit_turbulent_adv_wp3 ) then

       ! Turbulent advection of <w'^3> is being handled implicitly.

       ! <w'^4> = coef_wp4_implicit * <w'^2>^2.
       coef_wp4_implicit &
       = calc_coef_wp4_implicit( mixt_frac, F_w, &
                                 coef_sigma_w_1_sqd, &
                                 coef_sigma_w_2_sqd )

    else ! l_explicit_turbulent_adv_wp3

       ! Turbulent advection of <w'^3> is being handled explicitly.
       coef_wp4_implicit = zero

    endif ! .not. l_explicit_turbulent_adv_wp3

    if ( .not. l_explicit_turbulent_adv_xpyp ) then

       ! Turbulent advection of <rt'^2> and <thl'^2> is being handled
       ! implicitly.  Turbulent advection of <rt'thl'> is being handled
       ! semi-implicitly.

       ! <w'rt'^2> = coef_wprtp2_implicit * <rt'^2>
       coef_wprtp2_implicit &
       = calc_coef_wpxp2_implicit( wp2, rtp2, wprtp, sgn_wprtp, &
                                   mixt_frac, F_w, F_rt, &
                                   coef_sigma_w_1_sqd, &
                                   coef_sigma_w_2_sqd, &
                                   coef_sigma_rt_1_sqd, &
                                   coef_sigma_rt_2_sqd  )

       ! <w'thl'^2> = coef_wpthlp2_implicit * <thl'^2>
       coef_wpthlp2_implicit &
       = calc_coef_wpxp2_implicit( wp2, thlp2, wpthlp, sgn_wpthlp, &
                                   mixt_frac, F_w, F_thl, &
                                   coef_sigma_w_1_sqd, &
                                   coef_sigma_w_2_sqd, &
                                   coef_sigma_thl_1_sqd, &
                                   coef_sigma_thl_2_sqd  )

       ! <w'rt'thl'> = coef_wprtpthlp_implicit * <rt'thl'>
       !               + term_wprtpthlp_explicit
       call calc_coefs_wpxpyp_semiimpl( wp2, rtp2, thlp2, wprtp,       & ! In
                                        wpthlp, sgn_wprtp, sgn_wpthlp, & ! In
                                        mixt_frac, F_w, F_rt, F_thl,   & ! In
                                        coef_sigma_w_1_sqd  ,          & ! In
                                        coef_sigma_w_2_sqd,            & ! In
                                        coef_sigma_rt_1_sqd,           & ! In
                                        coef_sigma_rt_2_sqd,           & ! In
                                        coef_sigma_thl_1_sqd,          & ! In
                                        coef_sigma_thl_2_sqd,          & ! In
                                        coef_wprtpthlp_implicit,       & ! Out
                                        term_wprtpthlp_explicit        ) ! Out

    else ! l_explicit_turbulent_adv_xpyp

       ! Turbulent advection of <rt'^2>, <thl'^2>, and <rt'thl'> is being
       ! handled explicitly.
       coef_wprtp2_implicit = zero
       coef_wpthlp2_implicit = zero
       coef_wprtpthlp_implicit = zero
       term_wprtpthlp_explicit = zero

    endif ! .not. l_explicit_turbulent_adv_xpyp

    if ( .not. l_explicit_turbulent_adv_wpxp ) then

       ! Turbulent advection of <w'rt'> and <w'thl'> is being handled
       ! semi-implicitly.

       ! <w'^2 rt'> = coef_wp2rtp_implicit * <w'rt'> + term_wp2rtp_explicit
       call calc_coefs_wp2xp_semiimpl( wp2, rtp2, sgn_wprtp, & ! In
                                       mixt_frac, F_w, F_rt, & ! In
                                       coef_sigma_w_1_sqd,   & ! In
                                       coef_sigma_w_2_sqd,   & ! In
                                       coef_sigma_rt_1_sqd,  & ! In
                                       coef_sigma_rt_2_sqd,  & ! In
                                       coef_wp2rtp_implicit, & ! Out
                                       term_wp2rtp_explicit  ) ! Out

       ! <w'^2 thl'> = coef_wp2thlp_implicit * <w'thl'> + term_wp2thlp_explicit
       call calc_coefs_wp2xp_semiimpl( wp2, thlp2, sgn_wpthlp, & ! In
                                       mixt_frac, F_w, F_thl,  & ! In
                                       coef_sigma_w_1_sqd,     & ! In
                                       coef_sigma_w_2_sqd,     & ! In
                                       coef_sigma_thl_1_sqd,   & ! In
                                       coef_sigma_thl_2_sqd,   & ! In
                                       coef_wp2thlp_implicit,  & ! Out
                                       term_wp2thlp_explicit   ) ! Out

    else ! l_explicit_turbulent_adv_wpxp

       ! Turbulent advection of <w'rt'> and <w'thl'> is being handled
       ! explicitly.
       coef_wp2rtp_implicit = zero
       term_wp2rtp_explicit = zero
       coef_wp2thlp_implicit = zero
       term_wp2thlp_explicit = zero

    endif ! .not. l_explicit_turbulent_adv_wpxp

    ! Pack the implicit coefficients and explicit terms into a single type
    ! variable for output.
    pdf_implicit_coefs_terms%coef_wp4_implicit = coef_wp4_implicit
    pdf_implicit_coefs_terms%coef_wprtp2_implicit = coef_wprtp2_implicit
    pdf_implicit_coefs_terms%coef_wpthlp2_implicit = coef_wpthlp2_implicit
    pdf_implicit_coefs_terms%coef_wp2rtp_implicit = coef_wp2rtp_implicit
    pdf_implicit_coefs_terms%term_wp2rtp_explicit = term_wp2rtp_explicit
    pdf_implicit_coefs_terms%coef_wp2thlp_implicit = coef_wp2thlp_implicit
    pdf_implicit_coefs_terms%term_wp2thlp_explicit = term_wp2thlp_explicit
    pdf_implicit_coefs_terms%coef_wprtpthlp_implicit = coef_wprtpthlp_implicit
    pdf_implicit_coefs_terms%term_wprtpthlp_explicit = term_wprtpthlp_explicit


    return

  end subroutine new_pdf_driver

  !=============================================================================
  subroutine calc_responder_var( xm, xp2, sgn_wpxp, mixt_frac, & ! In
                                 coef_spread_DG_means_x,       & ! In
                                 exp_factor_x,                 & ! In
                                 max_Skx2_pos_Skx_sgn_wpxp,    & ! In
                                 max_Skx2_neg_Skx_sgn_wpxp,    & ! In
                                 Skx,                          & ! In/Out
                                 mu_x_1, mu_x_2,               & ! Out
                                 sigma_x_1_sqd, sigma_x_2_sqd, & ! Out
                                 coef_sigma_x_1_sqd,           & ! Out
                                 coef_sigma_x_2_sqd,           & ! Out
                                 F_x, min_F_x, max_F_x         ) ! Out

    ! Description:
    ! This is the sub-driver for a responder variable.  The upper limits of the
    ! magnitude of Skx are calculated, and Skx is clipped when its magnitude
    ! exceeds the upper limits.  The limits of the F_x parameter are calculated,
    ! and the value of F_x is set within those limits.  Then, the PDF parameters
    ! for responder variable x are calculated.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        zero    ! Variable(s)

    use new_pdf, only: &
        calc_limits_F_x_responder, & ! Procedure(s)
        calc_responder_params

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xm,           & ! Mean of x (overall)                       [units vary]
      xp2,          & ! Variance of x (overall)               [(units vary)^2]
      sgn_wpxp,     & ! Sign of the covariance of w and x                  [-]
      mixt_frac,    & ! Mixture fraction                                   [-]
      exp_factor_x    ! Factor of the form 1 - exp{}; reduces F_x          [-]

    real( kind = core_rknd ), intent(in) :: &
      coef_spread_DG_means_x    ! Coef.: spread betw. PDF comp. means of x   [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      max_Skx2_pos_Skx_sgn_wpxp, & ! Maximum Skx^2 when Skx*sgn(<w'x'>) >= 0 [-]
      max_Skx2_neg_Skx_sgn_wpxp    ! Maximum Skx^2 when Skx*sgn(<w'x'>) < 0  [-]

    ! Input/Output Variable
    real( kind = core_rknd ), dimension(gr%nz), intent(inout) :: &
      Skx    ! Skewness of x (overall)              [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)        [units vary]
      mu_x_2,        & ! Mean of x (2nd PDF component)        [units vary]
      sigma_x_1_sqd, & ! Variance of x (1st PDF component)    [(units vary)^2]
      sigma_x_2_sqd    ! Variance of x (2nd PDF component)    [(units vary)^2]

    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>    [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>    [-]

    ! Output only for recording statistics.
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      F_x,     & ! Param. for the spread betw. the PDF component means of x  [-]
      min_F_x, & ! Minimum allowable value of parameter F_x                  [-]
      max_F_x    ! Maximum allowable value of parameter F_x                  [-]


    ! Calculate the upper limit of the magnitude of Skx.
    where ( Skx * sgn_wpxp >= zero )
       where ( Skx**2 >= max_Skx2_pos_Skx_sgn_wpxp )
          where ( Skx >= zero )
             Skx = sqrt( 0.99_core_rknd * max_Skx2_pos_Skx_sgn_wpxp )
          elsewhere
             Skx = -sqrt( 0.99_core_rknd * max_Skx2_pos_Skx_sgn_wpxp )
          endwhere
       endwhere ! Skx^2 >= max_Skx2_pos_Skx_sgn_wpxp
    elsewhere ! Skx * sgn( <w'x'> ) < 0
       where ( Skx**2 >= max_Skx2_neg_Skx_sgn_wpxp )
          where ( Skx >= zero )
             Skx = sqrt( 0.99_core_rknd * max_Skx2_neg_Skx_sgn_wpxp )
          elsewhere
             Skx = -sqrt( 0.99_core_rknd * max_Skx2_neg_Skx_sgn_wpxp )
          endwhere
       endwhere ! Skx^2 >= max_Skx2_neg_Skx_sgn_wpxp
    endwhere ! Skx * sgn( <w'x'> ) >= 0

    call calc_limits_F_x_responder( mixt_frac, Skx, sgn_wpxp,  & ! In
                                    max_Skx2_pos_Skx_sgn_wpxp, & ! In
                                    max_Skx2_neg_Skx_sgn_wpxp, & ! In
                                    min_F_x, max_F_x )           ! Out

    ! F_x must have a value between min_F_x and max_F_x.
    F_x = calc_F_x_responder( coef_spread_DG_means_x, exp_factor_x, &
                              min_F_x, max_F_x )

    call calc_responder_params( xm, xp2, Skx, sgn_wpxp,       & ! In
                                F_x, mixt_frac,               & ! In
                                mu_x_1, mu_x_2,               & ! Out
                                sigma_x_1_sqd, sigma_x_2_sqd, & ! Out
                                coef_sigma_x_1_sqd,           & ! Out
                                coef_sigma_x_2_sqd            ) ! Out


    return

  end subroutine calc_responder_var

  !=============================================================================
  subroutine calc_F_x_zeta_x_setter( Skx,                          & ! In
                                     slope_coef_spread_DG_means_x, & ! In
                                     pdf_component_stdev_factor_x, & ! In
                                     lambda,                       & ! In
                                     F_x, zeta_x,                  & ! Out
                                     min_F_x, max_F_x              ) ! Out

    ! Description:
    ! Calculates the values of F_x and zeta_x for the setter variable (which is
    ! the variable that sets the mixture fraction).
    !
    ! The value of F_x is calculated between 0 (min_F_x) and 1 (max_F_x).  The
    ! equation is:
    !
    ! F_x = max_F_x + ( min_F_x - max_F_x )
    !                 * exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x };
    !
    ! which reduces to:
    !
    ! F_x = 1 - exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x };
    !
    ! where lambda > 0 and slope_coef_spread_DG_means_x > 0.  As |Skx| goes
    ! toward 0, the value of F_x goes toward 0, and as |Skx| becomes large, the
    ! value of F_x goes toward 1.  When slope_coef_spread_DG_means_x is small,
    ! the value of F_x tends toward 1, and when slope_coef_spread_DG_means_x is
    ! large, the value of F_x tends toward 0.  When lambda is small, the value
    ! of F_x is less dependent on Skx, and when lambda is large, the value of
    ! F_x is more dependent on Skx.
    !
    ! Mathematically, this equation will always produce a value of F_x that
    ! falls between min_F_x and max_F_x.  However, in order to prevent a value
    ! of F_x from being calculated outside the bounds of min_F_x and max_F_x
    ! owing to numerical underflow or loss of precision, this equation can be
    ! rewritten as:
    !
    ! F_x
    ! = min_F_x * exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x }
    !   + max_F_x * ( 1 - exp{ -|Skx|^lambda / slope_coef_spread_DG_means_x } ).
    !
    ! The value of zeta_x used to adjust the PDF component standard devations:
    !
    ! 1 + zeta_x = ( mixt_frac * sigma_x_1^2 )
    !              / ( ( 1 - mixt_frac ) * sigma_x_2^2 );
    !
    ! where zeta_x > -1.  The sign of zeta_x is used to easily determine if
    ! mixt_frac * sigma_x_1^2 is greater than ( 1 - mixt_frac ) * sigma_x_2^2
    ! (when zeta_x is positive), mixt_frac * sigma_x_1^2 is less than
    ! ( 1 - mixt_frac ) * sigma_x_2^2 (when zeta_x is negative), or
    ! mixt_frac * sigma_x_1^2 is equal to ( 1 - mixt_frac ) * sigma_x_2^2 (when
    ! zeta_x is 0).
    !
    ! In order to allow for a tunable parameter that is the pure ratio of
    ! mixt_frac * sigma_x_1^2 to ( 1 - mixt_frac ) * sigma_x_2^2, zeta_x is
    ! related to the parameter pdf_component_stdev_factor_x, where:
    !
    ! 1 + zeta_x = pdf_component_stdev_factor_x.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        one,  & ! Variable(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      Skx    ! Skewness of x (overall)              [-]

    real( kind = core_rknd ), intent(in) :: &
      slope_coef_spread_DG_means_x, & ! Slope coef: spread PDF comp. means x [-]
      pdf_component_stdev_factor_x, & ! Param.: PDF comp. standard devs.; x  [-]
      lambda                          ! Param. for Skx dependence            [-]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(out) :: &
      F_x,     & ! Parameter for the spread of the PDF component means of x  [-]
      zeta_x,  & ! Parameter for the PDF component variances of x            [-]
      min_F_x, & ! Minimum allowable value of parameter F_x                  [-]
      max_F_x    ! Maximum allowable value of parameter F_x                  [-]

    ! Local Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      exp_Skx_interp_factor    ! Function to interp. between min. and max.   [-]


    ! Set min_F_x to 0 and max_F_x to 1 for the setter variable.
    where ( abs( Skx ) > zero )
       min_F_x = 1.0e-3_core_rknd
    elsewhere
       min_F_x = zero
    endwhere
    max_F_x = one

    ! F_x must have a value between min_F_x and max_F_x.
    exp_Skx_interp_factor &
    = exp( -abs(Skx)**lambda / slope_coef_spread_DG_means_x )

    F_x = min_F_x * exp_Skx_interp_factor &
          + max_F_x * ( one - exp_Skx_interp_factor )

    ! The value of zeta_x must be greater than -1.
    zeta_x = pdf_component_stdev_factor_x - one


    return

  end subroutine calc_F_x_zeta_x_setter

  !=============================================================================
  function calc_F_x_responder( coef_spread_DG_means_x, exp_factor_x, &
                               min_F_x, max_F_x ) &
  result( F_x )

    ! Description:
    ! Calculates the value of F_x as a tunable function between min_F_x and
    ! max_F_x.
    !
    ! The value of F_x is calculated between min_F_x and max_F_x.  The equation
    ! is:
    !
    ! F_x = min_F_x + ( max_F_x - min_F_x )
    !                 * coef_spread_DG_means_x * exp_factor_x;
    !
    ! where 0 <= coef_spread_DG_means_x <= 1.  As coef_spread_DG_means_x
    ! goes toward 0, the value of F_x goes toward min_F_x, and as
    ! coef_spread_DG_means_x goes toward 1, the value of F_x goes toward
    ! max_F_x.  The exp_factor_x is a factor of the form 1 - exp{ }.  The range
    ! of values of exp_factor_x is 0 <= exp_factor_x <= 1.  Here, exp_factor_x
    ! is used to reduce the value of F_x under special conditions.
    !
    ! The equation for exp_factor_x depends on which responder variable is being
    ! solved for.  For rt, the F_rt equation is:
    !
    ! F_rt = min_F_rt + ( max_F_rt - min_F_rt )
    !                   * coef_spread_DG_means_rt * exp_factor_rt;
    !
    ! where exp_factor_rt is given by:
    !
    ! exp_factor_rt = 1 - exp{ -0.2 * ( adj_corr_rt_thl + 1 )^5 };
    !
    ! and where the adjusted (overall) correlation of rt and theta-l, denoted
    ! adj_corr_rt_thl, is given by:
    !
    ! adj_corr_rt_thl = <rt'thl'> / ( sqrt( <rt'^2> ) * sqrt( <thl'^2> ) )
    !                   * sgn( <w'rt'> ) * sgn( <w'thl'> ).
    !
    ! The range of values of the adjusted (overall) correlation of rt and
    ! theta-l is -1 <= adj_corr_rt_thl <= 1.  The adjusted (overall) correlation
    ! of rt and theta-l has the same absolute magnitude as the (overall)
    ! correlation of rt and theta-l, but the sign of the adjusted correlation is
    ! dependent on the agreement between the signs of <w'rt'> and <w'thl'>.
    ! When all three covariances are in sign agreement, which is when:
    !
    ! sgn( <w'rt'> ) * sgn( <w'thl'> ) * sgn( <rt'thl'> ) = 1,
    !
    ! the value of adj_corr_rt_thl is positive.  However, when the three
    ! covariances aren't consistent in sign, which is when:
    !
    ! sgn( <w'rt'> ) * sgn( <w'thl'> ) * sgn( <rt'thl'> ) = -1,
    !
    ! the value of adj_corr_rt_thl is negative.
    !
    ! For theta-l, the F_thl equation is:
    !
    ! F_thl = min_F_thl + ( max_F_thl - min_F_thl ) * coef_spread_DG_means_thl;
    !
    ! where exp_factor_thl is set to 1 (1 - exp{-inf} = 1) because reducing
    ! F_thl through the use of exp_factor_thl is not desired for theta-l.
    !
    ! The direction of the two PDF component means of x (mu_x_1 and mu_x_2) with
    ! respect to each other and the overall mean of x (<x>) is given by
    ! sgn( <w'x'> ).  When sgn( <w'x'> ) = 1, ( mu_x_1 - <x> ) >= 0 and
    ! ( mu_x_2 - <x> ) <= 0.  When sgn( <w'x'> ) = -1, ( mu_x_1 - <x> ) <= 0 and
    ! ( mu_x_2 - <x> ) >= 0.  This helps to promote realizability of the PDF
    ! component correlations (corr_w_x_1 and corr_w_x_2).
    !
    ! The realizability of the PDF component correlations corr_w_rt_1,
    ! corr_w_rt_2, corr_w_thl_1, and corr_w_thl_2 is promoted by using
    ! sgn( <w'rt'> ) and sgn( <w'thl'> ), respectively to choose the direction
    ! of the PDF component means of rt and theta-l.  However, when the three
    ! covariances (<w'rt'>, <w'thl'>, and <rt'thl'>) aren't consistent in sign
    ! (as shown above), it becomes difficult to keep the PDF component
    ! correlations corr_rt_thl_1 and corr_rt_thl_2 in the realizable range.
    ! However, in this situation, the above equation for exp_factor_rt brings
    ! F_rt closer to min_F_rt, which promotes realizability for corr_rt_thl_1
    ! and corr_rt_thl_2.
    !
    ! Mathematically, this equation will always produce a value of F_x that
    ! falls between min_F_x and max_F_x.  However, in order to prevent a value
    ! of F_x from being calculated outside the bounds of min_F_x and max_F_x
    ! owing to numerical underflow or loss of precision, this equation can be
    ! rewritten as:
    !
    ! F_x = min_F_x * ( 1 - coef_spread_DG_means_x * exp_factor_x )
    !       + max_F_x * coef_spread_DG_means_x * exp_factor_x.

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
    real( kind = core_rknd ), intent(in) :: &
      coef_spread_DG_means_x    ! Coef.: spread betw. PDF comp. means of x   [-]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      exp_factor_x, & ! Factor of the form 1 - exp{}; reduces F_x  [-]
      min_F_x,      & ! Minimum allowable value of parameter F_x   [-]
      max_F_x         ! Maximum allowable value of parameter F_x   [-]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      F_x    ! Parameter for the spread between the PDF component means of x [-]


    ! F_x must have a value between min_F_x and max_F_x.
    F_x = min_F_x * ( one - coef_spread_DG_means_x * exp_factor_x ) &
          + max_F_x * coef_spread_DG_means_x * exp_factor_x


    return

  end function calc_F_x_responder

  !=============================================================================

end module new_pdf_main
