! $Id$
!===============================================================================
module new_tsdadg_pdf

  ! Description:
  ! The new trivariate, skewness-dependent, analytic double Gaussian (TSDADG)
  ! PDF.

  implicit none

  public :: tsdadg_pdf_driver,      & ! Procedure(s)
            calc_setter_parameters, &
            calc_L_x_Skx_fnc

  private :: calc_respnder_parameters    ! Procedure(s)

  private  ! default scope

  contains

  !=============================================================================
  subroutine tsdadg_pdf_driver( wm, rtm, thlm, wp2, rtp2, thlp2,   & ! In
                                Skw, Skrt, Skthl, wprtp, wpthlp,   & ! In
                                mu_w_1, mu_w_2, mu_rt_1, mu_rt_2,  & ! Out
                                mu_thl_1, mu_thl_2, sigma_w_1_sqd, & ! Out
                                sigma_w_2_sqd, sigma_rt_1_sqd,     & ! Out
                                sigma_rt_2_sqd, sigma_thl_1_sqd,   & ! Out
                                sigma_thl_2_sqd, mixt_frac         ) ! Out


    ! Description:
    ! Selects which variable is used to set the mixture fraction for the PDF
    ! ("the setter") and which variables are handled after that mixture fraction
    ! has been set ("the responders").  Traditionally, w has been used to set
    ! the PDF.  However, here, the variable with the greatest magnitude of
    ! skewness is used to set the PDF.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable type(s)

    use constants_clubb, only: &
        one,     & ! Variable(s)
        zero,    &
        fstderr

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      wm,     & ! Mean of w (overall)                [m/s]
      rtm,    & ! Mean of rt (overall)               [kg/kg]
      thlm,   & ! Mean of thl (overall)              [K]
      wp2,    & ! Variance of w (overall)            [m^2/s^2]
      rtp2,   & ! Variance of rt (overall)           [kg^2/kg^2]
      thlp2,  & ! Variance of thl (overall)          [K^2]
      Skw,    & ! Skewness of w (overall)            [-]
      Skrt,   & ! Skewness of rt (overall)           [-]
      Skthl,  & ! Skewness of thl (overall)          [-]
      wprtp,  & ! Covariance of w and rt (overall)   [(m/s)kg/kg]
      wpthlp    ! Covariance of w and thl (overall)  [(m/s)K]

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

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nz) :: &
      big_L_w_1,   & ! Parameter for the 1st PDF comp. mean of w            [-]
      big_L_w_2,   & ! Parameter for the 2nd PDF comp. mean of w (setter)   [-]
      big_L_rt_1,  & ! Parameter for the 1st PDF comp. mean of rt           [-]
      big_L_rt_2,  & ! Parameter for the 2nd PDF comp. mean of rt (setter)  [-]
      big_L_thl_1, & ! Parameter for the 1st PDF comp. mean of thl          [-]
      big_L_thl_2    ! Parameter for the 2nd PDF comp. mean of thl (setter) [-]

    real( kind = core_rknd ) :: &
      small_l_w_1,   & ! Param. for the 1st PDF comp. mean of w            [-]
      small_l_w_2,   & ! Param. for the 2nd PDF comp. mean of w (setter)   [-]
      small_l_rt_1,  & ! Param. for the 1st PDF comp. mean of rt           [-]
      small_l_rt_2,  & ! Param. for the 2nd PDF comp. mean of rt (setter)  [-]
      small_l_thl_1, & ! Param. for the 1st PDF comp. mean of thl          [-]
      small_l_thl_2    ! Param. for the 2nd PDF comp. mean of thl (setter) [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      sgn_wprtp,  & ! Sign of the covariance of w and rt (overall)         [-]
      sgn_wpthlp, & ! Sign of the covariance of w and thl (overall)        [-]
      sgn_wp2       ! Sign of the variance of w (overall); always positive [-]

    real( kind = core_rknd ), dimension(gr%nz) :: &
      coef_sigma_w_1_sqd,   & ! sigma_w_1^2 = coef_sigma_w_1_sqd * <w'^2>    [-]
      coef_sigma_w_2_sqd,   & ! sigma_w_2^2 = coef_sigma_w_2_sqd * <w'^2>    [-]
      coef_sigma_rt_1_sqd,  & ! sigma_rt_1^2 = coef_sigma_rt_1_sqd * <rt'^2> [-]
      coef_sigma_rt_2_sqd,  & ! sigma_rt_2^2 = coef_sigma_rt_2_sqd * <rt'^2> [-]
      coef_sigma_thl_1_sqd, & ! sigma_thl_1^2=coef_sigma_thl_1_sqd*<thl'^2>  [-]
      coef_sigma_thl_2_sqd    ! sigma_thl_2^2=coef_sigma_thl_2_sqd*<thl'^2>  [-]

    integer :: k    ! Vertical level loop index


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

    ! The sign of the variance of w is always positive.
    sgn_wp2 = one

    small_l_w_1 = 0.75_core_rknd
    small_l_w_2 = 0.5_core_rknd
    small_l_rt_1 = 0.75_core_rknd
    small_l_rt_2 = 0.5_core_rknd
    small_l_thl_1 = 0.75_core_rknd
    small_l_thl_2 = 0.5_core_rknd

    do k = 1, gr%nz, 1

       call calc_L_x_Skx_fnc( Skw(k), sgn_wp2(k),        & ! In
                              small_l_w_1, small_l_w_2,  & ! In
                              big_L_w_1(k), big_L_w_2(k) ) ! Out

       call calc_L_x_Skx_fnc( Skrt(k), sgn_wprtp(k),       & ! In
                              small_l_rt_1, small_l_rt_2,  & ! In
                              big_L_rt_1(k), big_L_rt_2(k) ) ! Out

       call calc_L_x_Skx_fnc( Skthl(k), sgn_wpthlp(k),       & ! In
                              small_l_thl_1, small_l_thl_2,  & ! In
                              big_L_thl_1(k), big_L_thl_2(k) ) ! Out


       ! The variable with the greatest magnitude of skewness will be the setter
       ! variable and the other variables will be responder variables.
       if ( abs( Skw(k) ) >= abs( Skrt(k) ) &
            .and. abs( Skw(k) ) >= abs( Skthl(k) ) ) then

          ! The variable w has the greatest magnitude of skewness.

          call calc_setter_parameters( wm(k), wp2(k),                  & ! In
                                       Skw(k), sgn_wp2(k),             & ! In
                                       big_L_w_1(k), big_L_w_2(k),     & ! In
                                       mu_w_1(k), mu_w_2(k),           & ! Out
                                       sigma_w_1_sqd(k),               & ! Out
                                       sigma_w_2_sqd(k), mixt_frac(k), & ! Out
                                       coef_sigma_w_1_sqd(k),          & ! Out
                                       coef_sigma_w_2_sqd(k)           ) ! Out

          call calc_respnder_parameters( rtm(k), rtp2(k),             & ! In
                                         Skrt(k), sgn_wprtp(k),       & ! In
                                         mixt_frac(k), big_L_rt_1(k), & ! In
                                         mu_rt_1(k), mu_rt_2(k),      & ! Out
                                         sigma_rt_1_sqd(k),           & ! Out
                                         sigma_rt_2_sqd(k),           & ! Out
                                         coef_sigma_rt_1_sqd(k),      & ! Out
                                         coef_sigma_rt_2_sqd(k)       ) ! Out

          call calc_respnder_parameters( thlm(k), thlp2(k),            & ! In
                                         Skthl(k), sgn_wpthlp(k),      & ! In
                                         mixt_frac(k), big_L_thl_1(k), & ! In
                                         mu_thl_1(k), mu_thl_2(k),     & ! Out
                                         sigma_thl_1_sqd(k),           & ! Out
                                         sigma_thl_2_sqd(k),           & ! Out
                                         coef_sigma_thl_1_sqd(k),      & ! Out
                                         coef_sigma_thl_2_sqd(k)       ) ! Out


       elseif ( abs( Skrt(k) ) > abs( Skw(k) ) &
                .and. abs( Skrt(k) ) >= abs( Skthl(k) ) ) then

          ! The variable rt has the greatest magnitude of skewness.

          call calc_setter_parameters( rtm(k), rtp2(k),                 & ! In
                                       Skrt(k), sgn_wprtp(k),           & ! In
                                       big_L_rt_1(k), big_L_rt_2(k),    & ! In
                                       mu_rt_1(k), mu_rt_2(k),          & ! Out
                                       sigma_rt_1_sqd(k),               & ! Out
                                       sigma_rt_2_sqd(k), mixt_frac(k), & ! Out
                                       coef_sigma_rt_1_sqd(k),          & ! Out
                                       coef_sigma_rt_2_sqd(k)           ) ! Out

          call calc_respnder_parameters( wm(k), wp2(k),              & ! In
                                         Skw(k), sgn_wp2(k),         & ! In
                                         mixt_frac(k), big_L_w_1(k), & ! In
                                         mu_w_1(k), mu_w_2(k),       & ! Out
                                         sigma_w_1_sqd(k),           & ! Out
                                         sigma_w_2_sqd(k),           & ! Out
                                         coef_sigma_w_1_sqd(k),      & ! Out
                                         coef_sigma_w_2_sqd(k)       ) ! Out

          call calc_respnder_parameters( thlm(k), thlp2(k),            & ! In
                                         Skthl(k), sgn_wpthlp(k),      & ! In
                                         mixt_frac(k), big_L_thl_1(k), & ! In
                                         mu_thl_1(k), mu_thl_2(k),     & ! Out
                                         sigma_thl_1_sqd(k),           & ! Out
                                         sigma_thl_2_sqd(k),           & ! Out
                                         coef_sigma_thl_1_sqd(k),      & ! Out
                                         coef_sigma_thl_2_sqd(k)       ) ! Out


       else ! abs( Skthl ) > abs( Skw ) .and. abs( Skthl ) > abs( Skrt )

          ! The variable thl has the greatest magnitude of skewness.

          call calc_setter_parameters( thlm(k), thlp2(k),                & ! In
                                       Skthl(k), sgn_wpthlp(k),          & ! In
                                       big_L_thl_1(k), big_L_thl_2(k),   & ! In
                                       mu_thl_1(k), mu_thl_2(k),         & ! Out
                                       sigma_thl_1_sqd(k),               & ! Out
                                       sigma_thl_2_sqd(k), mixt_frac(k), & ! Out
                                       coef_sigma_thl_1_sqd(k),          & ! Out
                                       coef_sigma_thl_2_sqd(k)           ) ! Out

          call calc_respnder_parameters( wm(k), wp2(k),              & ! In
                                         Skw(k), sgn_wp2(k),         & ! In
                                         mixt_frac(k), big_L_w_1(k), & ! In
                                         mu_w_1(k), mu_w_2(k),       & ! Out
                                         sigma_w_1_sqd(k),           & ! Out
                                         sigma_w_2_sqd(k),           & ! Out
                                         coef_sigma_w_1_sqd(k),      & ! Out
                                         coef_sigma_w_2_sqd(k)       ) ! Out

          call calc_respnder_parameters( rtm(k), rtp2(k),             & ! In
                                         Skrt(k), sgn_wprtp(k),       & ! In
                                         mixt_frac(k), big_L_rt_1(k), & ! In
                                         mu_rt_1(k), mu_rt_2(k),      & ! Out
                                         sigma_rt_1_sqd(k),           & ! Out
                                         sigma_rt_2_sqd(k),           & ! Out
                                         coef_sigma_rt_1_sqd(k),      & ! Out
                                         coef_sigma_rt_2_sqd(k)       ) ! Out


       endif ! Find variable with the greatest magnitude of skewness.


       if ( sigma_w_1_sqd(k) < zero ) then
          write(fstderr,*) "WARNING:  New TSDADG PDF.  The variance of w in " &
                           // "the 1st PDF component is negative and is " &
                           // "being clipped to 0."
          write(fstderr,*) "sigma_w_1^2 (before clipping) = ", sigma_w_1_sqd(k)
          sigma_w_1_sqd(k) = zero
          coef_sigma_w_1_sqd(k) = zero
       endif ! sigma_w_1_sqd < 0

       if ( sigma_w_2_sqd(k) < zero ) then
          write(fstderr,*) "WARNING:  New TSDADG PDF.  The variance of w in " &
                           // "the 2nd PDF component is negative and is " &
                           // "being clipped to 0."
          write(fstderr,*) "sigma_w_2^2 (before clipping) = ", sigma_w_2_sqd(k)
          sigma_w_2_sqd(k) = zero
          coef_sigma_w_2_sqd(k) = zero
       endif ! sigma_w_2_sqd < 0

       if ( sigma_rt_1_sqd(k) < zero ) then
          write(fstderr,*) "WARNING:  New TSDADG PDF.  The variance of rt in " &
                           // "the 1st PDF component is negative and is " &
                           // "being clipped to 0."
          write(fstderr,*) "sigma_rt_1^2 (before clipping) = ", &
                           sigma_rt_1_sqd(k)
          sigma_rt_1_sqd(k) = zero
          coef_sigma_rt_1_sqd(k) = zero
       endif ! sigma_rt_1_sqd < 0

       if ( sigma_rt_2_sqd(k) < zero ) then
          write(fstderr,*) "WARNING:  New TSDADG PDF.  The variance of rt in " &
                           // "the 2nd PDF component is negative and is " &
                           // "being clipped to 0."
          write(fstderr,*) "sigma_rt_2^2 (before clipping) = ", &
                           sigma_rt_2_sqd(k)
          sigma_rt_2_sqd(k) = zero
          coef_sigma_rt_2_sqd(k) = zero
       endif ! sigma_rt_2_sqd < 0

       if ( sigma_thl_1_sqd(k) < zero ) then
          write(fstderr,*) "WARNING:  New TSDADG PDF.  The variance of thl " &
                           // "in the 1st PDF component is negative and is " &
                           // "being clipped to 0."
          write(fstderr,*) "sigma_thl_1^2 (before clipping) = ", &
                           sigma_thl_1_sqd(k)
          sigma_thl_1_sqd(k) = zero
          coef_sigma_thl_1_sqd(k) = zero
       endif ! sigma_thl_1_sqd < 0

       if ( sigma_thl_2_sqd(k) < zero ) then
          write(fstderr,*) "WARNING:  New TSDADG PDF.  The variance of thl " &
                           // "in the 2nd PDF component is negative and is " &
                           // "being clipped to 0."
          write(fstderr,*) "sigma_thl_2^2 (before clipping) = ", &
                           sigma_thl_2_sqd(k)
          sigma_thl_2_sqd(k) = zero
          coef_sigma_thl_2_sqd(k) = zero
       endif ! sigma_thl_2_sqd < 0

    enddo ! k = 1, gr%nz, 1


    return

  end subroutine tsdadg_pdf_driver

  !=============================================================================
  !
  ! DESCRIPTION OF THE METHOD FOR THE VARIABLE THAT SETS THE MIXTURE FRACTION
  ! =========================================================================
  !
  ! There are five PDF parameters that need to be calculated for the setting
  ! variable, which are mu_x_1 (the mean of x in the 1st PDF component), mu_x_2
  ! (the mean of x in the 2nd PDF component), sigma_x_1 (the standard deviation
  ! of x in the 1st PDF component), sigma_x_2 (the standard deviation of x in
  ! the 2nd PDF component), and mixt_frac (the mixture fraction).  In order to
  ! solve for these five parameters, five equations are needed.  These five
  ! equations are the equations for <x>, <x'^2>, and <x'^3> as found by
  ! integrating over the PDF.  Additionally, two more equations, which involve
  ! tunable parameters L_x_1 and L_x_2, and which are used to help control the
  ! spread of the PDF component means in the 1st PDF component and the 2nd PDF
  ! component, respectively, are used in this equation set.  The five equations
  ! are:
  !
  ! <x> = mixt_frac * mu_x_1 + ( 1 - mixt_frac ) * mu_x_2;
  !
  ! <x'^2> = mixt_frac * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
  !
  ! <x'^3> = mixt_frac * ( mu_x_1 - <x> )
  !                    * ( ( mu_x_1 - <x> )^2 + 3 * sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )
  !                              * ( ( mu_x_2 - <x> )^2 + 3 * sigma_x_2^2 );
  !
  ! mu_x_1 - <x> = L_x_1
  !                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> ); and
  !
  ! mu_x_2 - <x> = -L_x_2
  !                 * sqrt( ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                         / ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                 * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! where 0 <= L_x_1 <= 1, 0 <= L_x_2 <= 1, Skx is the skewness of x, such that
  ! Skx = <x'^3> / <x'^2>^(3/2), and sgn( <w'x'> ) is the sign of <w'x'>, such
  ! that:
  !
  ! sgn( <w'x'> ) = |  1, when <w'x'> >= 0;
  !                 | -1, when <w'x'> < 0.
  !
  ! The resulting equations for the five PDF parameters are:
  !
  ! mu_x_1 = <x> + L_x_1
  !                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! mu_x_2 = <x> - L_x_2
  !                * sqrt( ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                        / ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! mixt_frac = 1 / ( 1 + abs( mu_x_1_nrmlized / mu_x_2_nrmlized ) );
  !
  ! sigma_x_1 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
  !                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                     + ( 1 - mixt_frac )
  !                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                           - mu_x_1_nrmlized^2 / 3
  !                           + mu_x_2_nrmlized^2 / 3 ) )
  !                   * <x'^2> ); and
  !
  ! sigma_x_2 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
  !                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                     - mixt_frac
  !                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                           - mu_x_1_nrmlized^2 / 3
  !                           + mu_x_2_nrmlized^2 / 3 ) )
  !                   * <x'^2> ); where
  !
  ! mu_x_1_nrmlized = ( mu_x_1 - <x> ) / sqrt( <x'^2> ); and
  !
  ! mu_x_2_nrmlized = ( mu_x_2 - <x> ) / sqrt( <x'^2> ).
  !
  !
  ! Notes:
  !
  ! This method does NOT work for all values of L_x_1 and L_x_2 (where
  ! 0 <= L_x_1 <= 1 and 0 <= L_x_2 <= 1).  Only a subregion of this parameter
  ! space produces valid results.
  !
  ! When both L_x_1 = 0 and L_x_2 = 0, mu_x_1 = mu_x_2 = <x> (which can only
  ! happen when Skx = 0).  In this scenario, the above equations for mixt_frac,
  ! sigma_x_1, and sigma_x_2 are all undefined.  In this special case, the
  ! distribution reduces to a single Gaussian, so the following values are set:
  ! mixt_frac = 1/2, sigma_x_1 = sqrt( <x'^2> ), and sigma_x_2 = sqrt( <x'^2> ).
  !
  !
  ! Tunable parameters:
  !
  ! The parameter L_x_1 controls the 1st PDF component mean while L_x_2 controls
  ! the 2nd PDF component mean.  The equations involving the tunable parameters
  ! L_x_1 and L_x_2 (the mu_x_1 and mu_x_2 equations) are based on the values of
  ! mu_x_1 and mu_x_2 when sigma_x_1 = sigma_x_2 = 0.  In this scenario, the
  ! equation for mixture fraction reduces to:
  !
  ! mixt_frac = (1/2) * ( 1 +/- Skx / sqrt( 4 + Skx^2 ) ).
  !
  ! The +/- is dependent on the sign of ( mu_x_1 - <x> ) vs. ( mu_x_2 - <x> ).
  ! This is dependent on sgn( <w'x'> ), and the mixture fraction equation is
  ! written as:
  !
  ! mixt_frac = (1/2) * ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ).
  !
  ! Meanwhile, the equation for 1 - mixt_frac is:
  !
  ! 1 - mixt_frac = (1/2) * ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ).
  !
  ! When sigma_x_1 = sigma_x_2 = 0, the equations for mu_x_1 and mu_x_2 are:
  !
  ! mu_x_1 = <x> + sqrt( ( 1 - mixt_frac ) / mixt_frac ) * sqrt( <x'^2> )
  !                * sgn( <w'x'> ); and
  !
  ! mu_x_2 = <x> - sqrt( mixt_frac / ( 1 - mixt_frac ) ) * sqrt( <x'^2> )
  !                * sgn( <w'x'> ).
  !
  ! Substituting the equations for mixt_frac and 1 - mixt_frac into the
  ! equations for mu_x_1 and mu_x_2 (when sigma_x_1 = sigma_x_2 = 0), the
  ! equations for mu_x_1 and mu_x_2 become:
  !
  ! mu_x_1 = <x> + sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                      / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> ); and
  !
  ! mu_x_2 = <x> - sqrt( ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                      / ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> ).
  !
  ! These equations represent the maximum deviation of mu_x_1 and mu_x_2 from
  ! the overall mean, <x>.  The range of parameters of L_x_i is 0 <= L_x_i <= 1.
  ! When L_x_1 = L_x_2 = 0, the value of mu_x_1 = mu_x_2 = <x> (and the
  ! distribution becomes a single Gaussian).  When L_x_i = 1, the value of
  ! mu_x_i - <x> is at its maximum magnitude.
  !
  ! The values of L_x_1 and L_x_2 are also calculated by skewness functions.
  ! Those functions are:
  !
  ! L_x_1 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
  ! L_x_2 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 );
  !
  ! where both l_x_1 and l_x_2 are tunable parameters.
  !
  ! As previously stated, this method does not work for all combinations of
  ! L_x_1 and L_x_2, but rather only for a subregion of parameter space.  This
  ! applies to l_x_1 and l_x_2, as well.  The conditions on l_x_1 and l_x_2 are:
  !
  ! 2/3 < l_x_1 < 1 and 0 < l_x_2 < 1; when Skx * sgn( <w'x'> ) >= 0; and
  ! 0 < l_x_1 < 1 and 2/3 < l_x_2 < 1; when Skx * sgn( <w'x'> ) < 0.
  !
  ! The condition that l_x_1 > 2/3 prevents a negative PDF component variance
  ! when Skx = 0.
  !
  !
  ! Equations for PDF component standard deviations:
  !
  ! The equations for the PDF component standard deviations can also be written
  ! as:
  !
  ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
  !
  ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ); where
  !
  ! coef_sigma_x_1_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
  !                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                      + ( 1 - mixt_frac )
  !                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                            - mu_x_1_nrmlized^2 / 3
  !                            + mu_x_2_nrmlized^2 / 3 ); and
  !
  ! coef_sigma_x_2_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
  !                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                      - mixt_frac
  !                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                            - mu_x_1_nrmlized^2 / 3
  !                            + mu_x_2_nrmlized^2 / 3 ).
  !
  ! The above equations can be substituted into an equation for a variable that
  ! has been derived by integrating over the PDF.  Many variables like this are
  ! used in parts of the predictive equation set.  These substitutions allow
  ! some terms to solved implicitly or semi-implicitly in the predictive
  ! equations.
  !
  !
  !=============================================================================
  subroutine calc_setter_parameters( xm, xp2, Skx, sgn_wpxp,        & ! In
                                     big_L_x_1, big_L_x_2,          & ! In
                                     mu_x_1, mu_x_2, sigma_x_1_sqd, & ! Out
                                     sigma_x_2_sqd, mixt_frac,      & ! Out
                                     coef_sigma_x_1_sqd,            & ! Out
                                     coef_sigma_x_2_sqd             ) ! Out

    ! Description:
    ! Calculates the PDF component means, the PDF component standard deviations,
    ! and the mixture fraction for the variable that sets the PDF.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        four,     & ! Variable(s)
        three,    &
        one,      &
        one_half, &
        zero,     &
        eps

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,        & ! Mean of x (overall)                            [units vary]
      xp2,       & ! Variance of x (overall)                    [(units vary)^2]
      Skx,       & ! Skewness of x                                           [-]
      sgn_wpxp,  & ! Sign of the covariance of w and x (overall)             [-]
      big_L_x_1, & ! Parameter for the spread of the 1st PDF comp. mean of x [-]
      big_L_x_2    ! Parameter for the spread of the 2nd PDF comp. mean of x [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)              [units vary]
      mu_x_2,        & ! Mean of x (2nd PDF component)              [units vary]
      sigma_x_1_sqd, & ! Variance of x (1st PDF component)      [(units vary)^2]
      sigma_x_2_sqd, & ! Variance of x (2nd PDF component)      [(units vary)^2]
      mixt_frac        ! Mixture fraction                b                   [-]

    real( kind = core_rknd ), intent(out) :: &
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      mu_x_1_nrmlized, & ! Normalized mean of x (1st PDF component)          [-]
      mu_x_2_nrmlized    ! Normalized mean of x (2nd PDF component)          [-]

    real( kind = core_rknd ) :: &
      factor_plus,               &
      factor_minus,              &
      sqrt_factor_plus_ov_minus, &
      sqrt_factor_minus_ov_plus, &
      mu_x_1_nrmlized_thresh


    ! Calculate the factors in the PDF component mean equations.
    factor_plus = one + Skx * sgn_wpxp / sqrt( four + Skx**2 )

    factor_minus = one - Skx * sgn_wpxp / sqrt( four + Skx**2 )

    sqrt_factor_plus_ov_minus = sqrt( factor_plus / factor_minus )
    sqrt_factor_minus_ov_plus = sqrt( factor_minus / factor_plus )

    ! Calculate the normalized mean of x in the 1st PDF component.
    mu_x_1_nrmlized = big_L_x_1 * sqrt_factor_plus_ov_minus * sgn_wpxp

    ! Calculate the normalized mean of x in the 2nd PDF component.
    mu_x_2_nrmlized = -big_L_x_2 * sqrt_factor_minus_ov_plus * sgn_wpxp

    ! Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm + mu_x_1_nrmlized * sqrt( xp2 )

    ! Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm + mu_x_2_nrmlized * sqrt( xp2 )

    ! Calculate the mixture fraction.
    if ( abs( mu_x_1_nrmlized ) >= eps &
         .and. abs( mu_x_2_nrmlized ) >= eps ) then
       mixt_frac = one / ( one + abs( mu_x_1_nrmlized / mu_x_2_nrmlized ) )
    elseif ( abs( mu_x_1_nrmlized ) >= eps &
             .and. abs( mu_x_2_nrmlized ) < eps ) then
       mixt_frac = one / ( one + abs( mu_x_1_nrmlized / eps ) )
    elseif ( abs( mu_x_1_nrmlized ) < eps &
             .and. abs( mu_x_2_nrmlized ) >= eps ) then
       mixt_frac = one / ( one + abs( eps / mu_x_2_nrmlized ) )
    else ! abs( mu_x_1_nrmlized ) < eps and abs( mu_x_2_nrmlized ) < eps
       mixt_frac = one_half
    endif

    ! Use a minimum magnitude value of mu_x_1_nrmlized in the denominator of a
    ! term in the PDF component variance equations in order to prevent a
    ! divide-by-zero error.
    if ( mu_x_1_nrmlized >= zero ) then
       mu_x_1_nrmlized_thresh = max( mu_x_1_nrmlized, eps )
    else ! mu_x_1_nrmlized < 0
       mu_x_1_nrmlized_thresh = min( mu_x_1_nrmlized, -eps )
    endif ! mu_x_1_nrmlized >= 0

    ! Calculate the variance of x in the 1st PDF component.
    coef_sigma_x_1_sqd &
    = one - mixt_frac * mu_x_1_nrmlized**2 &
      - ( one - mixt_frac ) * mu_x_2_nrmlized**2 &
      + ( one - mixt_frac ) &
        * ( Skx / ( three * mixt_frac * mu_x_1_nrmlized_thresh ) &
            - mu_x_1_nrmlized**2 / three + mu_x_2_nrmlized**2 / three )

    sigma_x_1_sqd = coef_sigma_x_1_sqd * xp2

    ! Calculate the variance of x in the 2nd PDF component.
    coef_sigma_x_2_sqd & 
    = one - mixt_frac * mu_x_1_nrmlized**2 &
      - ( one - mixt_frac ) * mu_x_2_nrmlized**2 &
      - mixt_frac &
        * ( Skx / ( three * mixt_frac * mu_x_1_nrmlized_thresh ) &
            - mu_x_1_nrmlized**2 / three + mu_x_2_nrmlized**2 / three )

    sigma_x_2_sqd = coef_sigma_x_2_sqd * xp2


    return

  end subroutine calc_setter_parameters

  !=============================================================================
  subroutine calc_L_x_Skx_fnc( Skx, sgn_wpxp,            & ! In
                               small_l_x_1, small_l_x_2, & ! In
                               big_L_x_1, big_L_x_2      ) ! Out

    ! Description:
    ! Calculates the values of big_L_x_1 and big_L_x_2 as functions of Skx.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        four, & ! Variable(s)
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skx,         & ! Skewness of x (overall)                               [-]
      sgn_wpxp,    & ! Sign of the covariance of w and x (overall)           [-]
      small_l_x_1, & ! Param. for the spread of the 1st PDF comp. mean of x  [-]
      small_l_x_2    ! Param. for the spread of the 2nd PDF comp. mean of x  [-]

    ! Output Variable
    real( kind = core_rknd ), intent(out) :: &
      big_L_x_1, & ! Parameter for the spread of the 1st PDF comp. mean of x [-]
      big_L_x_2    ! Parameter for the spread of the 2nd PDF comp. mean of x [-]

    ! Local Variable
    real( kind = core_rknd ) :: &
      Skx_fnc_factor


    ! The values of L_x_1 and L_x_2 are calculated by skewness functions.
    ! Those functions are:
    !
    ! L_x_1 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
    ! L_x_2 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 ).
    !
    ! The conditions on l_x_1 and l_x_2 are:
    !
    ! 2/3 < l_x_1 < 1 and 0 < l_x_2 < 1; when Skx * sgn( <w'x'> ) >= 0; and
    ! 0 < l_x_1 < 1 and 2/3 < l_x_2 < 1; when Skx * sgn( <w'x'> ) < 0.
    !
    ! For simplicity, this can also be accomplished by setting 2/3 < l_x_1 < 1
    ! and 0 < l_x_2 < 1, and then using the following equations.
    !
    ! When Skx * sgn( <w'x'> ) >= 0:
    ! L_x_1 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
    ! L_x_2 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 );
    !
    ! otherwise, when Skx * sgn( <w'x'> ) < 0, switch l_x_1 and l_x_2:
    ! L_x_1 = l_x_2 * abs( Skx ) / sqrt( 4 + Skx^2 ); and
    ! L_x_2 = l_x_1 * abs( Skx ) / sqrt( 4 + Skx^2 ).
    Skx_fnc_factor = abs( Skx ) / sqrt( four + Skx**2 )

    if ( Skx * sgn_wpxp >= zero ) then
       big_L_x_1 = small_l_x_1 * Skx_fnc_factor
       big_L_x_2 = small_l_x_2 * Skx_fnc_factor
    else ! Skx * sgn_wpxp < 0
       big_L_x_1 = small_l_x_2 * Skx_fnc_factor
       big_L_x_2 = small_l_x_1 * Skx_fnc_factor
    endif


    return

  end subroutine calc_L_x_Skx_fnc

  !=============================================================================
  !
  ! DESCRIPTION OF THE METHOD FOR EACH RESPONDING VARIABLE
  ! ======================================================
  !
  ! In order to find equations for the four PDF parameters for each responding
  ! variable, which are mu_x_1, mu_x_2, sigma_x_1, and sigma_x_2 (where x stands
  ! for a responding variable here), four equations are needed.  These four
  ! equations are the equations for <x>, <x'^2>, and <x'^3> as found by
  ! integrating over the PDF.  Additionally, one more equation, which involves
  ! tunable parameter L_x_1, and which is used to help control the mean of the
  ! 1st PDF component, is used in this equation set.  The four equations are:
  !
  ! <x> = mixt_frac * mu_x_1 + ( 1 - mixt_frac ) * mu_x_2;
  !
  ! <x'^2> = mixt_frac * ( ( mu_x_1 - <x> )^2 + sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( ( mu_x_2 - <x> )^2 + sigma_x_2^2 );
  !
  ! <x'^3> = mixt_frac * ( mu_x_1 - <x> )
  !                    * ( ( mu_x_1 - <x> )^2 + 3 * sigma_x_1^2 )
  !          + ( 1 - mixt_frac ) * ( mu_x_2 - <x> )
  !                              * ( ( mu_x_2 - <x> )^2 + 3 * sigma_x_2^2 ); and
  !
  ! mu_x_1 - <x> = L_x_1
  !                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! where 0 <= L_x_1 <= 1, Skx is the skewness of x, such that
  ! Skx = <x'^3> / <x'^2>^(3/2), and sgn( <w'x'> ) is given by:
  !
  ! sgn( <w'x'> ) = |  1, when <w'x'> >= 0;
  !                 | -1, when <w'x'> < 0.
  !
  ! The resulting equations for the four PDF parameters are:
  !
  ! mu_x_1 = <x> + L_x_1
  !                * sqrt( ( 1 + Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) )
  !                        / ( 1 - Skx * sgn( <w'x'> ) / sqrt( 4 + Skx^2 ) ) )
  !                * sqrt( <x'^2> ) * sgn( <w'x'> );
  !
  ! mu_x_2 = <x> - ( mixt_frac / ( 1 - mixt_frac ) ) * ( mu_x_1 - <x> );
  !
  ! sigma_x_1 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
  !                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                     + ( 1 - mixt_frac )
  !                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                           - mu_x_1_nrmlized^2 / 3
  !                           + mu_x_2_nrmlized^2 / 3 ) )
  !                   * <x'^2> ); and
  !
  ! sigma_x_2 = sqrt( ( 1 - mixt_frac * mu_x_1_nrmlized^2
  !                     - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                     - mixt_frac
  !                       * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                           - mu_x_1_nrmlized^2 / 3
  !                           + mu_x_2_nrmlized^2 / 3 ) )
  !                   * <x'^2> ); where
  !
  ! mu_x_1_nrmlized = ( mu_x_1 - <x> ) / sqrt( <x'^2> ); and
  !
  ! mu_x_2_nrmlized = ( mu_x_2 - <x> ) / sqrt( <x'^2> ).
  !
  !
  ! Notes:
  !
  ! When L_x_1 = 0, mu_x_1 = mu_x_2 = <x>, sigma_x_1^2 = sigma_x_2^2 = <x'^2>,
  ! and the distribution reduces to a single Gaussian.
  !
  !
  ! Equations for PDF component standard deviations:
  !
  ! The equations for the PDF component standard deviations can also be written
  ! as:
  !
  ! sigma_x_1 = sqrt( coef_sigma_x_1_sqd * <x'^2> ); and
  !
  ! sigma_x_2 = sqrt( coef_sigma_x_2_sqd * <x'^2> ); where
  !
  ! coef_sigma_x_1_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
  !                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                      + ( 1 - mixt_frac )
  !                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                            - mu_x_1_nrmlized^2 / 3
  !                            + mu_x_2_nrmlized^2 / 3 ); and
  !
  ! coef_sigma_x_2_sqd = 1 - mixt_frac * mu_x_1_nrmlized^2
  !                      - ( 1 - mixt_frac ) * mu_x_2_nrmlized^2
  !                      - mixt_frac
  !                        * ( Skx / ( 3 * mixt_frac * mu_x_1_nrmlized )
  !                            - mu_x_1_nrmlized^2 / 3
  !                            + mu_x_2_nrmlized^2 / 3 ).
  !
  ! The above equations can be substituted into an equation for a variable that
  ! has been derived by integrating over the PDF.  Many variables like this are
  ! used in parts of the predictive equation set.  These substitutions allow
  ! some terms to solved implicitly or semi-implicitly in the predictive
  ! equations.
  !
  !
  !=============================================================================
  subroutine calc_respnder_parameters( xm, xp2, Skx, sgn_wpxp,       & ! In
                                       mixt_frac, big_L_x_1,         & ! In
                                       mu_x_1, mu_x_2,               & ! Out
                                       sigma_x_1_sqd, sigma_x_2_sqd, & ! Out
                                       coef_sigma_x_1_sqd,           & ! Out
                                       coef_sigma_x_2_sqd            ) ! Out

    ! Description:
    ! Calculates the PDF component means, the PDF component standard deviations,
    ! and the mixture fraction for the variable that sets the PDF.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        four,  & ! Variable(s)
        three, &
        one,   &
        zero,  &
        eps

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      xm,        & ! Mean of x (overall)                            [units vary]
      xp2,       & ! Variance of x (overall)                    [(units vary)^2]
      Skx,       & ! Skewness of x                                           [-]
      sgn_wpxp,  & ! Sign of the covariance of w and x (overall)             [-]
      mixt_frac, & ! Mixture fraction                                        [-]
      big_L_x_1    ! Parameter for the spread of the 1st PDF comp. mean of x [-]

    ! Output Variables
    real( kind = core_rknd ), intent(out) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)              [units vary]
      mu_x_2,        & ! Mean of x (2nd PDF component)              [units vary]
      sigma_x_1_sqd, & ! Variance of x (1st PDF component)      [(units vary)^2]
      sigma_x_2_sqd    ! Variance of x (2nd PDF component)      [(units vary)^2]

    real( kind = core_rknd ), intent(out) :: &
      coef_sigma_x_1_sqd, & ! sigma_x_1^2 = coef_sigma_x_1_sqd * <x'^2>      [-]
      coef_sigma_x_2_sqd    ! sigma_x_2^2 = coef_sigma_x_2_sqd * <x'^2>      [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      mu_x_1_nrmlized, & ! Normalized mean of x (1st PDF component)          [-]
      mu_x_2_nrmlized    ! Normalized mean of x (2nd PDF component)          [-]

    real( kind = core_rknd ) :: &
      factor_plus,               &
      factor_minus,              &
      sqrt_factor_plus_ov_minus, &
      mu_x_1_nrmlized_thresh


    ! Calculate the factors in the PDF component mean equations.
    factor_plus = one + Skx * sgn_wpxp / sqrt( four + Skx**2 )

    factor_minus = one - Skx * sgn_wpxp / sqrt( four + Skx**2 )

    sqrt_factor_plus_ov_minus = sqrt( factor_plus / factor_minus )

    ! Calculate the normalized mean of x in the 1st PDF component.
    mu_x_1_nrmlized = big_L_x_1 * sqrt_factor_plus_ov_minus * sgn_wpxp

    ! Calculate the normalized mean of x in the 2nd PDF component.
    mu_x_2_nrmlized = - ( mixt_frac / ( one - mixt_frac ) ) * mu_x_1_nrmlized

    ! Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm + mu_x_1_nrmlized * sqrt( xp2 )

    ! Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm + mu_x_2_nrmlized * sqrt( xp2 )

    ! Use a minimum magnitude value of mu_x_1_nrmlized in the denominator of a
    ! term in the PDF component variance equations in order to prevent a
    ! divide-by-zero error.
    if ( mu_x_1_nrmlized >= zero ) then
       mu_x_1_nrmlized_thresh = max( mu_x_1_nrmlized, eps )
    else ! mu_x_1_nrmlized < 0
       mu_x_1_nrmlized_thresh = min( mu_x_1_nrmlized, -eps )
    endif ! mu_x_1_nrmlized >= 0

    ! Calculate the variance of x in the 1st PDF component.
    coef_sigma_x_1_sqd &
    = one - mixt_frac * mu_x_1_nrmlized**2 &
      - ( one - mixt_frac ) * mu_x_2_nrmlized**2 &
      + ( one - mixt_frac ) &
        * ( Skx / ( three * mixt_frac * mu_x_1_nrmlized_thresh ) &
            - mu_x_1_nrmlized**2 / three + mu_x_2_nrmlized**2 / three )

    sigma_x_1_sqd = coef_sigma_x_1_sqd * xp2

    ! Calculate the variance of x in the 2nd PDF component.
    coef_sigma_x_2_sqd & 
    = one - mixt_frac * mu_x_1_nrmlized**2 &
      - ( one - mixt_frac ) * mu_x_2_nrmlized**2 &
      - mixt_frac &
        * ( Skx / ( three * mixt_frac * mu_x_1_nrmlized_thresh ) &
            - mu_x_1_nrmlized**2 / three + mu_x_2_nrmlized**2 / three )

    sigma_x_2_sqd = coef_sigma_x_2_sqd * xp2


    return

  end subroutine calc_respnder_parameters

  !=============================================================================

end module new_tsdadg_pdf
