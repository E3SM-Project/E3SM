! $Id$
!===============================================================================
module LY93_pdf

  ! Description:
  ! The multivariate, two-component PDF of Lewellen and Yoh (1993).

  ! References:
  ! Lewellen, W. S. and Yoh, S., 1993.  Binormal Model of Ensemble Partial
  ! Cloudiness.  J. Atmos. Sci., 50, 9, 1228--1237.
  !-------------------------------------------------------------------------

  implicit none

  public ::  LY93_driver,         & ! Procedure(s)
             calc_mixt_frac_LY93, &
             calc_params_LY93

  private ! default scope
    
  contains

  !=============================================================================
  subroutine LY93_driver( wm, rtm, thlm, wp2, rtp2,          & ! In
                          thlp2, Skw, Skrt, Skthl,           & ! In
                          mu_w_1, mu_w_2, mu_rt_1, mu_rt_2,  & ! Out
                          mu_thl_1, mu_thl_2, sigma_w_1_sqd, & ! Out
                          sigma_w_2_sqd, sigma_rt_1_sqd,     & ! Out
                          sigma_rt_2_sqd, sigma_thl_1_sqd,   & ! Out
                          sigma_thl_2_sqd, mixt_frac         ) ! Out

    ! Description:
    ! Calculates the mixture fraction and the PDF component means and PDF
    ! component variances of w, rt, and theta-l following Lewellen and Yoh.

    ! References:
    ! Lewellen, W. S. and Yoh, S., 1993.  Binormal Model of Ensemble Partial
    ! Cloudiness.  J. Atmos. Sci., 50, 9, 1228--1237.
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Type(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd), dimension(gr%nz), intent(in) :: &
      wm,    & ! Mean of w (overall)                      [m/s]
      wp2,   & ! Variance of w (overall)                  [m^2/s^2]
      Skw,   & ! Skewness of w (overall)                  [-]
      rtm,   & ! Mean of rt (overall)                     [kg/kg]
      rtp2,  & ! Variance of rt (overall)                 [kg^2/kg^2]
      Skrt,  & ! Skewness of rt (overall)                 [-]
      thlm,  & ! Mean of thl (overall)                    [K]
      thlp2, & ! Variance of thl (overall)                [K^2]
      Skthl    ! Skewness of thl (overall)                [-]

    ! Output Variables
    real( kind = core_rknd), dimension(gr%nz), intent(out) :: &
      mu_w_1,          & ! Mean of w (1st PDF component)          [m/s]
      mu_w_2,          & ! Mean of w (2nd PDF component)          [m/s]
      mu_rt_1,         & ! Mean of rt (1st PDF component)         [kg/kg]
      mu_rt_2,         & ! Mean of rt (2nd PDF component)         [kg/kg]
      mu_thl_1,        & ! Mean of thl (1st PDF component)        [K]
      mu_thl_2,        & ! Mean of thl (2nd PDF component)        [K]
      sigma_w_1_sqd,   & ! Variance of w (1st PDF component)      [m^2/s^2]
      sigma_w_2_sqd,   & ! Variance of w (2nd PDF component)      [m^2/s^2]
      sigma_rt_1_sqd,  & ! Variance of rt (1st PDF component)     [m^2/s^2]
      sigma_rt_2_sqd,  & ! Variance of rt (2nd PDF component)     [m^2/s^2]
      sigma_thl_1_sqd, & ! Variance of thl (1st PDF component)    [m^2/s^2]
      sigma_thl_2_sqd, & ! Variance of thl (2nd PDF component)    [m^2/s^2]
      mixt_frac          ! Mixture fraction                       [-]

    ! Local Variables
    real( kind = core_rknd), dimension(gr%nz) :: &
      Sk_max    ! Maximum of magnitudes of skewness        [-]


    ! Find the maximum of the magnitudes of skewness.
    Sk_max = max( abs( Skw ), abs( Skrt ), abs( Skthl ) )

    ! Calculate mixture fraction.
    mixt_frac = calc_mixt_frac_LY93( Sk_max )

    ! Calculate the PDF parameters for w.
    call calc_params_LY93( wm, wp2, Skw, mixt_frac,     & ! In
                           mu_w_1, mu_w_2,              & ! Out
                           sigma_w_1_sqd, sigma_w_2_sqd ) ! Out

    ! Calculate the PDF parameters for rt.
    call calc_params_LY93( rtm, rtp2, Skrt, mixt_frac,    & ! In
                           mu_rt_1, mu_rt_2,              & ! Out
                           sigma_rt_1_sqd, sigma_rt_2_sqd ) ! Out

    ! Calculate the PDF parameters for thl.
    call calc_params_LY93( thlm, thlp2, Skthl, mixt_frac,   & ! In
                           mu_thl_1, mu_thl_2,              & ! Out
                           sigma_thl_1_sqd, sigma_thl_2_sqd ) ! Out


    return

  end subroutine LY93_driver

  !=============================================================================
  function calc_mixt_frac_LY93( Sk_max ) &
  result( mixt_frac )

    ! Description:
    ! Calculates mixture fraction iteratively according to Lewellen and Yoh.

    ! References:
    ! Eq. (21) of Lewellen, W. S. and Yoh, S., 1993.  Binormal Model of Ensemble
    ! Partial Cloudiness.  J. Atmos. Sci., 50, 9, 1228--1237.
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Type(s)

    use constants_clubb, only: &
        one,           & ! Constant(s)
        three_fourths, &
        one_half,      &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variable
    real( kind = core_rknd), dimension(gr%nz), intent(in) :: &
      Sk_max    ! Maximum of magnitudes of skewness        [-]

    ! Return Variable
    real( kind = core_rknd), dimension(gr%nz) :: &
      mixt_frac    ! Mixture fraction                      [-]

    ! Local Variables
    real( kind = core_rknd) :: &
      mixt_frac_low,   & ! Low value of mixture frac. in iterative solver  [-]
      mixt_frac_high,  & ! High value of mixture frac.in iterative solver  [-]
      expr_equal_zero    ! Expr. mixt_frac^6 - Sk_max * ( 1 - mixt_frac )  [-]

    ! Tolerance for mixture fraction in solver [-]
    real( kind = core_rknd) :: &
      LY_mixt_frac_tol = 1.0e-4_core_rknd

    integer :: k    ! Vertical level index


    do k = 1, gr%nz, 1

       if ( Sk_max(k) > 0.84_core_rknd ) then

          mixt_frac_low = one_half
          mixt_frac_high = one

          do ! solve iteratively for mixture fraction

             mixt_frac(k) = one_half * ( mixt_frac_low + mixt_frac_high )

             expr_equal_zero &
             = mixt_frac(k)**6 - Sk_max(k)**2 * ( one - mixt_frac(k) )

             if ( abs( expr_equal_zero ) < LY_mixt_frac_tol ) then
                ! Mixture fraction has been solved for within the specificed
                ! tolerance.
                exit
             else
                if ( expr_equal_zero > zero ) then
                   mixt_frac_high = mixt_frac(k)
                else ! expr_equal_zero < 0
                   mixt_frac_low = mixt_frac(k)
                endif 
             endif

          enddo ! solve iteratively for mixture fraction

       else ! Sk_max <= 0.84

          mixt_frac(k) = three_fourths

       endif

    enddo ! k = 1, gr%nz, 1


    return

  end function calc_mixt_frac_LY93

  !=============================================================================
  subroutine calc_params_LY93( xm, xp2, Skx, mixt_frac,     & ! In
                               mu_x_1, mu_x_2,              & ! Out
                               sigma_x_1_sqd, sigma_x_2_sqd ) ! Out

    ! Description:
    ! Calculates the PDF component means and PDF component variances for
    ! variable x according to Lewellen and Yoh.

    ! References:
    ! Eq. (14), Eq. (15), Eq. (16), Eq. (17), and Eq. (18) of
    ! Lewellen, W. S. and Yoh, S., 1993.  Binormal Model of Ensemble Partial
    ! Cloudiness.  J. Atmos. Sci., 50, 9, 1228--1237.
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Type(s)

    use constants_clubb, only: &
        three,     & ! Constant(s)
        one,       &
        one_third, &
        zero

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    implicit none

    ! Input Variables
    real( kind = core_rknd), dimension(gr%nz), intent(in) :: &
      xm,        & ! Mean of x (overall)        [units vary]
      xp2,       & ! Variance of x (overall)    [(units vary)^2]
      Skx,       & ! Skewness of x (overall)    [-]
      mixt_frac    ! Mixture fraction           [-]

    ! Output Variables
    real( kind = core_rknd), dimension(gr%nz), intent(out) :: &
      mu_x_1,        & ! Mean of x (1st PDF component)        [units vary]
      mu_x_2,        & ! Mean of x (2nd PDF component)        [units vary]
      sigma_x_1_sqd, & ! Variance of x (1st PDF component)    [(units vary)^2]
      sigma_x_2_sqd    ! Variance of x (2nd PDF component)    [(units vary)^2]

    ! Local Variables
    real( kind = core_rknd), dimension(gr%nz) :: &
       sgn_Skx, & ! Sign of Skx                                   [-]
       B_x        ! Spread of the PDF component means function    [units vary]


    ! Find the sign of Skx
    where ( Skx >= zero )
       sgn_Skx = one
    elsewhere ! Skx < 0
       sgn_Skx = -one
    endwhere

    ! Calculate B_x, the LY function for the spread of the PDF component means.
    B_x = sgn_Skx * sqrt( xp2 ) &
          * ( abs( Skx ) / ( one - mixt_frac ) )**one_third

    ! Calculate the mean of x in the 1st PDF component.
    mu_x_1 = xm - B_x * ( one - mixt_frac )

    ! Calculate the mean of x in the 2nd PDF component.
    mu_x_2 = xm + B_x * mixt_frac

    ! Calculate the variance of x in the 1st PDF component.
    sigma_x_1_sqd = xp2 - B_x**2 * ( one - mixt_frac ) &
                          * ( one + mixt_frac + mixt_frac**2 ) &
                          / ( three * mixt_frac )

    ! Calculate the variance of x in the 2nd PDF component.
    sigma_x_2_sqd = xp2 + B_x**2 * ( one - mixt_frac )**2 / three


    return

  end subroutine calc_params_LY93

  !=============================================================================

end module LY93_pdf
