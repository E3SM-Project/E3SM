!-------------------------------------------------------------------------
! $Id$
!===============================================================================
module sigma_sqd_w_module

  implicit none

  public :: compute_sigma_sqd_w

  private ! Default scope

  contains

  !=============================================================================
  elemental function compute_sigma_sqd_w( gamma_Skw_fnc, wp2, thlp2, rtp2, &
                                          up2, vp2, wpthlp, wprtp, upwp, vpwp ) &
    result( sigma_sqd_w )

    ! Description:
    ! Compute the variable sigma_sqd_w (PDF width parameter).
    !
    ! The value of sigma_sqd_w is restricted in the ADG1 PDF in order to keep
    ! the marginal PDFs of all responder variables (variables that do not set
    ! the mixture fraction) valid.  The limits on sigma_sqd_w in order to keep
    ! the PDF of a responder variable, x, valid are:
    !
    ! 0 <= sigma_sqd_w <= 1 - <w'x'>^2 / ( <w'^2> * <x'^2> ).
    !
    ! The overall limits on sigma_sqd_w must be applied based on the most
    ! restrictive case so that all Double Gaussian PDF responder variables, x,
    ! have realizable PDFs.  The overall limits on sigma_sqd_w are:
    !
    ! 0 <= sigma_sqd_w <= 1 - max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all x).
    !
    ! The equation used for sigma_sqd_w is:
    !
    ! sigma_sqd_w = gamma_Skw_fnc
    !               * ( 1 - max( <w'x'>^2 / ( <w'^2> * <x'^2> ), for all x) );
    !
    ! where 0 <= gamma_Skw_fnc <= 1.

    ! References:
    !   Eqn 22 in ``Equations for CLUBB''
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,       & ! Constant(s)
        w_tol,     &
        rt_tol,    &
        thl_tol,   &
        w_tol_sqd

    use model_flags, only: &
        l_predict_upwp_vpwp    ! Variable(s)

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: min, max, sqrt

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      gamma_Skw_fnc, & ! Gamma as a function of skewness             [-]
      wp2,           & ! Variance of vertical velocity               [m^2/s^2]
      thlp2,         & ! Variance of liquid water potential temp.    [K^2]
      rtp2,          & ! Variance of total water mixing ratio        [kg^2/kg^2]
      up2,           & ! Variance of west-east horizontal velocity   [m^2/s^2]
      vp2,           & ! Variance of south-north horizontal velocity [m^2/s^2]
      wpthlp,        & ! Flux of liquid water potential temp.        [m/s K]
      wprtp,         & ! Flux of total water mixing ratio            [m/s kg/kg]
      upwp,          & ! Flux of west-east horizontal velocity       [m^2/s^2]
      vpwp             ! Flux of south-north horizontal velocity     [m^2/s^2]

    ! Output Variable
    real( kind = core_rknd ) :: sigma_sqd_w ! PDF width parameter      [-]

    ! Local Variable
    real( kind = core_rknd ) :: &
      max_corr_w_x_sqd    ! Max. val. of wpxp^2/(wp2*xp2) for all vars. x  [-]

    ! ---- Begin Code ----

    !----------------------------------------------------------------
    ! Compute sigma_sqd_w with new formula from Vince
    !----------------------------------------------------------------

    ! Find the maximum value of <w'x'>^2 / ( <w'^2> * <x'^2> ) for all
    ! variables x that are Double Gaussian PDF responder variables.  This
    ! includes rt and theta-l.  When l_predict_upwp_vpwp is enabled, u and v are
    ! also calculated as part of the PDF, and they are included as well.
    ! Additionally, when sclr_dim > 0, passive scalars (sclr) are also included.
    max_corr_w_x_sqd = max( ( wpthlp / ( sqrt( wp2 * thlp2 ) &
                            + 0.01_core_rknd * w_tol * thl_tol ) )**2, &
                            ( wprtp / ( sqrt( wp2 * rtp2 )  &
                            + 0.01_core_rknd * w_tol * rt_tol ) )**2 )

    if ( l_predict_upwp_vpwp ) then
       max_corr_w_x_sqd = max( max_corr_w_x_sqd, &
                               ( upwp / ( sqrt( up2 * wp2 ) &
                               + 0.01_core_rknd * w_tol_sqd ) )**2, &
                               ( vpwp / ( sqrt( vp2 * wp2 ) &
                               + 0.01_core_rknd * w_tol_sqd ) )**2 )
    endif ! l_predict_upwp_vpwp

    ! Calculate the value of sigma_sqd_w .
    sigma_sqd_w = gamma_Skw_fnc * ( one - min( max_corr_w_x_sqd, one ) )


    return

  end function compute_sigma_sqd_w

  !=============================================================================

end module sigma_sqd_w_module
