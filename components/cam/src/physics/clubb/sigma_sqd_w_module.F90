!-------------------------------------------------------------------------
! $Id: sigma_sqd_w_module.F90 6849 2014-04-22 21:52:30Z charlass@uwm.edu $
!===============================================================================
module sigma_sqd_w_module

  implicit none

  public :: compute_sigma_sqd_w

  private ! Default scope

  contains
!---------------------------------------------------------------------------------------------------
  elemental function compute_sigma_sqd_w( gamma_Skw_fnc, wp2, thlp2, rtp2, wpthlp, wprtp ) &
    result( sigma_sqd_w )
! Description:
!   Compute the variable sigma_sqd_w (PDF width parameter)
!
! References:
!   Eqn 22 in ``Equations for CLUBB''
!---------------------------------------------------------------------------------------------------
    use constants_clubb, only: &
      w_tol,  & ! Constant(s)
      rt_tol, &
      thl_tol

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: min, max, sqrt

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      gamma_Skw_fnc, & ! Gamma as a function of skewness   [-]
      wp2,           & ! Variance of vertical velocity     [m^2/s^2]
      thlp2,         & ! Variance of liquid pot. temp.     [K^2]
      rtp2,          & ! Variance of total water           [kg^2/kg^2]
      wpthlp,        & ! Flux of liquid pot. temp.         [m/s K]
      wprtp            ! Flux of total water               [m/s kg/kg]

    ! Output Variable
    real( kind = core_rknd ) :: sigma_sqd_w ! PDF width parameter      [-]

    ! ---- Begin Code ----

    !----------------------------------------------------------------
    ! Compute sigma_sqd_w with new formula from Vince
    !----------------------------------------------------------------

    sigma_sqd_w = gamma_Skw_fnc * &
      ( 1.0_core_rknd - min( &
                  max( ( wpthlp / ( sqrt( wp2 * thlp2 )  &
                      + 0.01_core_rknd * w_tol * thl_tol ) )**2, &
                       ( wprtp / ( sqrt( wp2 * rtp2 )  &
                      + 0.01_core_rknd * w_tol * rt_tol ) )**2 &
                     ), & ! max
             1.0_core_rknd ) & ! min - Known magic number (eq. 22 from "Equations for CLUBB")
       )

    return
  end function compute_sigma_sqd_w

end module sigma_sqd_w_module
