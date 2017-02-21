!-------------------------------------------------------------------------
!$Id: Skx_module.F90 7810 2015-07-09 13:34:25Z weberjk@uwm.edu $
!===============================================================================
module Skx_module

  implicit none

  private ! Default Scope

  public :: Skx_func, LG_2005_ansatz

  contains

!-------------------------------------------------------------------------------
  elemental function Skx_func( xp2, xp3, x_tol )  &
    result( Skx )

! Description:
!   Calculate the skewness of x

! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only:  &
      Skw_max_mag ! Max magnitude of skewness

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    use parameters_tunable, only: &
      Skw_denom_coef

    implicit none

    ! External
    intrinsic :: min, max

    ! Parameter Constants
    ! Whether to apply clipping to the final result
    logical, parameter ::  & 
      l_clipping_kluge = .false.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      xp2,  & ! x'^2
      xp3,  & ! x'^3
      x_tol   ! x tolerance value

    ! Output Variable
    real( kind = core_rknd ) :: & 
      Skx     ! Result Skw [-]

    ! ---- Begin Code ----

    !Skw = xp3 / ( max( xp2, x_tol**two ) )**1.5_core_rknd
    ! Calculation of skewness to help reduce the sensitivity of this value to
    ! small values of xp2.
    Skx = xp3 / ( xp2 + Skw_denom_coef * x_tol**2 )**1.5_core_rknd

    ! This is no longer needed since clipping is already
    ! imposed on wp2 and wp3 elsewhere in the code

    ! I turned clipping on in this local copy since thlp3 and rtp3 are not clipped
    if ( l_clipping_kluge ) then
      Skx = min( max( Skx, -Skw_max_mag ), Skw_max_mag )
    end if

    return
  end function Skx_func
!-----------------------------------------------------------------------

!-------------------------------------------------------------------------------
  elemental function LG_2005_ansatz( Skw, wpxp, wp2, xp2, beta, sigma_sqd_w, x_tol )  &
    result( Skx )

! Description:
!   Calculate the skewness of x using the diagnostic ansatz of Larson and Golaz (2005)

! References:
!   Vincent E. Larson and Jean-Christophe Golaz, 2005: Using Probability Density
!   Functions to Derive Consistent Closure Relationships among Higher-Order Moments.
!   Mon. Wea. Rev., 133, 1023–1042.
!-------------------------------------------------------------------------------

    use constants_clubb, only: &
      one, &
      w_tol_sqd, &
      eps

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skw, &        ! Normalized Skewness of w [-]
      wpxp,&        ! Turbulent flux of x
      wp2, &        ! Variance of w            [m^2/s^2]
      xp2, &        ! Variance of x
      beta,&        ! Tunable parameter
      sigma_sqd_w,& ! Normalized variance of w [-]
      x_tol

    ! Output Variable
    real( kind = core_rknd ) :: &
      Skx     ! Result Skw [-]

    real( kind = core_rknd ) :: &
      nlzd_corr_wx, & ! Normalized correlation of w and x
      nlzd_Skw        ! Normalized skewness of w

    ! ---- Begin Code ----
    ! weberjk, 8-July 2015. Commented this out for now. cgils was failing during some tests.

      nlzd_corr_wx = ( wpxp / ( sqrt(max(wp2,w_tol_sqd)) * sqrt( max(xp2,x_tol**2) ) ) ) &
                     / sqrt( one - sigma_sqd_w )

      nlzd_Skw = Skw * ( one - sigma_sqd_w) ** (-3.0_core_rknd / 2.0_core_rknd)

      ! Larson and Golaz (2005) eq. 33
      Skx = nlzd_Skw * nlzd_corr_wx * ( beta + (one - beta) * nlzd_corr_wx**2 )

      Skx = 0._core_rknd

    return
  end function LG_2005_ansatz
!-----------------------------------------------------------------------
end module Skx_module
