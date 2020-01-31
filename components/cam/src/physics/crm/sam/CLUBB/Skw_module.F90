!$Id: Skw_module.F90 5999 2012-12-18 23:53:13Z raut@uwm.edu $
!-------------------------------------------------------------------------------
module Skw_module

  implicit none

  private ! Default Scope

  public :: Skw_func

  contains

!-------------------------------------------------------------------------------
  elemental function Skw_func( wp2, wp3 )  &
    result( Skw )

! Description:
!   Calculate the skewness of w, Skw.

! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only:  &
      w_tol_sqd,  &! Constant for w_{_tol}^2, i.e. threshold for vertical velocity
      Skw_max_mag ! Max magnitude of skewness

    use clubb_precision, only: &
      core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: min, max

    ! Parameter Constants
    ! Factor to decrease sensitivity in the denominator of Skw calculation
    real( kind = core_rknd ), parameter :: &
      Skw_denom_coef = 8.0_core_rknd ! [-]

    ! Whether to apply clipping to the final result
    logical, parameter ::  & 
      l_clipping_kluge = .false.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: & 
      wp2,  & ! w'^2    [m^2/s^2]
      wp3     ! w'^3    [m^3/s^3]

    ! Output Variable
    real( kind = core_rknd ) :: & 
      Skw     ! Result Skw [-]

    ! ---- Begin Code ----

    !Skw = wp3 / ( max( wp2, w_tol_sqd ) )**1.5_core_rknd
    ! Calculation of skewness to help reduce the sensitivity of this value to
    ! small values of wp2.
    Skw = wp3 / ( wp2 + Skw_denom_coef * w_tol_sqd )**1.5_core_rknd

    ! This is no longer needed since clipping is already
    ! imposed on wp2 and wp3 elsewhere in the code
    if ( l_clipping_kluge ) then
      Skw = min( max( Skw, -Skw_max_mag ), Skw_max_mag )
    end if

    return
  end function Skw_func
!-----------------------------------------------------------------------

end module Skw_module
