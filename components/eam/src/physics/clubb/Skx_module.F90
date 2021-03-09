!-------------------------------------------------------------------------
!$Id$
!===============================================================================
module Skx_module

  implicit none

  private ! Default Scope

  public :: Skx_func, &
            LG_2005_ansatz, &
            xp3_LG_2005_ansatz

  contains

  !-----------------------------------------------------------------------------
  function Skx_func( xp2, xp3, x_tol )  &
  result( Skx )

    ! Description:
    ! Calculate the skewness of x

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three_halves      ! 3/2

    use parameters_tunable, only: &
        Skw_denom_coef, & ! Variable(s)
        Skw_max_mag       ! Max magnitude of skewness

    use clubb_precision, only: &
        core_rknd         ! Variable(s)

    use grid_class, only: &
        gr                ! Variable Type

    implicit none

    ! External
    intrinsic :: min, max

    ! Parameter Constants
    ! Whether to apply clipping to the final result
    logical, parameter ::  &
      l_clipping_kluge = .false.

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      x_tol    ! x tolerance value    [(x units)]

    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      xp2,   & ! <x'^2>               [(x units)^2]
      xp3      ! <x'^3>               [(x units)^3]

    ! Output Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      Skx      ! Skewness of x        [-]

    ! Local Variable
    real( kind = core_rknd ) :: &
      Skx_denom_tol

    ! ---- Begin Code ----

    Skx_denom_tol = Skw_denom_coef * x_tol**2

    !Skx = xp3 / ( max( xp2, x_tol**two ) )**three_halves
    ! Calculation of skewness to help reduce the sensitivity of this value to
    ! small values of xp2.
    Skx = xp3 / ( ( xp2 + Skx_denom_tol ) * sqrt( xp2 + Skx_denom_tol ) )

    ! This is no longer needed since clipping is already
    ! imposed on wp2 and wp3 elsewhere in the code

    ! I turned clipping on in this local copy since thlp3 and rtp3 are not clipped
    if ( l_clipping_kluge ) then
      Skx = min( max( Skx, -Skw_max_mag ), Skw_max_mag )
    end if

    return

  end function Skx_func

  !-----------------------------------------------------------------------------
  elemental function LG_2005_ansatz( Skw, wpxp, wp2, &
                                     xp2, beta, sigma_sqd_w, x_tol ) &
  result( Skx )

    ! Description:
    ! Calculate the skewness of x using the diagnostic ansatz of Larson and
    ! Golaz (2005).

    ! References:
    ! Vincent E. Larson and Jean-Christophe Golaz, 2005:  Using Probability
    ! Density Functions to Derive Consistent Closure Relationships among
    ! Higher-Order Moments.  Mon. Wea. Rev., 133, 1023–1042.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        three_halves, & ! Variable(s)
        one,          &
        w_tol_sqd

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt

    ! Input Variables
    real( kind = core_rknd ), intent(in) :: &
      Skw,         & ! Skewness of w                  [-]
      wpxp,        & ! Turbulent flux of x            [m/s (x units)]
      wp2,         & ! Variance of w                  [m^2/s^2]
      xp2,         & ! Variance of x                  [(x units)^2]
      beta,        & ! Tunable parameter              [-]
      sigma_sqd_w, & ! Normalized variance of w       [-]
      x_tol          ! Minimum tolerance of x         [(x units)]

    ! Output Variable
    real( kind = core_rknd ) :: &
      Skx            ! Skewness of x                  [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      nrmlzd_corr_wx, & ! Normalized correlation of w and x       [-]
      nrmlzd_Skw        ! Normalized skewness of w                [-]

    ! ---- Begin Code ----
    ! weberjk, 8-July 2015. Commented this out for now. cgils was failing during some tests.

    ! Larson and Golaz (2005) eq. 16
    nrmlzd_corr_wx &
    = wpxp / sqrt( max( wp2, w_tol_sqd ) * max( xp2, x_tol**2 ) * ( one - sigma_sqd_w ) )

    ! Larson and Golaz (2005) eq. 11
    nrmlzd_Skw = Skw / ( ( one - sigma_sqd_w) * sqrt( one - sigma_sqd_w ) )

    ! Larson and Golaz (2005) eq. 33
    Skx = nrmlzd_Skw * nrmlzd_corr_wx &
          * ( beta + ( one - beta ) * nrmlzd_corr_wx**2 )


    return

  end function LG_2005_ansatz

  !-----------------------------------------------------------------------------
  function xp3_LG_2005_ansatz( Skw_zt, wpxp_zt, wp2_zt, &
                               xp2_zt, sigma_sqd_w_zt, x_tol ) &
  result( xp3 )

    ! Description:
    ! Calculate <x'^3> after calculating the skewness of x using the ansatz of
    ! Larson and Golaz (2005).

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        gr    ! Variable Type

    use constants_clubb, only: &
        three_halves    ! Variable(s)

    use parameters_tunable, only: &
        beta,           & ! Variable(s)
        Skw_denom_coef

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    ! External
    intrinsic :: sqrt

    ! Input Variables
    real( kind = core_rknd ), dimension(gr%nz), intent(in) :: &
      Skw_zt,         & ! Skewness of w on thermodynamic levels   [-]
      wpxp_zt,        & ! Flux of x  (interp. to t-levs.)         [m/s(x units)]
      wp2_zt,         & ! Variance of w (interp. to t-levs.)      [m^2/s^2]
      xp2_zt,         & ! Variance of x (interp. to t-levs.)      [(x units)^2]
      sigma_sqd_w_zt    ! Normalized variance of w (interp. to t-levs.)   [-]

    real( kind = core_rknd ), intent(in) :: &
      x_tol    ! Minimum tolerance of x         [(x units)]

    ! Return Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      xp3    ! <x'^3> (thermodynamic levels)    [(x units)^3]

    ! Local Variable
    real( kind = core_rknd ), dimension(gr%nz) :: &
      Skx_zt, &    ! Skewness of x on thermodynamic levels    [-]
      Skx_denom_tol

    ! ---- Begin Code ----

    Skx_denom_tol = Skw_denom_coef * x_tol**2

    ! Calculate skewness of x using the ansatz of LG05.
    Skx_zt(1:gr%nz) &
    = LG_2005_ansatz( Skw_zt(1:gr%nz), wpxp_zt(1:gr%nz), wp2_zt(1:gr%nz), &
                      xp2_zt(1:gr%nz), beta, sigma_sqd_w_zt(1:gr%nz), x_tol )

    ! Calculate <x'^3> using the reverse of the special sensitivity reduction
    ! formula in function Skx_func above.
    xp3 = Skx_zt * ( xp2_zt + Skx_denom_tol ) * sqrt( xp2_zt + Skx_denom_tol )


    return

  end function xp3_LG_2005_ansatz

  !-----------------------------------------------------------------------------

end module Skx_module
