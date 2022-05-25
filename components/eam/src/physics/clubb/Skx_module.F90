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
  subroutine Skx_func( nz, ngrdcol, xp2, xp3, &
                       x_tol, Skw_denom_coef, Skw_max_mag, &
                       Skx )

    ! Description:
    ! Calculate the skewness of x

    ! References:
    ! None
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        core_rknd         ! Variable(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    ! External
    intrinsic :: min, max

    ! Parameter Constants
    ! Whether to apply clipping to the final result
    logical, parameter ::  &
      l_clipping_kluge = .false.

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      xp2,   & ! <x'^2>               [(x units)^2]
      xp3      ! <x'^3>               [(x units)^3]

    real( kind = core_rknd ), intent(in) :: &
      x_tol,          & ! x tolerance value                       [(x units)]
      Skw_denom_coef, & ! CLUBB tunable parameter Skw_denom_coef  [-]
      Skw_max_mag       ! Max magnitude of skewness               [-]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(out) :: &
      Skx      ! Skewness of x        [-]

    ! Local Variable
    real( kind = core_rknd ) :: &
      Skx_denom_tol
      
    integer :: i, k

    ! ---- Begin Code ----

    Skx_denom_tol = Skw_denom_coef * x_tol**2

    !Skx = xp3 / ( max( xp2, x_tol**two ) )**three_halves
    ! Calculation of skewness to help reduce the sensitivity of this value to
    ! small values of xp2.
    do k = 1, nz
      do i = 1, ngrdcol
        Skx(i,k) = xp3(i,k) / ( ( xp2(i,k) + Skx_denom_tol ) * sqrt( xp2(i,k) + Skx_denom_tol ) )
      end do
    end do

    ! This is no longer needed since clipping is already
    ! imposed on wp2 and wp3 elsewhere in the code

    ! I turned clipping on in this local copy since thlp3 and rtp3 are not clipped
    if ( l_clipping_kluge ) then
      Skx = min( max( Skx, -Skw_max_mag ), Skw_max_mag )
    end if

    return

  end subroutine Skx_func

  !-----------------------------------------------------------------------------
  subroutine LG_2005_ansatz( nz, ngrdcol, Skw, wpxp, wp2, &
                             xp2, beta, sigma_sqd_w, x_tol, &
                             Skx )

    ! Description:
    ! Calculate the skewness of x using the diagnostic ansatz of Larson and
    ! Golaz (2005).

    ! References:
    ! Vincent E. Larson and Jean-Christophe Golaz, 2005:  Using Probability
    ! Density Functions to Derive Consistent Closure Relationships among
    ! Higher-Order Moments.  Mon. Wea. Rev., 133, 1023–1042.
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one,          & ! Variable(s)
        w_tol_sqd

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none

    integer, intent(in) :: &
      nz, &
      ngrdcol

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      Skw,         & ! Skewness of w                  [-]
      wpxp,        & ! Turbulent flux of x            [m/s (x units)]
      wp2,         & ! Variance of w                  [m^2/s^2]
      xp2,         & ! Variance of x                  [(x units)^2]
      sigma_sqd_w    ! Normalized variance of w       [-]
      
    real( kind = core_rknd ), intent(in) :: &
      beta,        & ! Tunable parameter              [-]
      x_tol          ! Minimum tolerance of x         [(x units)]

    ! Output Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      Skx            ! Skewness of x                  [-]

    ! Local Variables
    real( kind = core_rknd ) :: &
      nrmlzd_corr_wx, & ! Normalized correlation of w and x       [-]
      nrmlzd_Skw        ! Normalized skewness of w                [-]
      
    integer :: i, k

    ! ---- Begin Code ----
    ! weberjk, 8-July 2015. Commented this out for now. cgils was failing during some tests.

    ! Larson and Golaz (2005) eq. 16
    do k = 1, nz
      do i = 1, ngrdcol
        nrmlzd_corr_wx = &
                wpxp(i,k) / sqrt( max( wp2(i,k), w_tol_sqd ) &
                             * max( xp2(i,k), x_tol**2 ) * ( one - sigma_sqd_w(i,k) ) )

        ! Larson and Golaz (2005) eq. 11
        nrmlzd_Skw = Skw(i,k) / ( ( one - sigma_sqd_w(i,k)) * sqrt( one - sigma_sqd_w(i,k) ) )

        ! Larson and Golaz (2005) eq. 33
        Skx(i,k) = nrmlzd_Skw * nrmlzd_corr_wx &
              * ( beta + ( one - beta ) * nrmlzd_corr_wx**2 )
      end do
    end do

    return

  end subroutine LG_2005_ansatz

  !-----------------------------------------------------------------------------
  subroutine xp3_LG_2005_ansatz( nz, ngrdcol, Skw_zt, wpxp_zt, wp2_zt, &
                                 xp2_zt, sigma_sqd_w_zt, &
                                 beta, Skw_denom_coef, x_tol, &
                                 xp3 )
    ! Description:
    ! Calculate <x'^3> after calculating the skewness of x using the ansatz of
    ! Larson and Golaz (2005).

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        grid ! Type

    use clubb_precision, only: &
        core_rknd ! Variable(s)

    implicit none
    
    integer, intent(in) :: &
      nz, &
      ngrdcol

    ! Input Variables
    real( kind = core_rknd ), dimension(ngrdcol,nz), intent(in) :: &
      Skw_zt,         & ! Skewness of w on thermodynamic levels   [-]
      wpxp_zt,        & ! Flux of x  (interp. to t-levs.)         [m/s(x units)]
      wp2_zt,         & ! Variance of w (interp. to t-levs.)      [m^2/s^2]
      xp2_zt,         & ! Variance of x (interp. to t-levs.)      [(x units)^2]
      sigma_sqd_w_zt    ! Normalized variance of w (interp. to t-levs.)   [-]

    real( kind = core_rknd ), intent(in) :: &
      beta,           & ! CLUBB tunable parameter beta            [-]
      Skw_denom_coef, & ! CLUBB tunable parameter Skw_denom_coef  [-]
      x_tol             ! Minimum tolerance of x                  [(x units)]

    ! Return Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      xp3    ! <x'^3> (thermodynamic levels)    [(x units)^3]

    ! Local Variable
    real( kind = core_rknd ), dimension(ngrdcol,nz) :: &
      Skx_zt, &    ! Skewness of x on thermodynamic levels    [-]
      Skx_denom_tol

    ! ---- Begin Code ----

    Skx_denom_tol = Skw_denom_coef * x_tol**2

    ! Calculate skewness of x using the ansatz of LG05.
    call LG_2005_ansatz( nz, ngrdcol, Skw_zt, wpxp_zt, wp2_zt, &
                         xp2_zt, beta, sigma_sqd_w_zt, x_tol, &
                         Skx_zt )

    ! Calculate <x'^3> using the reverse of the special sensitivity reduction
    ! formula in function Skx_func above.
    xp3 = Skx_zt * ( xp2_zt + Skx_denom_tol ) * sqrt( xp2_zt + Skx_denom_tol )


    return

  end subroutine xp3_LG_2005_ansatz

  !-----------------------------------------------------------------------------

end module Skx_module
