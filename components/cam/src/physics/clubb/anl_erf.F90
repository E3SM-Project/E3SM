!-----------------------------------------------------------------------
! $Id: anl_erf.F90 7269 2014-09-04 21:00:07Z raut@uwm.edu $
!===============================================================================
module anl_erf

  implicit none

  public :: dp_erf, &
            dp_erfc, &
            erf, &
            erfc

  private :: cr_erf, &
             cr_erfc

  ! The interfaces allow us to avoid a compiler warning about
  ! shadowing the intrinsic functions
  interface erf
    module procedure cr_erf
  end interface

  interface erfc
    module procedure cr_erfc
  end interface

  private ! Default Scope

  contains

  !=============================================================================
  pure function cr_erf( x ) result( erfx_core_rknd )
    ! Description:
    !   Calls dp_erf after casting x to double precision.
    !   This allows CLUBB to run erf even when core_rknd is in single precision.
    !
    !  Arguments:
    !    Input, real ( kind = dp ) x, the argument of ERF.
    !    Output, real ( kind = core_rknd ) erfx_core_rknd, the value of ERF(X).
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        dp, & ! Constants
        core_rknd

    implicit none

    ! Input Variables(s)
    real( kind = core_rknd), intent(in) :: x

    ! Return type
    real( kind = core_rknd ) :: erfx_core_rknd

    ! Local Variables
    real( kind = dp) :: x_dp, erfx_dp

    ! Cast the input to dp
    x_dp = real( x, kind = dp )

    ! Call the function with the correct argument
    erfx_dp = dp_erf( x_dp )

    ! Get the output in core_rknd
    erfx_core_rknd = real( erfx_dp, kind = core_rknd )

    return

  end function cr_erf

  !=============================================================================
  pure function dp_erf( x ) result( erfx )

    ! Description:
    !   DP_ERF evaluates the error function DP_ERF(X).
    !
    !   Original Author:
    !     William Cody,
    !     Mathematics and Computer Science Division,
    !     Argonne National Laboratory,
    !     Argonne, Illinois, 60439.
    !
    !   References:
    !     William Cody,
    !     "Rational Chebyshev approximations for the error function",
    !     Mathematics of Computation,
    !     1969, pages 631-638.
    !
    !  Arguments:
    !    Input, real ( kind = dp ) X, the argument of ERF.
    !    Output, real ( kind = dp ) ERFX, the value of ERF(X).
    !
    ! Modifications:
    !   kind = 8 was replaced by the more portable sp and dp by UWM.
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        dp, & ! Constants
        core_rknd

    implicit none

    ! Input Variables(s)
    real( kind = dp ), intent(in) :: x

    ! External
    intrinsic :: epsilon, exp, aint

    ! Local Constants
    real( kind = dp ), parameter, dimension( 5 ) ::  & 
    a = (/ 3.16112374387056560E+00_dp, &
           1.13864154151050156E+02_dp, &
           3.77485237685302021E+02_dp, &
           3.20937758913846947E+03_dp, &
           1.85777706184603153E-01_dp /)
    real( kind = dp ), parameter, dimension( 4 ) ::  & 
    b = (/ 2.36012909523441209E+01_dp, &
           2.44024637934444173E+02_dp, &
           1.28261652607737228E+03_dp, &
           2.84423683343917062E+03_dp /)
    real( kind = dp ), parameter, dimension( 9 ) ::  & 
    c = (/ 5.64188496988670089E-01_dp, &
           8.88314979438837594E+00_dp, &
           6.61191906371416295E+01_dp, &
           2.98635138197400131E+02_dp, &
           8.81952221241769090E+02_dp, &
           1.71204761263407058E+03_dp, &
           2.05107837782607147E+03_dp, &
           1.23033935479799725E+03_dp, &
           2.15311535474403846E-08_dp /)
    real( kind = dp ), parameter, dimension( 8 ) :: & 
    d = (/ 1.57449261107098347E+01_dp, &
           1.17693950891312499E+02_dp, &
           5.37181101862009858E+02_dp,  &
           1.62138957456669019E+03_dp, &
           3.29079923573345963E+03_dp, &
           4.36261909014324716E+03_dp, &
           3.43936767414372164E+03_dp, &
           1.23033935480374942E+03_dp /)
    real( kind = dp ), parameter, dimension( 6 ) ::  & 
    p = (/ 3.05326634961232344E-01_dp, &
           3.60344899949804439E-01_dp, &
           1.25781726111229246E-01_dp, &
           1.60837851487422766E-02_dp, &
           6.58749161529837803E-04_dp, &
           1.63153871373020978E-02_dp /)

    real( kind = dp ), parameter, dimension( 5 ) ::  & 
    q = (/ 2.56852019228982242E+00_dp, &
           1.87295284992346047E+00_dp, &
           5.27905102951428412E-01_dp, &
           6.05183413124413191E-02_dp, &
           2.33520497626869185E-03_dp /)

    real( kind = dp ), parameter ::  & 
    SQRPI  = 0.56418958354775628695E+00_dp, &
    THRESH = 0.46875E+00_dp, &
    XBIG   = 26.543E+00_dp

    ! Return type
    real( kind = dp ) :: erfx

    ! Local variables
    real( kind = dp ) ::  & 
    del, & 
    xabs, &
    xden, & 
    xnum, & 
    xsq

    integer :: i ! Index

    !-----------------------------------------------------------------------
    ! Get the abs value of xabs - schemena 20140827
    xabs = abs( x )

    !
    !  Evaluate ERF(X) for |X| <= 0.46875.
    !
    if ( xabs <= THRESH ) then

      if ( epsilon( xabs ) < xabs ) then
        xsq = xabs * xabs
      else
        xsq = 0.0E+00_dp
      end if

      xnum = a(5) * xsq
      xden = xsq
      do i = 1, 3
        xnum = ( xnum + a(i) ) * xsq
        xden = ( xden + b(i) ) * xsq
      end do

      erfx = x * ( xnum + a(4) ) / ( xden + b(4) )
      !
      !  Evaluate ERFC(X) for 0.46875 <= |X| <= 4.0.
      !
    else if ( xabs <= 4.0E+00_dp ) then

      xnum = c(9) * xabs
      xden = xabs
      do i = 1, 7
        xnum = ( xnum + c(i) ) * xabs
        xden = ( xden + d(i) ) * xabs
      end do

      erfx = ( xnum + c(8) ) / ( xden + d(8) )
      xsq = aint( xabs * 16.0E+00_dp ) / 16.0E+00_dp
      del = ( xabs - xsq ) * ( xabs + xsq )
      ! xsq * xsq in the exponential was changed to xsq**2.
      ! This seems to decrease runtime by about a half a percent.
      ! ~~EIHoppe//20090622
      erfx = exp( - xsq**2 ) * exp( - del ) * erfx

      erfx = ( 0.5E+00_dp - erfx ) + 0.5E+00_dp

      if ( x < 0.0E+00_dp ) then
        erfx = - erfx
      end if
      !
      !  Evaluate ERFC(X) for 4.0 < |X|.
      !
    else

      if ( XBIG <= xabs ) then

        if ( 0.0E+00_dp < real(x, kind=dp) ) then
          erfx = 1.0E+00_dp
        else
          erfx = -1.0E+00_dp
        end if

      else

        xsq = 1.0E+00_dp / ( xabs * xabs )

        xnum = p(6) * xsq
        xden = xsq
        do i = 1, 4
          xnum = ( xnum + p(i) ) * xsq
          xden = ( xden + q(i) ) * xsq
        end do

        erfx = xsq * ( xnum + p(5) ) / ( xden + q(5) )
        erfx = ( SQRPI -  erfx ) / xabs
        xsq = aint( xabs * 16.0E+00_dp ) / 16.0E+00_dp
        del = ( xabs - xsq ) * ( xabs + xsq )
        erfx = exp( - xsq * xsq ) * exp( - del ) * erfx

        erfx = ( 0.5E+00_dp - erfx ) + 0.5E+00_dp
        if ( x < 0.0E+00_dp ) then
          erfx = - erfx
        end if

      end if

    end if

    return

  end function dp_erf

  !=============================================================================
  pure function cr_erfc( x ) result( erfcx_core_rknd )
    ! Description:
    !   Calls dp_erfc after casting x to double precision.
    !   This allows CLUBB to run erfc even when core_rknd is in single precision.
    !
    !  Arguments:
    !    Input, real ( kind = core_rknd ) x, the argument of ERFC.
    !    Output, real ( kind = core_rknd ) erfcx_core_rknd, the value of ERFC(X).
    !-----------------------------------------------------------------------

    use clubb_precision, only: &
        dp, & ! Constants
        core_rknd

    implicit none

    ! Input Variables(s)
    real( kind = core_rknd), intent(in) :: x

    ! Return type
    real( kind = core_rknd ) :: erfcx_core_rknd

    ! Local Variables
    real( kind = dp) :: x_dp, erfcx_dp

    ! Cast the input to dp
    x_dp = real( x, kind = dp )

    ! Call the function with the correct argument
    erfcx_dp = dp_erfc( x_dp )

    ! Get the output in core_rknd
    erfcx_core_rknd = real( erfcx_dp, kind = core_rknd )

    return

  end function cr_erfc

  !=============================================================================
  pure function dp_erfc( x ) result( erfcx )

    ! Description:
    ! The complimentary error function of x:
    !
    ! erfc(x) = 1 - erf(x);
    !
    ! where:
    !
    ! erf(x) = ( 2 / sqrt(pi) ) INT(0:x) e^-t^2 dt;
    !
    ! and
    !
    ! erfc(x) = ( 2 / sqrt(pi) ) INT(x:inf) e^-t^2 dt.

    ! References:
    !-----------------------------------------------------------------------

    use constants_clubb, only: &
        one_dp  ! Constant(s)

    use clubb_precision, only: &
        dp    ! Variable(s)

    implicit none

    ! Input Variable
    real( kind = dp ), intent(in) :: x

    ! Return Variable
    real( kind = dp ) :: erfcx


    erfcx = one_dp - dp_erf( x )


    return

  end function dp_erfc

!===============================================================================

end module anl_erf
