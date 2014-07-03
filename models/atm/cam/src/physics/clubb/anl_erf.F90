! $Id: anl_erf.F90 5324 2011-07-27 21:05:45Z dschanen@uwm.edu $
module anl_erf

  implicit none

  public :: erf

  interface erf
    module procedure dp_erf, sp_erf
  end interface

  private :: dp_erf, sp_erf

  private ! Default Scope

  contains

  function dp_erf( x ) result( erfx )
!-----------------------------------------------------------------------
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
!    Input, real ( kind = 8 ) X, the argument of ERF.
!    Output, real ( kind = 8 ) ERFX, the value of ERF(X).
!-----------------------------------------------------------------------

    implicit none

    ! Input Variables(s)
    double precision, intent(in) :: x

    ! External
    intrinsic :: epsilon, exp, aint

    ! Local Constants
    real( kind = 8 ), parameter, dimension( 5 ) ::  & 
    a = (/ 3.16112374387056560D+00, & 
           1.13864154151050156D+02, & 
           3.77485237685302021D+02, & 
           3.20937758913846947D+03, & 
           1.85777706184603153D-01 /)
    real( kind = 8 ), parameter, dimension( 4 ) ::  & 
    b = (/ 2.36012909523441209D+01, & 
           2.44024637934444173D+02, & 
           1.28261652607737228D+03, & 
           2.84423683343917062D+03 /)
    real( kind = 8 ), parameter, dimension( 9 ) ::  & 
    c = (/ 5.64188496988670089D-01, & 
           8.88314979438837594D+00, & 
           6.61191906371416295D+01, & 
           2.98635138197400131D+02, & 
           8.81952221241769090D+02, & 
           1.71204761263407058D+03, & 
           2.05107837782607147D+03, & 
           1.23033935479799725D+03, & 
           2.15311535474403846D-08 /)
    real( kind = 8 ), parameter, dimension( 8 ) :: & 
    d = (/ 1.57449261107098347D+01, & 
           1.17693950891312499D+02, & 
           5.37181101862009858D+02,  & 
           1.62138957456669019D+03, & 
           3.29079923573345963D+03, & 
           4.36261909014324716D+03, & 
           3.43936767414372164D+03, & 
           1.23033935480374942D+03 /)
    real( kind = 8 ), parameter, dimension( 6 ) ::  & 
    p = (/ 3.05326634961232344D-01, & 
           3.60344899949804439D-01, & 
           1.25781726111229246D-01, & 
           1.60837851487422766D-02, & 
           6.58749161529837803D-04, & 
           1.63153871373020978D-02 /)

    real( kind = 8 ), parameter, dimension( 5 ) ::  & 
    q = (/ 2.56852019228982242D+00, & 
           1.87295284992346047D+00, & 
           5.27905102951428412D-01, & 
           6.05183413124413191D-02, & 
           2.33520497626869185D-03 /)

    real( kind = 8 ), parameter ::  & 
    SQRPI  = 0.56418958354775628695D+00, & 
    THRESH = 0.46875D+00, & 
    XBIG   = 26.543D+00

    ! Return type
    real( kind = 8 ) :: erfx

    ! Local variables
    real( kind = 8 ) ::  & 
    del, & 
    xabs, & 
    xden, & 
    xnum, & 
    xsq

    integer :: i ! Index

!-------------------------------------------------------------------------------
    xabs = abs( x )

    !
    !  Evaluate ERF(X) for |X| <= 0.46875.
    !
    if ( xabs <= THRESH ) then

      if ( epsilon( xabs ) < xabs ) then
        xsq = xabs * xabs
      else
        xsq = 0.0D+00
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
    else if ( xabs <= 4.0D+00 ) then

      xnum = c(9) * xabs
      xden = xabs
      do i = 1, 7
        xnum = ( xnum + c(i) ) * xabs
        xden = ( xden + d(i) ) * xabs
      end do

      erfx = ( xnum + c(8) ) / ( xden + d(8) )
      xsq = aint( xabs * 16.0D+00 ) / 16.0D+00
      del = ( xabs - xsq ) * ( xabs + xsq )
      ! xsq * xsq in the exponential was changed to xsq**2.
      ! This seems to decrease runtime by about a half a percent.
      ! ~~EIHoppe//20090622
      erfx = exp( - xsq**2 ) * exp( - del ) * erfx

      erfx = ( 0.5D+00 - erfx ) + 0.5D+00

      if ( x < 0.0D+00 ) then
        erfx = - erfx
      end if
      !
      !  Evaluate ERFC(X) for 4.0 < |X|.
      !
    else

      if ( XBIG <= xabs ) then

        if ( 0.0D+00 < x ) then
          erfx = 1.0D+00
        else
          erfx = -1.0D+00
        end if

      else

        xsq = 1.0D+00 / ( xabs * xabs )

        xnum = p(6) * xsq
        xden = xsq
        do i = 1, 4
          xnum = ( xnum + p(i) ) * xsq
          xden = ( xden + q(i) ) * xsq
        end do

        erfx = xsq * ( xnum + p(5) ) / ( xden + q(5) )
        erfx = ( SQRPI -  erfx ) / xabs
        xsq = aint( xabs * 16.0D+00 ) / 16.0D+00
        del = ( xabs - xsq ) * ( xabs + xsq )
        erfx = exp( - xsq * xsq ) * exp( - del ) * erfx

        erfx = ( 0.5D+00 - erfx ) + 0.5D+00
        if ( x < 0.0D+00 ) then
          erfx = - erfx
        end if

      end if

    end if

    return
  end function dp_erf

!-----------------------------------------------------------------------
  function sp_erf( x ) result( erfx )

! Description:
!   Return a truncation of the 64bit approx. of the error function.
!   Ideally we would probably use a 32bit table for our approx.

! References:
!   None
!-----------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: real

    ! Input Variables
    real( kind=4 ), intent(in) :: x

    ! Return type
    real( kind=4 ) :: erfx

    erfx = real( dp_erf( real(x, kind=8) ), kind=4 )

    return
  end function sp_erf

end module anl_erf
