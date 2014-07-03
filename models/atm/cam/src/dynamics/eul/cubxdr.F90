subroutine cubxdr(pidim   ,ibeg    ,len     ,dx      ,f       ,  &
                  fxl     ,fxr     )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute Lagrangian cubic derivative estimates for data on an equally
! spaced grid.
!
! Method: 
! Compute Lagrangian cubic derivative estimates for data on an equally
! spaced grid.  Suppose grid interval i is centered in a 4 point
! stencil consisting of grid points i-1, i, i+1, and i+2.  Then the
! derivative at the left edge of the interval (i.e., grid point i)
! is stored in fxl(i), and the derivative at the right edge of the
! interval (i.e., grid point i+1) is stored in fxr(i).  Note that
! fxl(i) is not necessarily equal to fxr(i-1) even though both of
! these values are estimates of the derivative at grid point i.
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
  integer, intent(in) :: pidim             ! dimension
  integer, intent(in) :: ibeg              ! starting index to perform computation
  integer, intent(in) :: len               ! length over which to perform comp.
!
  real(r8), intent(in) :: dx               ! grid interval
  real(r8), intent(in) :: f(pidim)         ! input field values
!
! Output arguments
!
  real(r8), intent(out) :: fxl(pidim)           ! left  derivative of interval i in "f"
  real(r8), intent(out) :: fxr(pidim)           ! right derivative of interval i in "f"
!-----------------------------------------------------------------------
!
!  pidim   Length of f, fxl, and fxr.
!  ibeg    First interval of grid for which derivatives are computed.
!  len     Number of grid intervals for which derivatives are computed.
!          (There are pidim - 1 intervals between the pidim gridpoints
!          represented in f, fxl, and fxr.)
!  dx      Value of grid spacing.
!  f       Values on equally spaced grid for which derivatives are
!          computed.
!  fxl     fxl(i) is the derivative at the left  edge of interval i.
!  fxr     fxr(i) is the derivative at the right edge of interval i.
!
!---------------------------Local variables-----------------------------
!
  integer i                 ! index
  integer iend              ! index denoting end of computation
!
  real(r8) rdx6                 ! normalization weight
!
!-----------------------------------------------------------------------
!
  iend = ibeg + len - 1
  rdx6 = 1._r8/(6._r8*dx)
!
  do i = ibeg,iend
     fxl(i) = ( -2._r8*f(i-1) - 3._r8*f(i) + 6._r8*f(i+1) -    f(i+2) )*rdx6
     fxr(i) = (     f(i-1) - 6._r8*f(i) + 3._r8*f(i+1) + 2._r8*f(i+2) )*rdx6
  end do
!
  return
end subroutine cubxdr

