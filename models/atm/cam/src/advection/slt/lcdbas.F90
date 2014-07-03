
subroutine lcdbas(grd     ,dbas2   ,dbas3   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate weights used to evaluate derivative estimates at the
! inner grid points of a four point stencil based on Lagrange
! cubic polynomial through four unequally spaced points.
! 
! Method: 
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

!------------------------------Arguments--------------------------------
  real(r8), intent(in) :: grd(4)    ! grid stencil
  real(r8), intent(out):: dbas2(4)  ! derivatives at grid point 2.
  real(r8), intent(out):: dbas3(4)  ! derivatives at grid point 3.
!
!  grd    Coordinate values of four points in stencil.
!  dbas2  Derivatives of the four basis functions at grid point 2.
!  dbas3  Derivatives of the four basis functions at grid point 3.
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  real(r8) x1                   !  |
  real(r8) x2                   !  |- grid values
  real(r8) x3                   !  |
  real(r8) x4                   !  |
  real(r8) x1mx2                !  |
  real(r8) x1mx3                !  |
  real(r8) x1mx4                !  |- differences of grid values
  real(r8) x2mx3                !  |
  real(r8) x2mx4                !  |
  real(r8) x3mx4                !  |
!-----------------------------------------------------------------------
!
  x1 = grd(1)
  x2 = grd(2)
  x3 = grd(3)
  x4 = grd(4)
  x1mx2 = x1 - x2
  x1mx3 = x1 - x3
  x1mx4 = x1 - x4
  x2mx3 = x2 - x3
  x2mx4 = x2 - x4
  x3mx4 = x3 - x4

  dbas2(1) =   x2mx3 * x2mx4 / ( x1mx2 * x1mx3 * x1mx4 )
  dbas2(2) =   -1._r8/x1mx2 + 1._r8/x2mx3 + 1._r8/x2mx4
  dbas2(3) = - x1mx2 * x2mx4 / ( x1mx3 * x2mx3 * x3mx4 )
  dbas2(4) =   x1mx2 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )

  dbas3(1) = - x2mx3 * x3mx4 / ( x1mx2 * x1mx3 * x1mx4 )
  dbas3(2) =   x1mx3 * x3mx4 / ( x1mx2 * x2mx3 * x2mx4 )
  dbas3(3) =   -1._r8/x1mx3 - 1._r8/x2mx3 + 1._r8/x3mx4
  dbas3(4) = - x1mx3 * x2mx3 / ( x1mx4 * x2mx4 * x3mx4 )

  return
end subroutine lcdbas

