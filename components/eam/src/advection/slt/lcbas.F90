
subroutine lcbas (grd, bas1, bas2)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Evaluate the partial Lagrangian cubic basis functions (denominator
! only ) for the grid points and gather grid values
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
  real(r8), intent(in) :: grd(4)  ! grid stencil
  real(r8), intent(out):: bas1(4) ! grid values on stencil
  real(r8), intent(out):: bas2(4) ! lagrangian basis functions
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  real(r8) x0mx1     ! |
  real(r8) x0mx2     ! |
  real(r8) x0mx3     ! |- grid value differences used in weights
  real(r8) x1mx2     ! |  
  real(r8) x1mx3     ! |
  real(r8) x2mx3     ! |
!-----------------------------------------------------------------------
!
  x0mx1   = grd(1) - grd(2)
  x0mx2   = grd(1) - grd(3)
  x0mx3   = grd(1) - grd(4)
  x1mx2   = grd(2) - grd(3)
  x1mx3   = grd(2) - grd(4)
  x2mx3   = grd(3) - grd(4)

  bas1(1) = grd(1)
  bas1(2) = grd(2)
  bas1(3) = grd(3)
  bas1(4) = grd(4)

  bas2(1) =  1._r8/ ( x0mx1 * x0mx2 * x0mx3 )
  bas2(2) = -1._r8/ ( x0mx1 * x1mx2 * x1mx3 )
  bas2(3) =  1._r8/ ( x0mx2 * x1mx2 * x2mx3 )
  bas2(4) = -1._r8/ ( x0mx3 * x1mx3 * x2mx3 )

  return
end subroutine lcbas

