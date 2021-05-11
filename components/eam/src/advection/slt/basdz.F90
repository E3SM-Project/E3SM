
subroutine basdz(pkdim   ,sig     ,lbasdz  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute weights for the calculation of derivative estimates at two
! center points of the four point stencil for each interval in the
! unequally spaced vertical grid (as defined by the array sig).
! Estimates are from differentiating a Lagrange cubic polynomial
! through the four point stencil.
! 
! Method: 
!  pkdim   Number of grid points in vertical grid.
!  sig     Sigma values in the vertical grid.
!  lbasdz  Weights for derivative estimates based on Lagrange cubic
!          polynomial on the unequally spaced vertical grid.
!          If grid interval j is surrounded by a 4 point stencil,
!          then the derivative at the "top" of the interval (smaller
!          sigma value) uses the weights lbasdz(1,1,j),lbasdz(2,1,j),
!          lbasdz(3,1,j), and lbasdz(4,1,j).  The derivative at the
!          "bottom" of the interval uses lbasdz(1,2,j), lbasdz(2,2,j),
!          lbasdz(3,2,j), and lbasdz(4,2,j).  (Recall the vertical
!          level indices increase from the top of the atmosphere
!          towards the bottom.)
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: pkdim      ! vertical dimension
  real(r8), intent(in) :: sig(pkdim) ! sigma levels (actually a generic vert. coord)
  real(r8), intent(out):: lbasdz(4,2,pkdim) ! vertical interpolation weights
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer kk                ! index
!-----------------------------------------------------------------------
!
  do kk = 2,pkdim-2
     call lcdbas( sig(kk-1), lbasdz(1,1,kk), lbasdz(1,2,kk) )
  end do
!
  return
end subroutine basdz

