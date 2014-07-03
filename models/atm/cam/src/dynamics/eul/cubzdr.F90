
subroutine cubzdr(nlon    ,pkdim   ,f       ,lbasdz  ,dfz1    ,  &
                  dfz2    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Vertical derivative estimates for a vertical slice using Lagrangian
! cubic formulas.  
!
! Method: 
! Derivatives are set to zero at the top and bottom.
! At the "inner nodes" of the top and bottom intervals, a "one sided"
! estimate is used.
!
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, March 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon
!-----------------------------------------------------------------------
   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments 
!
   integer, intent(in) :: nlon              ! number of longitudes
   integer, intent(in) :: pkdim             ! vertical     dimension
!
   real(r8), intent(in) :: f(plon,pkdim)       ! constituent field
   real(r8), intent(in) :: lbasdz(4,2,pkdim)   ! vertical interpolation weights
!
! Output arguments
!
   real(r8), intent(out) :: dfz1(plon,pkdim)    ! derivative at top of interval
   real(r8), intent(out) :: dfz2(plon,pkdim)    ! derivative at bot of interval
!-----------------------------------------------------------------------
!
!  nlon    Number of longitudes
!  pkdim   Vertical dimension of arrays.
!  f       Vertical slice of data for which derivative estimates are
!          made
!  lbasdz  Lagrangian cubic basis functions for evaluating the
!          derivatives on the unequally spaced vertical grid.
!  dfz1    dfz1 contains derivative estimates at the "top" edges of the
!          intervals in the f array.
!  dfz2    dfz2 contains derivative estimates at the "bottom" edges of
!          the intervals in the f array.
!
!---------------------------Local variables-----------------------------
!
   integer i,k               ! indices
!
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=2,pkdim-2
      do i=1,nlon
!
! Lagrangian derivative estimates (cubic) for the two center nodes in a
! four node stencil.
!
         dfz1(i,k) = lbasdz(1,1,k)*f(i,k-1) +  &
            lbasdz(2,1,k)*f(i,k)   +  &
            lbasdz(3,1,k)*f(i,k+1) +  &
            lbasdz(4,1,k)*f(i,k+2)
!
         dfz2(i,k) = lbasdz(1,2,k)*f(i,k-1) +  &
            lbasdz(2,2,k)*f(i,k)   +  &
            lbasdz(3,2,k)*f(i,k+1) +  &
            lbasdz(4,2,k)*f(i,k+2)
      end do
   end do
!
! Constrain derivatives to zero at top and bottom of vertical grid.
! At the interior nodes of the intervals at the top and bottom of the
! vertical grid, use the derivative estimate at that same node for the
! adjacent interval.  (This is a "one-sided" estimate for that node.)
!
   do i=1,nlon
      dfz1(i,1) = 0.0_r8
      dfz2(i,1) = dfz1(i,2)
      dfz1(i,pkdim-1) = dfz2(i,pkdim-2)
      dfz2(i,pkdim-1) = 0.0_r8
   end do
!
   return
end subroutine cubzdr

