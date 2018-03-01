
subroutine basdy(phi     ,lbasdy  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute weights for the calculation of derivative estimates at the two
! center points of the four point stencil for each interval in the
! unequally spaced latitude grid. Estimates are from differentiating
! a Lagrange cubic polynomial through the four point stencil.
! 
! Method: 
!  phi     Latitude values in the extended grid.
!  lbasdy  Weights for derivative estimates based on Lagrange cubic
!          polynomial on the unequally spaced latitude grid.
!          If grid interval j (in extended grid) is surrounded by
!          a 4 point stencil, then the derivative at the "bottom"
!          of the interval uses the weights lbasdy(1,1,j),
!          lbasdy(2,1,j), lbasdy(3,1,j), and lbasdy(4,1,j).
!          The derivative at the "top" of the interval
!          uses lbasdy(1,2,j), lbasdy(2,2,j), lbasdy(3,2,j),
!          and lbasdy(4,2,j).
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use scanslt,      only: nxpt, platd
  implicit none

!------------------------------Parameters-------------------------------
  integer, parameter ::  jfirst = nxpt + 1          ! first index to be computed
  integer, parameter ::  jlast  = platd - nxpt - 1  ! last  index to be computed
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
  real(r8), intent(in)  :: phi(platd)          ! latitude coordinates of model grid
  real(r8), intent(out) :: lbasdy(4,2,platd)   ! derivative estimate weights
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer jj                ! index
!-----------------------------------------------------------------------
!
  do jj = jfirst,jlast
     call lcdbas( phi(jj-1), lbasdy(1,1,jj), lbasdy(1,2,jj) )
  end do
!
  return
end subroutine basdy

