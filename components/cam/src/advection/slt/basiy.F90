
subroutine basiy(phi     ,lbasiy  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute weights used in Lagrange cubic polynomial interpolation in 
! the central interval of a four point stencil.  Done for each interval
! in the unequally spaced latitude grid.
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
  use shr_kind_mod, only: r8 => shr_kind_r8
  use scanslt,      only: nxpt, platd
  implicit none

!------------------------------Parameters-------------------------------
  integer, parameter ::  jfirst = nxpt + 1          ! first index to be computed
  integer, parameter ::  jlast  = platd - nxpt - 1  ! last  index to be computed
!-----------------------------------------------------------------------

!------------------------------Arguments--------------------------------
  real(r8), intent(in)  :: phi(platd)           ! grid values in extended grid
  real(r8), intent(out) :: lbasiy(4,2,platd)    ! Weights for Lagrange cubic interp
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer jj                ! index
!-----------------------------------------------------------------------
!
  do jj = jfirst,jlast
     call lcbas( phi(jj-1),lbasiy(1,1,jj),lbasiy(1,2,jj) )
  end do
!
  return
end subroutine basiy

