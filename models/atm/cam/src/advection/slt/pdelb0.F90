
subroutine pdelb0(ps      ,pdelb   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the pressure intervals between the interfaces for the "B"
! (surface pressure dependent) portion of the hybrid grid only. 
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
  use pmgrid,       only: plev, plon
  use hycoef,       only: hybd
  implicit none
  
!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon               ! longitude dimension
  real(r8), intent(in) :: ps(plon)           ! surface Pressure
  real(r8), intent(out):: pdelb(plon,plev)   ! pressure difference between interfaces
                                             ! (pressure defined using the "B" part  
                                             ! of the hybrid grid only)
!-----------------------------------------------------------------------

!---------------------------Local workspace-----------------------------
  integer i,k               ! longitude, level indices
!-----------------------------------------------------------------------
!
! Compute del P(B)
!
!$OMP PARALLEL DO PRIVATE (K, I)
  do k = 1,plev
     do i = 1,nlon
        pdelb(i,k) = hybd(k)*ps(i)
     end do
  end do

  return
end subroutine pdelb0

