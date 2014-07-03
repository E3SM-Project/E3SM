subroutine etadtn(lat   ,nlon    ,d       ,edstar  ,etadot  )
!
!-----------------------------------------------------------------------
!
! Purpose:
! to compute (1/ps)etadot(dp/deta) for time steps subsequent to nstep = 0.
! This routine completes a partial computation which began in the previous
! time step (SCANDYN)
!
! NOTE:  Computing only the bottom "pver" levels of a half-level field.
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use hycoef, only : hypd, hybi, hypi, nprlev
  implicit none

!------------------------------Arguments--------------------------------
!
  integer, intent(in) :: lat                     ! latitude
  integer, intent(in) :: nlon                    ! number of longitudes

  real(r8), intent(in)   :: d     (plon,plev)    ! divergence
  real(r8), intent(in)   :: edstar(plon,plev)    ! partially computed etadot
  real(r8), intent(out)  :: etadot(plon,plev)    ! vertical velocity in eta coordinates
!
!-----------------------------------------------------------------------
!
! local workspace
!
  integer i,k               ! indices
  real(r8) ddpn(plon)       ! divergence integral
  real(r8) ddpk(plon)       ! divergence accumulator
  real(r8) dpr (plev)       ! pressure coefficient
!
!-----------------------------------------------------------------------
!
  do i = 1,nlon
     ddpn(i) = 0._r8
     ddpk(i) = 0._r8
  end do
!
  do k = 1,plev
     dpr(k) = hypd(k)/hypi(plevp)
  end do
!
  do k = 1,plev
     do i = 1,nlon
        ddpn(i) = ddpn(i) + d(i,k)*dpr(k)
     end do
  end do
!
  do k = 1,plev-1
     do i = 1,nlon
        ddpk(i) = ddpk(i) + d(i,k)*dpr(k)
        etadot(i,k) = edstar(i,k) - ddpk(i)
     end do
     if(k.ge.nprlev) then
        do i = 1,nlon
           etadot(i,k) = etadot(i,k) + hybi(k+1)*ddpn(i)
        end do
     endif
  end do
!
! Finish bottom level
!
  do i = 1,nlon
     etadot(i,plev) = 0._r8
  end do
!
  return
end subroutine etadtn

