subroutine stats(lat     ,pint    ,pdel    ,pstar   , &
                 div     ,t       ,q       ,nlon    )
!-----------------------------------------------------------------------
!
! Purpose:
! Accumulation of diagnostic statistics for 1 latitude.
!
! Author:  J. Rosinski
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use commap

  implicit none

#include <comsta.h>

!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: lat                ! latitude index (S->N)

  real(r8), intent(in)   :: pint (plon,plevp)  ! pressure at model interfaces
  real(r8), intent(in)   :: pdel (plon,plev)   ! pdel(k) = pint(k+1) - pint(k)
  real(r8), intent(in)   :: pstar(plon)        ! ps + psr (surface pressure)
  real(r8), intent(in)   :: div  (plon,plev)   ! divergence
  real(r8), intent(in)   :: t    (plon,plev)   ! temperature
  real(r8), intent(in)   :: q    (plon,plev)   ! moisture
  integer , intent(in)   :: nlon               ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  real(r8) prat            ! pdel(i,k)/pint(i,plevp)
  integer  i,k             ! longitude, level indices
!
!-----------------------------------------------------------------------
!
! Compute statistics for current latitude line
!
  rmsz (lat) = 0._r8
  rmsd (lat) = 0._r8
  rmst (lat) = 0._r8
  stq  (lat) = 0._r8
  psurf(lat) = 0._r8

  do i = 1,nlon
     psurf(lat) = psurf(lat) + pstar(i)
  end do
  do k = 1,plev
     do i = 1,nlon
        prat      = pdel(i,k)/pint(i,plevp)
        rmsd(lat) = rmsd(lat) + div(i,k)*div(i,k)*prat
        rmst(lat) = rmst(lat) + (t(i,k)**2)*prat
        stq (lat) = stq(lat) + q(i,k)*pdel(i,k)
     end do
  end do
  psurf(lat) = w(lat)*psurf(lat)/nlon
  rmsd (lat) = w(lat)*rmsd(lat)/nlon
  rmst (lat) = w(lat)*rmst(lat)/nlon
  stq  (lat) = w(lat)*stq(lat)/nlon
! 
  return
end subroutine stats

