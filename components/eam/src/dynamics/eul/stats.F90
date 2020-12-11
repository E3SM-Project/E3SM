subroutine stats(lat     ,pint    ,pdel    ,pstar   , &
                 vort    ,div     ,t       ,q       ,nlon    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Accumulation of diagnostic statistics for 1 latitude.
! 
! Method: 
! 
! Author: 
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, J. Hack, August 1992
! Reviewed:          D. Williamson, March 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, plevp, plat
   use pspect
   use commap

   implicit none

#include <comsta.h>
!
! Input arguments
!
   integer, intent(in) :: lat              ! latitude index (S->N)
   integer, intent(in) :: nlon

   real(r8), intent(in) :: pint(plon,plevp)   ! pressure at model interfaces
   real(r8), intent(in) :: pdel(plon,plev)    ! pdel(k) = pint(k+1) - pint(k)
   real(r8), intent(in) :: pstar(plon)        ! ps + psr (surface pressure)
   real(r8), intent(in) :: vort(plon,plev)    ! vorticity
   real(r8), intent(in) :: div(plon,plev)     ! divergence
   real(r8), intent(in) :: t(plon,plev)       ! temperature
   real(r8), intent(in) :: q(plon,plev)       ! moisture
!
!---------------------------Local workspace-----------------------------
!
   real(r8) prat                ! pdel(i,k)/pint(i,plevp)

   integer i,k                  ! longitude, level indices
   integer ifld                 ! field index
!
!-----------------------------------------------------------------------
!
! Compute statistics for current latitude line
!
   psurf(lat) = 0._r8
   do i=1,nlon
      psurf(lat) = psurf(lat) + pstar(i)
   end do
   psurf(lat)= w(lat)*psurf(lat)/nlon

!$OMP PARALLEL DO PRIVATE (IFLD, K, I, PRAT)
   do ifld=1,4
      if     (ifld == 1) then

         rmsz (lat) = 0._r8
         do k=1,plev
            do i=1,nlon
               prat = pdel(i,k)/pint(i,plevp)
               rmsz(lat) = rmsz(lat) + vort(i,k)*vort(i,k)*prat
            end do
         end do
         rmsz(lat) = w(lat)*rmsz(lat)/nlon

      elseif (ifld == 2) then

         rmsd (lat) = 0._r8
         do k=1,plev
            do i=1,nlon
               prat = pdel(i,k)/pint(i,plevp)
               rmsd(lat) = rmsd(lat) + div(i,k)*div(i,k)*prat
            end do
         end do
         rmsd(lat) = w(lat)*rmsd(lat)/nlon

      elseif (ifld == 3) then

         rmst (lat) = 0._r8
         do k=1,plev
            do i=1,nlon
               prat = pdel(i,k)/pint(i,plevp)
               rmst(lat) = rmst(lat) + (t(i,k)**2)*prat
            end do
         end do
         rmst(lat) = w(lat)*rmst(lat)/nlon

      else

         stq  (lat) = 0._r8
         do k=1,plev
            do i=1,nlon
               prat = pdel(i,k)/pint(i,plevp)
               stq(lat) = stq(lat) + q(i,k)*pdel(i,k)
            end do
         end do
         stq (lat) = w(lat)*stq(lat)/nlon

      endif
   enddo
! 
   return
end subroutine stats
