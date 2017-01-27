! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module module_ra_rrtmg_vinterp
 use mpas_kind_types

 implicit none
 private
 public:: vinterp_ozn


 contains


!=================================================================================================================
 subroutine vinterp_ozn (lchnk, ncol, pcols, pver, pmid_in, pin_in, levsiz, ozmix_in, o3vmr_out)
!--------------------------------------------------------------------------------------------------
! 
! Purpose: Interpolate ozone from current time-interpolated values to model levels
! 
! Method: Use pressure values to determine interpolation levels
! 
! Author: Bruce Briegleb
! 
!-----------------------------------------------------------------------------------------------------------------
!  use shr_kind_mod, only: r8 => shr_kind_r8
!  use ppgrid
!  use phys_grid,     only: get_lat_all_p, get_lon_all_p
!  use comozp
!  use abortutils, only: endrun
!-----------------------------------------------------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------------------------------------------------
!
! Arguments
!
   integer, intent(in) :: lchnk               ! chunk identifier
   integer, intent(in) :: pcols, pver
   integer, intent(in) :: ncol                ! number of atmospheric columns
   integer, intent(in) :: levsiz              ! number of ozone layers

!ldf begin:
!  real(kind=RKIND), intent(in) :: pmid(pcols,pver)   ! level pressures (mks)
!  real(kind=RKIND), intent(in) :: pin(levsiz)        ! ozone data level pressures (mks)
!  real(kind=RKIND), intent(in) :: ozmix(pcols,levsiz)! ozone mixing ratio

!  real(kind=RKIND), intent(out) :: o3vmr(pcols,pver) ! ozone volume mixing ratio

   real(kind=RKIND), intent(in) :: pmid_in(pcols,pver)   ! level pressures (mks)
   real(kind=RKIND), intent(in) :: pin_in(levsiz)        ! ozone data level pressures (mks)
   real(kind=RKIND), intent(in) :: ozmix_in(pcols,levsiz)! ozone mixing ratio

   real(kind=RKIND), intent(out) :: o3vmr_out(pcols,pver)! ozone volume mixing ratio
!ldf end.
!
! local storage
!
   integer i                   ! longitude index
   integer k, kk, kkstart      ! level indices
   integer kupper(pcols)       ! Level indices for interpolation
   integer kount               ! Counter
   integer lats(pcols)         ! latitude indices
   integer lons(pcols)         ! latitude indices

   real(kind=RKIND) dpu        ! upper level pressure difference
   real(kind=RKIND) dpl        ! lower level pressure difference

!ldf begin:
   real(kind=RKIND):: pmid(pcols,pver)                ! level pressures (mks)
   real(kind=RKIND):: pin(levsiz)                     ! ozone data level pressures (mks)
   real(kind=RKIND):: ozmix(pcols,levsiz)             ! ozone mixing ratio
   real(kind=RKIND):: o3vmr(pcols,pver) ! ozone volume mixing ratio
!
! Initialize latitude indices
!
!  call get_lat_all_p(lchnk, ncol, lats)
!  call get_lon_all_p(lchnk, ncol, lons)
!

!ldf begin:
   do k = 1,levsiz
      pin(k) = pin_in(k)
   enddo
   do i = 1,pcols
      do k = 1,levsiz
         ozmix(i,k) = ozmix_in(i,k)
      enddo
   enddo
   do i = 1,pcols
      do k = 1,pver
         kk = pver-k+1
         pmid(i,kk) = pmid_in(i,k)
      enddo
   enddo
!ldf end.

! Initialize index array
!
   do i=1,ncol
      kupper(i) = 1
   end do

   do k=1,pver
!
! Top level we need to start looking is the top level for the previous k
! for all longitude points
!
      kkstart = levsiz
      do i=1,ncol
         kkstart = min0(kkstart,kupper(i))
      end do
      kount = 0
!
! Store level indices for interpolation
!
      do kk=kkstart,levsiz-1
         do i=1,ncol
            if (pin(kk).lt.pmid(i,k) .and. pmid(i,k).le.pin(kk+1)) then
               kupper(i) = kk
               kount = kount + 1
            end if
         end do
!
! If all indices for this level have been found, do the interpolation and
! go to the next level
!
         if (kount.eq.ncol) then
            do i=1,ncol
               dpu = pmid(i,k) - pin(kupper(i))
               dpl = pin(kupper(i)+1) - pmid(i,k)
               o3vmr(i,k) = (ozmix(i,kupper(i))*dpl + &
                             ozmix(i,kupper(i)+1)*dpu)/(dpl + dpu)
            end do
            goto 35
         end if
      end do
!
! If we've fallen through the kk=1,levsiz-1 loop, we cannot interpolate and
! must extrapolate from the bottom or top ozone data level for at least some
! of the longitude points.
!
      do i=1,ncol
         if (pmid(i,k) .lt. pin(1)) then
            o3vmr(i,k) = ozmix(i,1)*pmid(i,k)/pin(1)
         else if (pmid(i,k) .gt. pin(levsiz)) then
            o3vmr(i,k) = ozmix(i,levsiz)
         else
            dpu = pmid(i,k) - pin(kupper(i))
            dpl = pin(kupper(i)+1) - pmid(i,k)
            o3vmr(i,k) = (ozmix(i,kupper(i))*dpl + &
                          ozmix(i,kupper(i)+1)*dpu)/(dpl + dpu)
         end if
      end do

      if (kount.gt.ncol) then
!        call endrun ('VINTERP_OZN: Bad ozone data: non-monotonicity suspected')
      end if
35    continue
   end do

!ldf begin:
   do i = 1,pcols
   do k = 1,pver
      kk = pver-k+1
      o3vmr_out(i,kk) = o3vmr(i,k)
   enddo
   enddo
!ldf end.

   return
end subroutine vinterp_ozn

!=================================================================================================================
 end module module_ra_rrtmg_vinterp
!=================================================================================================================
