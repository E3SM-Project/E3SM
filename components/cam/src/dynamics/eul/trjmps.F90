subroutine trjmps(dt      ,upr     ,vpr     ,phimp   ,lampr   , &
                  phipr   ,nlon    )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Estimate mid-point interval of parcel trajectory (global spherical
! coordinates).
! 
! Method: 
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev
!-----------------------------------------------------------------------
   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   real(r8), intent(in) :: dt                    ! time step (seconds)
   real(r8), intent(in) :: upr  (plon,plev)      ! u-comp of wind at midpoint
   real(r8), intent(in) :: vpr  (plon,plev)      ! v-comp of wind at midpoint
   real(r8), intent(in) :: phimp(plon,plev)      ! lat coord at midpoint

   integer, intent(in) :: nlon
!
! Output arguments
!
   real(r8), intent(out) :: lampr(plon,plev)      ! relative long coord of midpoint
   real(r8), intent(out) :: phipr(plon,plev)      ! relative lat  coord of midpoint
!
!-----------------------------------------------------------------------
!
!  dt      Time interval that corresponds to the parcel trajectory.
!  upr     u-coordinate of velocity corresponding to the most recent
!          estimate of the trajectory mid-point.
!  vpr     v-coordinate of velocity corresponding to the most recent
!          estimate of the trajectory mid-point.
!  phimp   Phi value of trajectory midpoint (most recent estimate).
!  lampr   Longitude coordinate of trajectory mid-point relative to the
!          arrival point.
!  phipr   Latitude  coordinate of trajectory mid-point relative to the
!          arrival point.
!
!---------------------------Local variables-----------------------------
!
   integer i,k                 ! index
!
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE (K, I)
   do k=1,plev
      do i = 1,nlon
         lampr(i,k) = -.5_r8*dt* upr(i,k) / cos( phimp(i,k) )
         phipr(i,k) = -.5_r8*dt* vpr(i,k)
      end do
   end do
!
   return
end subroutine trjmps
