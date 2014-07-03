
!----------------------------------------------------------------------- 
!BOP
! !ROUTINE:  initcom --- Initialize model commons
!
! !INTERFACE:
subroutine initcom

! !USES:
   use shr_kind_mod,  only: r8 => shr_kind_r8
   use physconst,     only: pi, rair
   use pmgrid
   use pspect
   use commap,        only: w, clat, clon, w_staggered, clat
   use commap,        only: clat_staggered, clon, latdeg, londeg
   use commap,        only: latdeg_st, londeg_st
   use abortutils,    only: endrun
   use rgrid,         only: fullgrid
   use spmd_utils,    only: masterproc
   use cam_logfile,   only: iulog
   implicit none

!-----------------------------------------------------------------------

! !DESCRIPTION:
!
!   Initialize Model commons.
! 
! !REVISION HISTORY:
!
!   92.06.01      Bath          Creation from CCM1
!   96.02.01      Buja          Modifications
!   01.01.19      Lin           incorporated Rasch's bug fix for the 
!                               "Gaussian" weights
!   02.04.04      Sawyer        Removed comspe
!
!EOP
!-----------------------------------------------------------------------
!BOC
!
! !LOCAL VARIABLES:

   real(r8), parameter ::  D0_0                     =  0.0_r8
   real(r8), parameter ::  D1_0                     =  1.0_r8
   real(r8), parameter ::  D1_5                     =  1.5_r8
   real(r8), parameter ::  D2_0                     =  2.0_r8
   real(r8), parameter ::  D90_0                    =  90.0_r8
   real(r8), parameter ::  D180_0                   = 180.0_r8
   real(r8), parameter ::  D360_0                   = 360.0_r8
   real(r8), parameter ::  D1EM8                    =  1.0e-8_r8

   integer i           ! longitude index
   integer j           ! Latitude index
   integer k           ! Level index
   integer m           ! Index for legendre array
   integer irow        ! Latitude pair index
   integer lat         ! Latitude index

   real(r8) dp             ! Spacing between latitudes
   real(r8) sum
!
!-----------------------------------------------------------------------

! L-R dynamics uses a regular latitude distribution (not gausian).
! The algorithm below is a bastardized version of LSM: map.F.

   dp = D180_0/(plat-1)
   do lat = 1, plat
      latdeg(lat) = -D90_0 + (lat-1)*dp
      clat(lat) = latdeg(lat)*pi/D180_0
   end do

! Calculate latitudes for the staggered grid

   do lat = 1, plat-1
      clat_staggered(lat) = (clat(lat) + clat(lat+1)) / D2_0
      latdeg_st     (lat) = clat_staggered(lat)*D180_0/pi
   end do

! Weights are defined as cos(phi)*(delta-phi)
! For a sanity check, the sum of w across all lats should be 2, or 1 across
! half of the latitudes.

   do lat = 2, plat-1
      w(lat) = sin(clat_staggered(lat)) - sin(clat_staggered(lat-1))
   end do
   w(1) = sin(clat_staggered(1)) + D1_0
   w(plat) = w(1)

   sum = D0_0
   do lat=1,plat
      if (masterproc) write(iulog,*) 'initcom: lat, clat, w ', lat, clat(lat), w(lat)
      sum = sum + w(lat)
   end do

   if (abs(sum - D2_0) > D1EM8) then
      write(iulog,*) 'INITCOM 1: weights do not sum to 2. sum=',sum
      call endrun
   end if

   dp = pi / real(plat-1,r8)
   do lat = 1, plat-1
      w_staggered(lat) = sin(clat(lat+1)) - sin(clat(lat))
   end do

   sum = D0_0
   do lat=1,plat-1
      sum = sum + w_staggered(lat)
   end do

   if (abs(sum - D2_0) > D1EM8) then
      write(iulog,*) 'INITCOM 2: weights do not sum to 2. sum=',sum
      call endrun
   end if
!
! Determine whether full or reduced grid
!
   fullgrid = .true.
   if (masterproc) write(iulog,*) 'Number of longitudes per latitude = ', plon
!
! Longitude array
!
   do lat=1,plat
      do i=1,plon
         londeg(i,lat) = (i-1)*D360_0/plon
         clon(i,lat)   = (i-1)*D2_0*pi/plon
      end do
   end do
   do lat=1,plat
      do i=1,splon
         londeg_st(i,lat) = (i-D1_5)*D360_0/splon
      end do
   end do
!
! Set flag indicating dynamics grid is now defined.
! NOTE: this ASSUMES initcom is called after spmdinit.  The setting of nlon done here completes
! the definition of the dynamics grid.
!
   return
!EOC
end subroutine initcom
!-----------------------------------------------------------------------
