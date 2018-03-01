
module prognostics

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Prognostic variables held in-core for convenient access.
! q3 is specific humidity (water vapor) and other constituents.
! 
! Author: G. Grant
! 
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, beglat, endlat
   use infnan,       only: posinf, assignment(=)
   use constituents, only: pcnst


   implicit none

   private

   public ps, u3, v3, t3, q3, qminus, vort, div, dpsl, dpsm, dps, omga, phis, hadv, pdeld
   public n3, n3m1, n3m2, ptimelevels
   public initialize_prognostics
   public shift_time_indices

   integer, parameter :: ptimelevels = 3  ! number of time levels in the dycore
   integer :: n3   = 3
   integer :: n3m1 = 2
   integer :: n3m2 = 1

   real(r8), allocatable, target :: ps(:,:,:)
   real(r8), allocatable, target :: u3(:,:,:,:)
   real(r8), allocatable, target :: v3(:,:,:,:)
   real(r8), allocatable, target :: t3(:,:,:,:)
   real(r8), allocatable, target :: pdeld(:,:,:,:)
   real(r8), allocatable, target :: q3(:,:,:,:,:)
   real(r8), allocatable :: qminus(:,:,:,:)
   real(r8), allocatable :: hadv  (:,:,:,:)

   real(r8), allocatable, target :: vort(:,:,:,:)   ! vorticity
   real(r8), allocatable, target :: div(:,:,:,:)    ! divergence

   real(r8), allocatable, target :: dpsl(:,:)       ! longitudinal pressure gradient
   real(r8), allocatable, target :: dpsm(:,:)       ! meridional pressure gradient
   real(r8), allocatable, target :: dps(:,:)        ! pressure gradient
   real(r8), allocatable, target :: phis(:,:)       ! surface geopotential
   real(r8), allocatable, target :: omga(:,:,:)     ! vertical velocity

CONTAINS

   subroutine initialize_prognostics
!
! Purpose:  Allocate and initialize the prognostic arrays.
!

      allocate (ps    (plon           ,beglat:endlat    ,ptimelevels))
      allocate (u3    (plon,plev      ,beglat:endlat,ptimelevels))
      allocate (v3    (plon,plev      ,beglat:endlat,ptimelevels))
      allocate (t3    (plon,plev      ,beglat:endlat,ptimelevels))
      allocate (q3    (plon,plev,pcnst,beglat:endlat,ptimelevels))
      allocate (qminus(plon,plev,pcnst,beglat:endlat  ))
      allocate (hadv  (plon,plev,pcnst,beglat:endlat  ))

      allocate (vort  (plon,plev,beglat:endlat,ptimelevels))   
      allocate (div   (plon,plev,beglat:endlat,ptimelevels))    

      allocate (dpsl  (plon,beglat:endlat))        
      allocate (dpsm  (plon,beglat:endlat))        
      allocate (dps   (plon,beglat:endlat))         
      allocate (phis  (plon,beglat:endlat))        
      allocate (omga  (plon,plev,beglat:endlat))    
      allocate (pdeld (plon,plev,beglat:endlat,ptimelevels))

      ps(:,:,:)       = posinf
      u3(:,:,:,:)     = posinf
      v3(:,:,:,:)     = posinf
      t3(:,:,:,:)     = posinf
      pdeld(:,:,:,:)  = posinf
      q3(:,:,:,:,:)   = posinf
      qminus(:,:,:,:) = posinf
      hadv  (:,:,:,:) = posinf

      vort(:,:,:,:)   = posinf
      div (:,:,:,:)   = posinf

      dpsl  (:,:) = posinf
      dpsm  (:,:) = posinf
      dps   (:,:) = posinf
      phis  (:,:) = posinf
      omga  (:,:,:) = posinf

      return
   end subroutine initialize_prognostics

   subroutine shift_time_indices
!
! Purpose: 
! Shift the indices that keep track of which index stores
! the relative times (current time, previous, time before previous etc).
!
      integer :: itmp

      itmp = n3m2

      n3m2 = n3m1
      n3m1 = n3
      n3   = itmp
   end subroutine shift_time_indices

end module prognostics
