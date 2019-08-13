
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
   use pmgrid
   use infnan,       only: nan, assignment(=)
   use constituents, only: pcnst

   implicit none
   save

   integer, parameter :: ptimelevels = 2  ! number of time levels in the dycore
   integer :: n3   = 2
   integer :: n3m1 = 1

   real(r8), allocatable, target :: ps(:,:,:)
   real(r8), allocatable, target :: u3(:,:,:,:)
   real(r8), allocatable, target :: v3(:,:,:,:)
   real(r8), allocatable, target :: t3(:,:,:,:)
   real(r8), allocatable, target :: q3(:,:,:,:,:)
   real(r8), allocatable, target :: qm1(:,:,:,:)
   real(r8), allocatable, target :: tarrsld(:,:,:)
   real(r8), allocatable, target :: parrsld(:,:,:)
   real(r8), allocatable, target :: etadot(:,:,:)

   real(r8), allocatable, target :: div(:,:,:,:)    ! divergence

   real(r8), allocatable, target :: urhs(:,:,:)
   real(r8), allocatable, target :: vrhs(:,:,:)
   real(r8), allocatable, target :: trhs(:,:,:)
   real(r8), allocatable, target :: prhs(:,:,:)

   real(r8), allocatable, target :: ql(:,:,:)
   real(r8), allocatable, target :: qm(:,:,:)
   real(r8), allocatable, target :: ed1(:,:,:)
   real(r8), allocatable, target :: tlm1(:,:,:)
   real(r8), allocatable, target :: tl(:,:,:)
   real(r8), allocatable, target :: tmm1(:,:,:)
   real(r8), allocatable, target :: tm(:,:,:)
   real(r8), allocatable, target :: omga(:,:,:)     ! vertical velocity
   real(r8), allocatable, target :: dpsl(:,:)       ! longitudinal pressure gradient
   real(r8), allocatable, target :: dpslm1(:,:)     ! longitudinal pressure gradient
   real(r8), allocatable, target :: dpslp1(:,:)     ! longitudinal pressure gradient
   real(r8), allocatable, target :: dpsm(:,:)       ! meridional pressure gradient
   real(r8), allocatable, target :: dpsmm1(:,:)     ! meridional pressure gradient
   real(r8), allocatable, target :: dpsmp1(:,:)     ! meridional pressure gradient
   real(r8), allocatable, target :: dps(:,:)        ! pressure gradient
   real(r8), allocatable, target :: phis(:,:)       ! surface geopotential
   real(r8), allocatable, target :: phisl(:,:)      ! surface geopotential
   real(r8), allocatable, target :: phism(:,:)      ! surface geopotential

CONTAINS

  subroutine initialize_prognostics
!
! Purpose: Allocate and initialize the prognostic variables
!
    allocate (ps(plon                 ,beglat:endlat,ptimelevels))
    allocate (u3(plon,plev            ,beglat:endlat,ptimelevels))
    allocate (v3(plon,plev            ,beglat:endlat,ptimelevels))
    allocate (t3(plon,plev            ,beglat:endlat,ptimelevels))
    allocate (q3(plon,plev,pcnst,beglat:endlat,ptimelevels))
    allocate (qm1(plon,plev,pcnst,beglat:endlat))
    allocate (etadot(plon,plevp,beglat:endlat))
    allocate (tarrsld(plon,plev,beglat:endlat))
    allocate (parrsld(plon,plev,beglat:endlat))

    allocate (div   (plon,plev,beglat:endlat,ptimelevels))

    allocate (urhs(plon,plev,beglat:endlat))
    allocate (vrhs(plon,plev,beglat:endlat))
    allocate (trhs(plon,plev,beglat:endlat))
    allocate (prhs(plon,plev,beglat:endlat))

    allocate (ql     (plon,plev,beglat:endlat))
    allocate (qm     (plon,plev,beglat:endlat))
    allocate (ed1    (plon,plev,beglat:endlat))
    allocate (tlm1   (plon,plev,beglat:endlat))
    allocate (tl     (plon,plev,beglat:endlat))
    allocate (tmm1   (plon,plev,beglat:endlat))
    allocate (tm     (plon,plev,beglat:endlat))
    allocate (omga   (plon,plev,beglat:endlat))
    allocate (dpsl   (plon,beglat:endlat))
    allocate (dpslm1 (plon,beglat:endlat))
    allocate (dpslp1 (plon,beglat:endlat))
    allocate (dpsm   (plon,beglat:endlat))
    allocate (dpsmm1 (plon,beglat:endlat))
    allocate (dpsmp1 (plon,beglat:endlat))
    allocate (dps    (plon,beglat:endlat))
    allocate (phis   (plon,beglat:endlat))
    allocate (phisl  (plon,beglat:endlat))
    allocate (phism  (plon,beglat:endlat))

    ps(:,:,:)     = nan
    u3(:,:,:,:)   = nan
    v3(:,:,:,:)   = nan
    t3(:,:,:,:)   = nan
    q3(:,:,:,:,:) = nan
    qm1(:,:,:,:) = nan
    tarrsld(:,:,:) = nan
    parrsld(:,:,:) = nan
    etadot(:,:,:) = nan

    div(:,:,:,:) = nan

    urhs(:,:,:) = nan
    vrhs(:,:,:) = nan
    trhs(:,:,:) = nan
    prhs(:,:,:) = nan

    ql     (:,:,:) = nan
    qm     (:,:,:) = nan
    ed1    (:,:,:) = nan
    tlm1   (:,:,:) = nan
    tl     (:,:,:) = nan
    tmm1   (:,:,:) = nan
    tm     (:,:,:) = nan
    omga   (:,:,:) = nan
    dpsl   (:,:) = nan
    dpslm1 (:,:) = nan
    dpslp1 (:,:) = nan
    dpsm   (:,:) = nan
    dpsmm1 (:,:) = nan
    dpsmp1 (:,:) = nan
    dps    (:,:) = nan
    phis   (:,:) = nan
    phisl  (:,:) = nan
    phism  (:,:) = nan

    return
  end subroutine initialize_prognostics

  subroutine shift_time_indices
!
! Purpose: Shift the time indices that track the current and previous times.
!
     integer :: itmp

     itmp = n3m1
     n3m1 = n3
     n3   = itmp
   end subroutine shift_time_indices

end module prognostics
