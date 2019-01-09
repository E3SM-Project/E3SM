MODULE parmsld
! This module declares parameters that describe the main configuration of CRMs.
! vert_dimension & channel_dimension must be given.

!------------------------------------------------------
!  Define the number of tracer species, tie this to CAM
!------------------------------------------------------
  USE constituents, only: pcnst
  USE ppgrid,       only: pver

IMPLICIT NONE

!----------------------------------------
!  Define the parameters for vGCM points
!----------------------------------------
   INTEGER, PARAMETER :: nLevel = pver          ! The numer of vertical layers of vGCM 
   INTEGER, PARAMETER :: nhalo_vGCM = 2         ! The size of vGCM halo
   INTEGER, PARAMETER :: ntracer = pcnst - 6    ! The numer of tracers

!----------------------------------------
!  Define the parameters for VVM
!----------------------------------------   
  ! Define the vertical dimension of VVM
   INTEGER, PARAMETER :: nk1 = 30                ! The numer of vertical layers of VVM (active)

!JUNG_TEST_VERT   
!   INTEGER, PARAMETER :: nk1 = 36
   
   INTEGER, PARAMETER :: nk2 = nk1+1             ! The numer of vertical interfaces (active) 
   INTEGER, PARAMETER :: nk3 = nk1+2

   INTEGER, PARAMETER :: num_ext_layers = 6                ! number of layers above nk2   
               
!JUNG_TEST_VERT   
!   INTEGER, PARAMETER :: num_ext_layers = 0
   
   INTEGER, PARAMETER :: nk1_ext = nk1 + num_ext_layers    ! total number of layers     
   INTEGER, PARAMETER :: nk2_ext = nk1_ext + 1             ! total number of interfaces        
   
   ! Define the numer of VVM grids across the channel (not including ghost points)
   INTEGER, PARAMETER :: channel_w = 1           

   ! Define the depth of VVM halo region
   INTEGER, PARAMETER :: nhalo_cal = 1
   INTEGER, PARAMETER :: nhalo_adv = 2
   INTEGER, PARAMETER :: nhalo = nhalo_adv+1     ! The size of CRM halo

   

   INTEGER :: nVGCM_seg        ! The numer of vGCM grids in a segment
   INTEGER :: nVGCM            ! The numer of vGCM grids in a channel
   INTEGER :: netsz            ! The numer of VVM grids in a vGCM cell
   INTEGER :: netsz_sq         ! NETSZ square
   INTEGER :: channel_l        ! The numer of VVM grids along a channel
   INTEGER :: channel_seg_l    ! The numer of VVM grids along a channel-segment

   public :: parmsld_set_sizes

CONTAINS

  SUBROUTINE PARMSLD_SET_SIZES(nVGCM_seg_in, netsz_in)
    INTEGER, INTENT(IN) :: nVGCM_seg_in
    INTEGER, INTENT(IN) :: netsz_in

    nVGCM_seg     = nVGCM_seg_in       ! The numer of vGCM grids in a segment
    nVGCM         = nVGCM_seg * 4      ! The numer of vGCM grids in a channel

    netsz         = netsz_in           ! The numer of VVM grids in a vGCM cell
    netsz_sq      = netsz * netsz      ! netsz**2

    channel_l     = netsz * nVGCM      ! The numer of VVM grids along a channel
    channel_seg_l = netsz * nVGCM_seg      ! The numer of VVM grids along a channel-segment

  end subroutine parmsld_set_sizes

END MODULE parmsld
