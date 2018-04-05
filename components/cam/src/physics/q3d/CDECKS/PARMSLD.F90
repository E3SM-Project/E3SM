MODULE parmsld
! This module declares parameters that describe the main configuration of CRMs.
! vert_dimension & channel_dimension must be given.

!! Temp comment-out   USE ppgrid, only: pver
   
IMPLICIT NONE

!  nVGCM, ntracer, channel_seg_l (or channel_l) can be given by the main code like "pver"
 
!----------------------------------------
!  Define the number of vGCM points
!----------------------------------------
   INTEGER, PARAMETER :: nVGCM = 400             ! The numer of vGCM grids in a channel
   INTEGER, PARAMETER :: nVGCM_seg = nVGCM/4     ! The numer of vGCM grids in a segment

   INTEGER, PARAMETER :: nhalo_vGCM = 2          ! The size of vGCM halo   
   
!----------------------------------------
!  Define the number of tracer species
!----------------------------------------
   INTEGER, PARAMETER :: ntracer = 0             ! The numer of tracers

!----------------------------------------   
!  Define the vertical dimension
!----------------------------------------
!!   INTEGER, PARAMETER :: nk1 = pver              ! The numer of vertical levels (q-point) of CRM 
   INTEGER, PARAMETER :: nk1 = 20
   
   INTEGER, PARAMETER :: nk2 = nk1+1,nk3=nk1+2
   
   INTEGER, PARAMETER :: nLevel = nk1            ! The numer of vertical levels of vGCM
!----------------------------------------   
!  Define the channel segment dimension
!----------------------------------------
   INTEGER, PARAMETER :: channel_w = 1           ! The numer of CRM grids across the channel 
                                                 ! (including only prediction points)
   
   INTEGER, PARAMETER :: channel_seg_l = 2500    ! The numer of CRM grids along a channel-segment 
                                                 ! (Temporary setting for now) 
   
   INTEGER, PARAMETER :: channel_l = channel_seg_l*4    ! The numer of CRM grids along a channel
   
   INTEGER, PARAMETER :: netsz = channel_l/nVGCM        ! The numer of CRM grids in a vGCM cell
   INTEGER, PARAMETER :: netsz_sq = netsz*netsz         ! NETSZ square
!----------------------------------------
!  Define the depth of the CRM halo region
!----------------------------------------
   INTEGER, PARAMETER :: nhalo_cal = 1
   INTEGER, PARAMETER :: nhalo_adv = 2 
   INTEGER, PARAMETER :: nhalo = nhalo_adv+1     ! The size of CRM halo

END MODULE parmsld
