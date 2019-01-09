  Module MPHYSICS_VARIABLES
! This module declares variables used in microphysical process.
! (Used in physics_interface.f90 & physics_v10d.f90)

  USE shr_kind_mod, only: r8 => shr_kind_r8
  USE parmsld,      only: channel_l,nk1

  IMPLICIT NONE

  ! Vertical Profiles
  real (kind=r8) :: RHO_INT(0:nk1)
  real (kind=r8) :: RHOL(nk1)
  real (kind=r8) :: DZL(nk1)
  real (kind=r8) :: PL0(nk1)
  real (kind=r8) :: PIL0(nk1)

  ! module main_variables
  real (kind=r8), allocatable :: THETA(:,:)
  real (kind=r8), allocatable :: QV(:,:)
  real (kind=r8), allocatable :: QC(:,:)
  real (kind=r8), allocatable :: QI(:,:)
  real (kind=r8), allocatable :: QR(:,:)
  real (kind=r8), allocatable :: QS(:,:)
  real (kind=r8), allocatable :: QG(:,:)
  
  real (kind=r8), allocatable :: RHO_INT_loc(:,:)  ! JUNG_LocalP
  real (kind=r8), allocatable :: RHOL_loc(:,:)     ! JUNG_LocalP     
  real (kind=r8), allocatable :: PL0_loc(:,:)      ! JUNG_LocalP
  real (kind=r8), allocatable :: PIL0_loc(:,:)     ! JUNG_LocalP

  ! module physics_tendencies
  real (kind=r8), allocatable :: TENDENCY_MICROPHYSICS_THETA(:,:)
  real (kind=r8), allocatable :: TENDENCY_MICROPHYSICS_QV(:,:)
  real (kind=r8), allocatable :: TENDENCY_MICROPHYSICS_QC(:,:)
  real (kind=r8), allocatable :: TENDENCY_MICROPHYSICS_QI(:,:)
  real (kind=r8), allocatable :: TENDENCY_MICROPHYSICS_QR(:,:)
  real (kind=r8), allocatable :: TENDENCY_MICROPHYSICS_QS(:,:)
  real (kind=r8), allocatable :: TENDENCY_MICROPHYSICS_QG(:,:)

  real (kind=r8), allocatable :: TENDENCY_RAIN(:,:)
  real (kind=r8), allocatable :: TENDENCY_SNOW(:,:)
  real (kind=r8), allocatable :: TENDENCY_GRAUPEL(:,:)

  real (kind=r8), allocatable :: VTR_INT(:,:)
  real (kind=r8), allocatable :: VTS_INT(:,:)
  real (kind=r8), allocatable :: VTG_INT(:,:)

  real (kind=r8), allocatable :: SURFACE_RAIN(:)
  real (kind=r8), allocatable :: SURFACE_SNOW(:)
  real (kind=r8), allocatable :: SURFACE_GRAUPEL(:)

  real (kind=r8), allocatable :: LATENT_HEATING_RATE(:,:)

  ! module topography_parameters
  integer, allocatable        :: HX_SUB(:)

  public :: allocate_mphysics_vars

CONTAINS

!===================================================================================
  SUBROUTINE allocate_mphysics_vars ()
!===================================================================================
   
    ! module main_variables
    allocate(THETA(channel_l,nk1))
    allocate(QV(channel_l,nk1))
    allocate(QC(channel_l,nk1))
    allocate(QI(channel_l,nk1))
    allocate(QR(channel_l,nk1))
    allocate(QS(channel_l,nk1))
    allocate(QG(channel_l,nk1))
    
    allocate(RHO_INT_loc(channel_l,0:nk1))  ! JUNG_LocalP
    allocate(RHOL_loc(channel_l,nk1))       ! JUNG_LocalP     
    allocate(PL0_loc(channel_l,nk1))        ! JUNG_LocalP
    allocate(PIL0_loc(channel_l,nk1))       ! JUNG_LocalP
    
    ! module physics_tendencies
    allocate(TENDENCY_MICROPHYSICS_THETA(channel_l,nk1))
    allocate(TENDENCY_MICROPHYSICS_QV(channel_l,nk1))
    allocate(TENDENCY_MICROPHYSICS_QC(channel_l,nk1))
    allocate(TENDENCY_MICROPHYSICS_QI(channel_l,nk1))
    allocate(TENDENCY_MICROPHYSICS_QR(channel_l,nk1))
    allocate(TENDENCY_MICROPHYSICS_QS(channel_l,nk1))
    allocate(TENDENCY_MICROPHYSICS_QG(channel_l,nk1))
    
    allocate(TENDENCY_RAIN(channel_l,nk1))
    allocate(TENDENCY_SNOW(channel_l,nk1))
    allocate(TENDENCY_GRAUPEL(channel_l,nk1))
    
    allocate(VTR_INT(channel_l,0:nk1))
    allocate(VTS_INT(channel_l,0:nk1))
    allocate(VTG_INT(channel_l,0:nk1))
    
    allocate(SURFACE_RAIN(channel_l))
    allocate(SURFACE_SNOW(channel_l))
    allocate(SURFACE_GRAUPEL(channel_l))
    
    allocate(LATENT_HEATING_RATE(channel_l,nk1))
    
    ! module topography_parameters
    allocate(HX_SUB(channel_l))
    
  END SUBROUTINE allocate_mphysics_vars

  end module mphysics_variables
