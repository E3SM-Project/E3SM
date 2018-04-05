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
  real (kind=r8) :: THETA(channel_l,nk1)
  real (kind=r8) :: QV(channel_l,nk1)
  real (kind=r8) :: QC(channel_l,nk1)
  real (kind=r8) :: QI(channel_l,nk1)
  real (kind=r8) :: QR(channel_l,nk1)
  real (kind=r8) :: QS(channel_l,nk1)
  real (kind=r8) :: QG(channel_l,nk1)
 
  ! module physics_tendencies 
  real (kind=r8) :: TENDENCY_MICROPHYSICS_THETA(channel_l,nk1)
  real (kind=r8) :: TENDENCY_MICROPHYSICS_QV(channel_l,nk1)
  real (kind=r8) :: TENDENCY_MICROPHYSICS_QC(channel_l,nk1)
  real (kind=r8) :: TENDENCY_MICROPHYSICS_QI(channel_l,nk1)
  real (kind=r8) :: TENDENCY_MICROPHYSICS_QR(channel_l,nk1)
  real (kind=r8) :: TENDENCY_MICROPHYSICS_QS(channel_l,nk1)
  real (kind=r8) :: TENDENCY_MICROPHYSICS_QG(channel_l,nk1)
    
  real (kind=r8) :: TENDENCY_RAIN(channel_l,nk1)
  real (kind=r8) :: TENDENCY_SNOW(channel_l,nk1)
  real (kind=r8) :: TENDENCY_GRAUPEL(channel_l,nk1)
    
  real (kind=r8) :: VTR_INT(channel_l,0:nk1)
  real (kind=r8) :: VTS_INT(channel_l,0:nk1)
  real (kind=r8) :: VTG_INT(channel_l,0:nk1)
    
  real (kind=r8) :: SURFACE_RAIN(channel_l)
  real (kind=r8) :: SURFACE_SNOW(channel_l)
  real (kind=r8) :: SURFACE_GRAUPEL(channel_l)
    
  real (kind=r8) :: LATENT_HEATING_RATE(channel_l,nk1)
       
  ! module topography_parameters
  integer :: HX_SUB(channel_l)     
    
  end module mphysics_variables  
