  Module MPHYSICS_VARIABLES
! This module declares variables used in microphysical process. 
! (Used in physics_interface.f90 & physics_v10d.f90)

  USE shr_kind_mod, only: dbl_kind => shr_kind_r8
  USE parmsld,      only: channel_l,nk1
  
  IMPLICIT NONE

  ! Vertical Profiles
  real (kind=dbl_kind) :: RHO_INT(0:nk1)
  real (kind=dbl_kind) :: RHOL(nk1)
  real (kind=dbl_kind) :: DZL(nk1)
  real (kind=dbl_kind) :: PL0(nk1)
  real (kind=dbl_kind) :: PIL0(nk1)
         
  ! module main_variables
  real (kind=dbl_kind) :: THETA(channel_l,nk1)
  real (kind=dbl_kind) :: QV(channel_l,nk1)
  real (kind=dbl_kind) :: QC(channel_l,nk1)
  real (kind=dbl_kind) :: QI(channel_l,nk1)
  real (kind=dbl_kind) :: QR(channel_l,nk1)
  real (kind=dbl_kind) :: QS(channel_l,nk1)
  real (kind=dbl_kind) :: QG(channel_l,nk1)
 
  ! module physics_tendencies 
  real (kind=dbl_kind) :: TENDENCY_MICROPHYSICS_THETA(channel_l,nk1)
  real (kind=dbl_kind) :: TENDENCY_MICROPHYSICS_QV(channel_l,nk1)
  real (kind=dbl_kind) :: TENDENCY_MICROPHYSICS_QC(channel_l,nk1)
  real (kind=dbl_kind) :: TENDENCY_MICROPHYSICS_QI(channel_l,nk1)
  real (kind=dbl_kind) :: TENDENCY_MICROPHYSICS_QR(channel_l,nk1)
  real (kind=dbl_kind) :: TENDENCY_MICROPHYSICS_QS(channel_l,nk1)
  real (kind=dbl_kind) :: TENDENCY_MICROPHYSICS_QG(channel_l,nk1)
    
  real (kind=dbl_kind) :: TENDENCY_RAIN(channel_l,nk1)
  real (kind=dbl_kind) :: TENDENCY_SNOW(channel_l,nk1)
  real (kind=dbl_kind) :: TENDENCY_GRAUPEL(channel_l,nk1)
    
  real (kind=dbl_kind) :: VTR_INT(channel_l,0:nk1)
  real (kind=dbl_kind) :: VTS_INT(channel_l,0:nk1)
  real (kind=dbl_kind) :: VTG_INT(channel_l,0:nk1)
    
  real (kind=dbl_kind) :: SURFACE_RAIN(channel_l)
  real (kind=dbl_kind) :: SURFACE_SNOW(channel_l)
  real (kind=dbl_kind) :: SURFACE_GRAUPEL(channel_l)
    
  real (kind=dbl_kind) :: LATENT_HEATING_RATE(channel_l,nk1)
       
  ! module topography_parameters
  integer, pointer :: HX_SUB(channel_l)     
    
  end module mphysics_variables  
