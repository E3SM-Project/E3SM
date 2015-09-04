!! This module defines constants used by the CARMA code. Where possible, it uses constants
!! already defined in CAM.
!!
!! NOTE: CARMA constants are typically in cgs units, while CAM constants are typically in
!! mks units, so some unit conversion needs to be performed.
!!
!! NOTE: This file is adapted for use within CAM and is different than the file of the
!! same name that is part of the standard CARMA distribution.
module carma_constants_mod

  use carma_precision_mod
  use shr_const_mod,       only: SHR_CONST_TKTRIP, SHR_CONST_RHOICE, SHR_CONST_CDAY
  use physconst,           only: p_pi=>pi, avogad, boltz, r_universal, cpair, rhoh2o, latvap, latice
  use radconstants,        only: nswbands, nlwbands
  use cam_history_support, only: fillvalue

  implicit none

  !--
  ! Physical constants
  
  ! Meter-Kilogram-Second (MKS) convention for units
  ! This convention is different from CARMA's original 
  !  Centimeter-Gram-Second (CGS) convention.  Be wary of
  !  this conversion to the new convention.
  
  ! Use the _f for all literal constants, e.g. 1.2e_f.
  ! If you omit the _f in the initialization, a compiler may cast this
  !  number into single precision and then store it as _f precision.
  
  !! Define triple-point temperature (K)
  real(kind=f), parameter :: T0 = SHR_CONST_TKTRIP
  
  ! Define constants for circles and trig  
  real(kind=f), parameter :: PI = p_pi 
  real(kind=f), parameter :: DEG2RAD = PI / 180._f
  real(kind=f), parameter :: RAD2DEG = 180._f / PI
  
  !! Define avogadro's number [ # particles / mole ]
  real(kind=f), parameter :: AVG = avogad / 1000._f
   
  !! Define Boltzmann's constant [ erg / deg_K ]
  real(kind=f), parameter :: BK = boltz * 1e7_f
   
  !! Define Loschmidt's number [ mole / cm^3, @ STP ]
  real(kind=f), parameter :: ALOS = 2.68719e+19_f
  
  !! Define reference pressure, e.g. for potential temp calcs [ dyne / cm^2 ]
  real(kind=f), parameter :: PREF = 1000.e+3_f
      
  !! Define conversion factor for mb to cgs [ dyne / cm^2 ] units
  real(kind=f), parameter :: RMB2CGS = 1000.e+0_f
  
  !! Define conversion factor for Pa to cgs [ dyne / cm^2 ] units
  real(kind=f), parameter :: RPA2CGS = 10.e+0_f

  !! Define conversion factor for m to cgs [ cm ] units
  real(kind=f), parameter :: RM2CGS = 100.0_f
  
  !! Define universal gas constant [ erg / deg_K / mole ]
  real(kind=f), parameter :: RGAS = r_universal * 1e7_f / 1000._f
  
  !! Define number of seconds per the planet's day [ s / d ]
  real(kind=f), parameter :: SCDAY = SHR_CONST_CDAY
   
  !! Define specific heat at constant pres of dry air [ cm^2 / s^2 / deg_K ]
  real(kind=f), parameter :: CP = cpair * 1e7_f / 1000._f

  !! Define mass density of liquid water [ g / cm^3 ]
  real(kind=f), parameter :: RHO_W = rhoh2o / 1000._f
  
  !! Define mass density of water ice [ g / cm^3 ]
  real(kind=f), parameter :: RHO_I = SHR_CONST_RHOICE / 1000._f

  !! Latent heat of evaporation for gas [cm^2/s^2]
  real(kind=f), parameter :: RLHE_CNST = latvap * 1e4_f
  
  !! Latent heat of ice melting for gas [cm^2/s^2]
  real(kind=f), parameter :: RLHM_CNST = latice * 1e4_f
  
  
  !! The dimension of THETD, ELTRMX, CSTHT, PI, TAU, SI2THT.
  !! IT must correspond exactly to the second dimension of ELTRMX.
  integer, parameter      :: IT = 1
  
 !! String length of names 
  integer, parameter      :: CARMA_NAME_LEN = 255

  !! String length of short names 
  integer, parameter      :: CARMA_SHORT_NAME_LEN = 6
  
  !! Fill value indicating no value is being returned
  real(kind=f), parameter :: CAM_FILL = fillvalue

  
  !! Define particle number concentration [ # / cm3 ]
  !! used to decide whether to bypass microphysical processes.
  real(kind=f), parameter :: FEW_PC = 1e-6_f
  
  !! Define small particle number concentration
  !! [ # / x_units / y_units / z_units ]
  !!
  !! NOTE: Currently has mass conservation errors on the order
  !! of one part in 10^8. Keep this small enough, so that you can
  !! filter out all of the small stuff around SMALL_PC without
  !! getting neat FEW_PC.
  !!
  !! For degree/degree/hybrid coordinates, the metric is on the
  !! order of 1e20.
! real(kind=f), parameter :: SMALL_PC = 1e-50_f
  real(kind=f), parameter :: SMALL_PC = FEW_PC * 1e20 * 1e-30
  
  !!  Define core fraction (for core mass and second moment) used
  !!  when particle number concentrations are limited to SMALL_PC
  real(kind=f), parameter :: FIX_COREF = 0.1_f
    
  !! Minimum Cloud Fraction 
  real(kind=f), parameter :: CLDFRC_MIN = 0.009_f
  
  !! Incloud Cloud Fraction Threshold for statistics
  real(kind=f), parameter :: CLDFRC_INCLOUD = 0.01_f
  
  !! NWAVE should be the total number of bands CAM supports.
  integer, public, parameter      :: NWAVE           = nlwbands+nswbands ! Number of wavelength bands


  
  
  
  !! These are constants per CARMA's definition, but are set dynamically in CAM and thus
  !! can not be set as constants. They must be initialized as variables in carma_init.
  
  !! Acceleration of gravity near Earth surface [ cm/s^2 ]
  real(kind=f)            :: GRAV
  
  !! Define planet equatorial radius [ cm ]
  real(kind=f)            :: REARTH
  
  !! Define molecular weight of dry air [ g / mole ]
  real(kind=f)            :: WTMOL_AIR
  
  !! Define molecular weight of water [ g / mole ]
  real(kind=f)            :: WTMOL_H2O
  
  !! Define gas constant for dry air [ erg / deg_K / mole ]
  real(kind=f)            :: R_AIR
  
  !! Define ratio of gas constant for dry air and specific heat
  real(kind=f)            :: RKAPPA
  

end module 
