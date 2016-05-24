module carma_constants_mod

  use carma_precision_mod

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
	real(kind=f), parameter :: T0 = 273.16_f
	
	! Define constants for circles and trig  
	real(kind=f), parameter :: PI = 3.14159265358979_f 
	real(kind=f), parameter :: DEG2RAD = PI / 180._f
	real(kind=f), parameter :: RAD2DEG = 180._f / PI
	
	!! Acceleration of gravity near Earth surface [ cm/s^2 ]
	real(kind=f), parameter :: GRAV = 980.6_f 
	
	!! Define planet equatorial radius [ cm ]
	real(kind=f), parameter :: REARTH  = 6.37e+8_f 
	
	!! Define avogadro's number [ # particles / mole ]
	real(kind=f), parameter :: AVG = 6.02252e+23_f
	 
	!! Define Boltzmann's constant [ erg / deg_K ]
	real(kind=f), parameter :: BK = 1.38054e-16_f 
	 
	!! Define Loschmidt's number [ mole / cm^3, @ STP ]
	real(kind=f), parameter :: ALOS = 2.68719e+19_f
	
	!! Define molecular weight of dry air [ g / mole ]
	real(kind=f), parameter :: WTMOL_AIR = 28.966_f
	
	!! Define molecular weight of water vaor [ g / mole ]
	real(kind=f), parameter :: WTMOL_H2O = 18.016_f
	
  !! Define reference pressure, e.g. for potential temp calcs [ dyne / cm^2 ]
  real(kind=f), parameter :: PREF = 1000.e+3_f
      
	!! Define conversion factor for mb to cgs [ dyne / cm^2 ] units
	real(kind=f), parameter :: RMB2CGS = 1000.e+0_f
	
	!! Define conversion factor for Pa to cgs [ dyne / cm^2 ] units
	real(kind=f), parameter :: RPA2CGS = 10.e+0_f

	!! Define conversion factor for m to cgs [ cm ] units
	real(kind=f), parameter :: RM2CGS = 100.0_f
	
	!! Define universal gas constant [ erg / deg_K / mole ]
	real(kind=f), parameter :: RGAS = 8.31430e+07_f 
	
	!! Define gas constant for dry air [ erg / deg_K / mole ]
	real(kind=f), parameter :: R_AIR = RGAS / WTMOL_AIR
	
	!! Define number of seconds per the planet's day [ s / d ]
	real(kind=f), parameter :: SCDAY = 86400._f
	 
	!! Define specific heat at constant pres of dry air [ cm^2 / s^2 / deg_K ]
	real(kind=f), parameter :: CP = 1.004e7_f 
	 
	!! Define ratio of gas constant for dry air and specific heat
	real(kind=f), parameter :: RKAPPA = R_AIR / CP
	
	!! Define mass density of liquid water [ g / cm^3 ]
	real(kind=f), parameter :: RHO_W = 1._f
	
	!! Define mass density of water ice [ g / cm^3 ]
	real(kind=f), parameter :: RHO_I = 0.93_f

	!! Latent heat of evaporation for gas [cm^2/s^2]
	real(kind=f), parameter :: RLHE_CNST = 2.501e10_f
	
	!! Latent heat of ice melting for gas [cm^2/s^2]
	real(kind=f), parameter :: RLHM_CNST = 3.337e9_f
	
	!! The dimension of THETD, ELTRMX, CSTHT, PI, TAU, SI2THT.
  !! IT must correspond exactly to the second dimension of ELTRMX.
  integer, parameter      :: IT = 1
  
  !! String length of names 
  integer, parameter      ::  CARMA_NAME_LEN = 255

  !! String length of short names 
  integer, parameter      ::  CARMA_SHORT_NAME_LEN = 6
  
  !! Fill value indicating no value is being returned
  integer, parameter      :: CAM_FILL = -999

	
	!! Define small particle number concentration
	!! [ # / x_units / y_units / z_units ]
	real(kind=f), parameter :: SMALL_PC = 1e-50_f
!	real(kind=f), parameter :: SMALL_PC = tiny( ONE )
	
	!!  Define particle number concentration [ # / ? ]
	!!  used to decide whether to bypass microphysical processes.
	!!  
	!!  Set it to SMALL_PC/xmet/ymet/zmet to never bypass the calculations.
	
	real(kind=f), parameter :: FEW_PC = SMALL_PC * 1e6_f
!	real(kind=f), parameter :: FEW_PC = tiny(ONE) * 1e6_f
	
	!!  Define core fraction (for core mass and second moment) used
	!!  when particle number concentrations are limited to SMALL_PC
	real(kind=f), parameter :: FIX_COREF = 0.1_f
	
	!! Minimum Cloud Fraction 
	real(kind=f), parameter :: CLDFRC_MIN = 1e-4_f
	
	!! Incloud Cloud Fraction Threshold for statistics
	real(kind=f), parameter :: CLDFRC_INCLOUD = 0.10_f
end module 
