!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_constants
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!>
!
!-----------------------------------------------------------------------

module seaice_constants

  use mpas_derived_types
  use shr_const_mod

  private
  save

!echmod move constants in gradually
#ifdef addconstants
#ifdef CCSMCOUPLED
      real (kind=RKIND), parameter :: &
         secday    = SHR_CONST_CDAY   ,&! seconds in calendar day
         rhoi      = SHR_CONST_RHOICE ,&! density of ice (kg/m^3)
         rhow      = SHR_CONST_RHOSW  ,&! density of seawater (kg/m^3)
         cp_air    = SHR_CONST_CPDAIR ,&! specific heat of air (J/kg/K)

         ! Emissivity has been changed to unity here so that coupling is
         ! physically correct - instantaneous radiative coupling in CIME
         emissivity= 1.0_RKIND     ,&! emissivity of snow and ice

         cp_ice    = SHR_CONST_CPICE  ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = SHR_CONST_CPSW   ,&! specific heat of ocn    (J/kg/K)
                                        ! freshwater value needed for enthalpy
#ifdef RASM_MODS
         dragio    = 0.00962_RKIND ,&! ice-ocn drag coefficient for RASM as temporary measure
#else
         dragio    = 0.00536_RKIND ,&! ice-ocn drag coefficient
#endif
         gravit    = SHR_CONST_G      ,&! gravitational acceleration (m/s^2)
         rhofresh  = SHR_CONST_RHOFW  ,&! density of fresh water (kg/m^3)
         zvir      = SHR_CONST_ZVIR   ,&! rh2o/rair - 1.0
         vonkar    = SHR_CONST_KARMAN ,&! von Karman constant
         cp_wv     = SHR_CONST_CPWV   ,&! specific heat of water vapor (J/kg/K)
         stefan_boltzmann = SHR_CONST_STEBOL ,&!  W/m^2/K^4
         Tffresh   = SHR_CONST_TKFRZ  ,&! freezing temp of fresh ice (K)
         Lsub      = SHR_CONST_LATSUB ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = SHR_CONST_LATVAP ,&! latent heat, vaporization freshwater (J/kg)
         Lfresh    = SHR_CONST_LATICE ,&! latent heat of melting of fresh ice (J/kg)
         Timelt    = SHR_CONST_TKFRZ-SHR_CONST_TKFRZ ,&! melting temp. ice top surface  (C)
         Tsmelt    = SHR_CONST_TKFRZ-SHR_CONST_TKFRZ ,&! melting temp. snow top surface (C)
         ice_ref_salinity = SHR_CONST_ICE_REF_SAL ,&! (psu)
!        ocn_ref_salinity = SHR_CONST_OCN_REF_SAL ,&! (psu)

         snowpatch = 0.005_RKIND   ,&! parameter for fractional snow area (m)
         R_gC2molC = SHR_CONST_MWC    ,&! molar mass of carbon
         pi        = SHR_CONST_PI       ! pi

#else

      real (kind=RKIND), parameter, public :: &
         secday    = 86400.0_RKIND ,&! seconds in calendar day
         rhoi      = 917.0_RKIND   ,&! density of ice (kg/m^3)
         rhow      = 1026.0_RKIND  ,&! density of seawater (kg/m^3)
         cp_air    = 1005.0_RKIND  ,&! specific heat of air (J/kg/K)
         ! (Briegleb JGR 97 11475-11485  July 1992)
         emissivity= 0.95_RKIND    ,&! emissivity of snow and ice
!echmod         emissivity= 0.985_RKIND    ,&! emissivity of snow and ice
         cp_ice    = 2106._RKIND   ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = 4218._RKIND   ,&! specific heat of ocn    (J/kg/K)
                                        ! freshwater value needed for enthalpy
         dragio    = 0.00536_RKIND ,&! ice-ocn drag coefficient
         gravit    = 9.80616_RKIND ,&! gravitational acceleration (m/s^2)
         rhofresh  = 1000.0_RKIND  ,&! density of fresh water (kg/m^3)
         zvir      = 0.606_RKIND   ,&! rh2o/rair - 1.0
         vonkar    = 0.4_RKIND     ,&! von Karman constant
         cp_wv     = 1.81e3_RKIND  ,&! specific heat of water vapor (J/kg/K)
         stefan_boltzmann = 567.0e-10_RKIND,&!  W/m^2/K^4
         Tffresh   = 273.15_RKIND  ,&! freezing temp of fresh ice (K)
         Lsub      = 2.835e6_RKIND ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501e6_RKIND ,&! latent heat, vaporization freshwater (J/kg)
         Lfresh    = Lsub-Lvap        ,&! latent heat of melting of fresh ice (J/kg)
         Timelt    = 0.0_RKIND     ,&! melting temperature, ice top surface  (C)
         Tsmelt    = 0.0_RKIND     ,&! melting temperature, snow top surface (C)
         ice_ref_salinity = 4._RKIND ,&! (ppt)
         ! ocn_ref_salinity = 34.7_RKIND,&! (ppt)

         snowpatch = 0.02_RKIND, & ! parameter for fractional snow area (m)
         R_gC2molC = 12.0107_RKIND, & ! molar mass of carbon
         pi     = 3.14159265358979323846_RKIND
#endif

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=RKIND), parameter, public :: &
         rhos      = 330.0_RKIND  ,&! density of snow (kg/m^3)
         depressT  = 0.054_RKIND   ,&! Tf:brine salinity ratio (C/ppt)
         albocn    = 0.06_RKIND    ,&! ocean albedo
         viscosity_dyn = 1.79e-3_RKIND, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz   = -1.8_RKIND    ,&! freezing temp of seawater (C), used 
                                        ! as Tsfcn for open water only when 
                                        ! tfrz_option is 'minus1p8' or null

         iceruf   = 0.0005_RKIND   ,&! ice surface roughness (m)
         cprho    = cp_ocn*rhow       ,&! for ocean mixed layer (J kg / K m^3)

         ! for ice strength
         Cf       = 17._RKIND      ,&! ratio of ridging work to PE change in ridging 
         Cp       = 0.5_RKIND*gravit*(rhow-rhoi)*rhoi/rhow ,&! proport const for PE 
         Pstar    = 2.75e4_RKIND   ,&! constant in Hibler strength formula 
                                        ! (kstrength = 0) 
         Cstar    = 20._RKIND      ,&! constant in Hibler strength formula 
                                        ! (kstrength = 0) 

         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
         kappav = 1.4_RKIND ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
         !kappan = 17.6_RKIND,&! vis extnctn coef in ice, wvlngth<700nm (1/m)

         ! kice is not used for mushy thermo
         kice   = 2.03_RKIND  ,&! thermal conductivity of fresh ice(W/m/deg)
         ! kseaice is used only for zero-layer thermo
         kseaice= 2.00_RKIND  ,&! thermal conductivity of sea ice (W/m/deg)
                                   ! (used in zero layer thermodynamics option)
         ksno   = 0.30_RKIND  ,&! thermal conductivity of snow  (W/m/deg)
         zref   = 10._RKIND   ,&! reference height for stability (m)
         hs_min = 1.e-4_RKIND ,&! min snow thickness for computing zTsn (m)

         ! biogeochemistry
         sk_l = 0.03_RKIND      ! (m) skeletal layer thickness

      integer, parameter, public :: & 
         nspint = 3              ,& ! number of solar spectral intervals
         nspint_5bd = 5            ! number of solar spectral intervals used in SNICAR

      ! weights for albedos 
      ! 4 Jan 2007 BPB  Following are appropriate for complete cloud
      ! in a summer polar atmosphere with 1.5m bare sea ice surface:
      ! .636/.364 vis/nir with only 0.5% direct for each band.
      real (kind=RKIND), parameter, public :: &           ! currently used only
         awtvdr = 0.00318_RKIND, &! visible, direct  ! for history and
         awtidr = 0.00182_RKIND, &! near IR, direct  ! diagnostics
         awtvdf = 0.63282_RKIND, &! visible, diffuse
         awtidf = 0.36218_RKIND   ! near IR, diffuse

      real (kind=RKIND), parameter, public :: &
         qqqice  = 11637800._RKIND   ,&! for qsat over ice
         TTTice  = 5897.8_RKIND      ,&! for qsat over ice
         qqqocn  = 627572.4_RKIND    ,&! for qsat over ocn
         TTTocn  = 5107.4_RKIND        ! for qsat over ocn

      ! orbital parameters
      integer, public :: iyear_AD  ! Year to calculate orbit for
      real(kind=RKIND),public :: eccen  ! Earth's orbital eccentricity
      real(kind=RKIND),public :: obliqr ! Earth's obliquity in radians
      real(kind=RKIND),public :: lambm0 ! Mean longitude of perihelion at the
                                           ! vernal equinox (radians)
      real(kind=RKIND),public :: mvelpp ! Earth's moving vernal equinox longitude
                                           ! of perihelion + pi (radians)
      real(kind=RKIND),public :: obliq  ! obliquity in degrees
      real(kind=RKIND),public :: mvelp  ! moving vernal equinox long
      real(kind=RKIND),public :: decln  ! solar declination angle in radians
      real(kind=RKIND),public :: eccf   ! earth orbit eccentricity factor
      logical,public :: log_print ! Flags print of status/error

      ! snow parameters
      real (kind=RKIND), parameter, public :: &
         rhosmin   = 100.0_RKIND    ! minimum snow density (kg/m^3)

      !-----------------------------------------------------------------
      ! numbers used in column package
      !-----------------------------------------------------------------

      real (kind=RKIND), parameter, public :: &
        c0   = 0.0_RKIND, &
        c1   = 1.0_RKIND, &
        c1p5 = 1.5_RKIND, &
        c2   = 2.0_RKIND, &
        c3   = 3.0_RKIND, &
        c4   = 4.0_RKIND, &
        c5   = 5.0_RKIND, &
        c6   = 6.0_RKIND, &
        c8   = 8.0_RKIND, &
        c10  = 10.0_RKIND, &
        c15  = 15.0_RKIND, &
        c16  = 16.0_RKIND, &
        c20  = 20.0_RKIND, &
        c25  = 25.0_RKIND, &
        c100 = 100.0_RKIND, &
        c1000= 1000.0_RKIND, &
        p001 = 0.001_RKIND, &
        p01  = 0.01_RKIND, &
        p1   = 0.1_RKIND, &
        p2   = 0.2_RKIND, &
        p4   = 0.4_RKIND, &
        p5   = 0.5_RKIND, &
        p6   = 0.6_RKIND, &
        p05  = 0.05_RKIND, &
        p15  = 0.15_RKIND, &
        p25  = 0.25_RKIND, &
        p75  = 0.75_RKIND, &
        p333 = c1/c3, &
        p666 = c2/c3, &
        puny   = 1.0e-11_RKIND, &
        bignum = 1.0e+30_RKIND, &
        pih    = p5*pi
#endif

#ifdef CCSMCOUPLED
      real (kind=RKIND), parameter, public :: &
  ! fundamental constants
       seaiceSecondsPerDay  = SHR_CONST_CDAY     ! seconds in calendar day

      real (kind=RKIND), public :: &
       seaiceGravity        = SHR_CONST_G      ,&! gravitational acceleration (m/s^2)

  ! physical constants
       seaiceDensityIce     = SHR_CONST_RHOICE ,&! density of ice (kg/m^3)
       seaiceDensitySeaWater= SHR_CONST_RHOSW  ,&! density of seawater (kg/m^3)

  ! thermodynamic constants
       seaiceStefanBoltzmann           = SHR_CONST_STEBOL ,&!  W/m^2/K^4
       seaiceIceSnowEmissivity         = 1.0_RKIND        ,&! emissivity of snow and ice
       seaiceFreshWaterFreezingPoint   = SHR_CONST_TKFRZ  ,&! freezing temp of fresh ice (K)
       seaiceFreshIceSpecificHeat      = SHR_CONST_CPICE  ,&! specific heat of fresh ice (J/kg/K)
       seaiceAirSpecificHeat           = SHR_CONST_CPDAIR ,&! specific heat of air (J/kg/K)
       seaiceWaterVaporSpecificHeat    = SHR_CONST_CPWV   ,&! specific heat of water vapor (J/kg/K)
       seaiceZvir                      = SHR_CONST_ZVIR   ,&! rh2o/rair - 1.0
       seaiceLatentHeatSublimation     = SHR_CONST_LATSUB ,&! latent heat, sublimation freshwater (J/kg)
       seaiceIceSurfaceMeltingTemperature= SHR_CONST_TKFRZ-SHR_CONST_TKFRZ ,&! melting temp. ice top surface  (C)
       seaiceSnowSurfaceMeltingTemperature= SHR_CONST_TKFRZ-SHR_CONST_TKFRZ ,&! melting temp. snow top surface (C)
       seaiceVonKarmanConstant         = SHR_CONST_KARMAN ,&! von Karman constant
       seaiceSeaWaterSpecificHeat      = SHR_CONST_CPSW   ,&! specific heat of ocn    (J/kg/K)
                                                            ! freshwater value needed for enthalpy
       seaiceLatentHeatVaporization    = SHR_CONST_LATVAP ,&! latent heat, vaporization freshwater (J/kg)
       seaiceReferenceSalinity         = SHR_CONST_ICE_REF_SAL,&! ice reference salinity (ppt)
       seaiceSnowPatchiness            = 0.005_RKIND      ,&! parameter for fractional snow area (m)

#ifdef RASM_MODS
       seaiceIceOceanDragCoefficient   = 0.00962_RKIND      ! ice-ocn drag coefficient for RASM as temporary measure
#else
       seaiceIceOceanDragCoefficient   = 0.00536_RKIND      ! ice-ocn drag coefficient
#endif
!         rhofresh  = SHR_CONST_RHOFW  ,&! density of fresh water (kg/m^3)
!         Lfresh    = SHR_CONST_LATICE ,&! latent heat of melting of fresh ice (J/kg)
!         R_gC2molC = SHR_CONST_MWC    ,&! molar mass of carbon
!         pi        = SHR_CONST_PI       ! pi

#else
      real (kind=RKIND), parameter, public :: &
  ! fundamental constants
       seaiceSecondsPerDay  = 24.0_RKIND * 3600.0_RKIND
 
      real (kind=RKIND), public :: &
       seaiceGravity        = 9.80616_RKIND            ,&! gravitational acceleration (m/s^2)

  ! physical constants
       seaiceDensityIce                = 917.0_RKIND   ,&! density of ice (kg/m^3)
       seaiceDensitySeaWater           = 1026.0_RKIND  ,&! density of seawater (kg/m^3)

  ! thermodynamic constants
       seaiceStefanBoltzmann           = 567.0e-10_RKIND,&! W/m^2/K^4
         ! (Briegleb JGR 97 11475-11485  July 1992)
       seaiceIceSnowEmissivity         = 0.95_RKIND    ,&! emissivity of snow and ice
!echmod         emissivity= 0.985_RKIND
       seaiceFreshWaterFreezingPoint   = 273.15_RKIND  ,&! freezing temp of fresh ice (K)
       seaiceFreshIceSpecificHeat      = 2106._RKIND   ,&! specific heat of fresh ice (J/kg/K)
       seaiceAirSpecificHeat           = 1005.0_RKIND  ,&! specific heat of air (J/kg/K)
       seaiceWaterVaporSpecificHeat    = 1.81e3_RKIND  ,&! specific heat of water vapor (J/kg/K)
       seaiceZvir                      = 0.606_RKIND   ,&! rh2o/rair - 1.0
       seaiceLatentHeatSublimation     = 2.835e6_RKIND ,&! latent heat, sublimation freshwater (J/kg)
       seaiceIceSurfaceMeltingTemperature  = 0.0_RKIND ,&! melting temperature, ice top surface  (C)
       seaiceSnowSurfaceMeltingTemperature = 0.0_RKIND ,&! melting temperature, snow top surface  (C)
       seaiceVonKarmanConstant         = 0.4_RKIND     ,&! von Karman constant
       seaiceSeaWaterSpecificHeat      = 4218._RKIND   ,&! specific heat of ocn    (J/kg/K)
                                                         ! freshwater value needed for enthalpy
       seaiceLatentHeatVaporization    = 2.501e6_RKIND ,&! latent heat, vaporization freshwater (J/kg)
       seaiceReferenceSalinity         = 4._RKIND      ,&! ice reference salinity (ppt)
       seaiceSnowPatchiness            = 0.02_RKIND    ,&! parameter for fractional snow area (m)

       seaiceIceOceanDragCoefficient   = 0.00536_RKIND   ! ice-ocn drag coefficient

!         rhofresh  = 1000.0_RKIND  ,&! density of fresh water (kg/m^3)
!         Lfresh    = Lsub-Lvap        ,&! latent heat of melting of fresh ice (J/kg)
!         R_gC2molC = 12.0107_RKIND, & ! molar mass of carbon
!         pi     = 3.14159265358979323846_RKIND
#endif

  ! fundamental constants
  real (kind=RKIND), parameter, public :: &
       pii = 3.141592653589793_RKIND, &          ! echmod - why is this different from pi above?
       seaiceDegreesToRadians = pii / 180.0_RKIND, &
       seaiceRadiansToDegrees = 180.0_RKIND / pii, &
       seaiceSecondsPerYear = 24.0_RKIND * 3600.0_RKIND * 365.0_RKIND, &
       seaiceDaysPerSecond = 1.0_RKIND/seaiceSecondsPerDay

  real (kind=RKIND), public :: &
       seaicePi                           ! pi

  ! Earth constants
  real (kind=RKIND), parameter, public :: &
       omega         = 7.29212e-5_RKIND   ! angular rotation rate of the Earth [s-1]

  character (len=*), public, parameter :: &
       coupleAlarmID = 'coupling'

  real(kind=RKIND), public :: &
       seaicePuny

  ! physical constants
  real(kind=RKIND), public :: &
       seaiceDensitySnow, &     ! density of snow (kg/m^3)
       seaiceDensityFreshwater  ! density of freshwater (kg/m^3)

  ! thermodynamic constants
  real(kind=RKIND), public :: &
       seaiceLatentHeatMelting, &       ! latent heat of melting of fresh ice (J/kg)
       seaiceMeltingTemperatureDepression, &  ! melting temp. depression factor (C/ppt)
       seaiceOceanAlbedo, &             ! Ocean albedo
       seaiceIceSurfaceRoughness, &     ! ice surface roughness (m)
       seaiceMaximumSalinity, &         ! ice maximum salinity (ppt)
       seaiceStabilityReferenceHeight   ! stability reference height (m)

  ! dynamics constants
  real(kind=RKIND), public :: &
       seaiceIceStrengthConstantHiblerP, & ! P* constant in Hibler strength formulation
       seaiceIceStrengthConstantHiblerC    ! C* constant in Hibler strength formulation

  ! minimum sea ice area
  real(kind=RKIND), public :: &
       iceAreaMinimum, &
       iceThicknessMinimum, &
       snowThicknessMinimum

  ! biogeochemistry constants
  real(kind=RKIND), public :: &
       skeletalLayerThickness, &
       gramsCarbonPerMolCarbon   ! g carbon per mol carbon

end module seaice_constants
