!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  seaice_constants
!
!> \brief
!> \author Adrian K. Turner, LANL
!> \date 2013-2014
!> \details
!> Values in this module have been copied from the two ice_constants_colpkg.F90 
!  modules in the column physics package (colpkg, in column/).  They will remain 
!  duplicates until the column physic package is deprecated.
!-----------------------------------------------------------------------

module seaice_constants

  use mpas_derived_types
  use shr_const_mod

  private
  save

!echmod move constants in gradually
#ifdef addconstants

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=RKIND), parameter, public :: &
         viscosity_dyn = 1.79e-3_RKIND, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz   = -1.8_RKIND    ,&! freezing temp of seawater (C), used 
                                        ! as Tsfcn for open water only when 
                                        ! tfrz_option is 'minus1p8' or null

         cprho    = cp_ocn*rhow       ,&! for ocean mixed layer (J kg / K m^3)

         ! for ice strength
         Cf       = 17._RKIND      ,&! ratio of ridging work to PE change in ridging 
         Cp       = 0.5_RKIND*gravit*(rhow-rhoi)*rhoi/rhow ,&! proport const for PE 

         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
         kappav = 1.4_RKIND ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
         !kappan = 17.6_RKIND,&! vis extnctn coef in ice, wvlngth<700nm (1/m)

         ! kice is not used for mushy thermo
         kice   = 2.03_RKIND  ,&! thermal conductivity of fresh ice(W/m/deg)
         ! kseaice is used only for zero-layer thermo
         kseaice= 2.00_RKIND  ,&! thermal conductivity of sea ice (W/m/deg)
                                   ! (used in zero layer thermodynamics option)
         ksno   = 0.30_RKIND  ,&! thermal conductivity of snow  (W/m/deg)
         hs_min = 1.e-4_RKIND ,&! min snow thickness for computing zTsn (m)

         ! biogeochemistry
!         sk_l = 0.03_RKIND      ! (m) skeletal layer thickness

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
        puny   
        bignum = 1.0e+30_RKIND, &
        pih    = p5*pi
#endif

#ifdef CCSMCOUPLED
      real (kind=RKIND), parameter, public :: &
  ! fundamental constants
       pi                              = SHR_CONST_PI     ,&! pi
       seaiceSecondsPerDay             = SHR_CONST_CDAY     ! seconds in calendar day

      real (kind=RKIND), public :: &
       seaiceGravity                   = SHR_CONST_G      ,&! gravitational acceleration (m/s^2)

  ! physical constants
       seaiceDensityIce                = SHR_CONST_RHOICE ,&! density of ice (kg/m^3)
       seaiceDensitySeaWater           = SHR_CONST_RHOSW  ,&! density of seawater (kg/m^3)
       seaiceDensityFreshwater         = SHR_CONST_RHOFW  ,&! density of fresh water (kg/m^3)

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
       seaiceLatentHeatMelting         = SHR_CONST_LATICE ,&! latent heat of melting of fresh ice (J/kg)
       seaiceReferenceSalinity         = SHR_CONST_ICE_REF_SAL,&! ice reference salinity (ppt)
       seaiceSnowPatchiness            = 0.005_RKIND      ,&! parameter for fractional snow area (m)

#ifdef RASM_MODS
       seaiceIceOceanDragCoefficient   = 0.00962_RKIND      ! ice-ocn drag coefficient for RASM as temporary measure
#else
       seaiceIceOceanDragCoefficient   = 0.00536_RKIND      ! ice-ocn drag coefficient
#endif

!         R_gC2molC = SHR_CONST_MWC    ,&! molar mass of carbon

#else
      real (kind=RKIND), parameter, public :: &
  ! fundamental constants
       pi                              = 3.14159265358979323846_RKIND
       seaiceSecondsPerDay             = 24.0_RKIND * 3600.0_RKIND
 
      real (kind=RKIND), public :: &
       seaiceGravity                   = 9.80616_RKIND ,&! gravitational acceleration (m/s^2)

  ! physical constants
       seaiceDensityIce                = 917.0_RKIND   ,&! density of ice (kg/m^3)
       seaiceDensitySeaWater           = 1026.0_RKIND  ,&! density of seawater (kg/m^3)
       seaiceDensityFreshwater         = 1000.0_RKIND  ,&! density of fresh water (kg/m^3)

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
       seaiceLatentHeatMelting         = seaiceLatentHeatSublimation & ! latent heat of melting of fresh ice (J/kg)
                                       - seaiceLatentHeatVaporization, &
       seaiceReferenceSalinity         = 4._RKIND      ,&! ice reference salinity (ppt)
       seaiceSnowPatchiness            = 0.02_RKIND    ,&! parameter for fractional snow area (m)

       seaiceIceOceanDragCoefficient   = 0.00536_RKIND   ! ice-ocn drag coefficient

!         R_gC2molC = 12.0107_RKIND, & ! molar mass of carbon
#endif

  ! fundamental constants
  real (kind=RKIND), parameter, public :: &
       pii = 3.141592653589793_RKIND, & ! echmod - why is this different from pi above?
       seaiceDegreesToRadians = pii / 180.0_RKIND, &
       seaiceRadiansToDegrees = 180.0_RKIND / pii, &
!echmod - BFB but 20% slower???
!       pii = pi, &
!       seaiceDegreesToRadians = pi / 180.0_RKIND, &
!       seaiceRadiansToDegrees = 180.0_RKIND / pi, &
       seaiceSecondsPerYear = 24.0_RKIND * 3600.0_RKIND * 365.0_RKIND, & !echmod - incorrect for leap years
       seaiceDaysPerSecond = 1.0_RKIND/seaiceSecondsPerDay

  real (kind=RKIND), public :: &
       seaicePi
!echmod - BFB but 20% slower???
!       seaicePi                           = pi ! pi

  ! Earth constants
  real (kind=RKIND), parameter, public :: &
       omega                              = 7.29212e-5_RKIND ! angular rotation rate of the Earth [s-1]

  character (len=*), public, parameter :: &
       coupleAlarmID = 'coupling'

  real(kind=RKIND), public :: &
       seaicePuny                         = 1.0e-11_RKIND   ! a small value

  ! physical constants
  real(kind=RKIND), public :: &
       seaiceDensitySnow                  = 330.0_RKIND     ! density of snow (kg/m^3)

  ! thermodynamic constants
  real(kind=RKIND), public :: &
       seaiceMeltingTemperatureDepression = 0.054_RKIND   ,&! melting temp. depression factor (C/ppt)
       seaiceOceanAlbedo                  = 0.06_RKIND    ,&! ocean albedo
       seaiceIceSurfaceRoughness          = 0.0005_RKIND  ,&! ice surface roughness (m)
       seaiceMaximumSalinity              = 3.2_RKIND     ,&! max salinity at ice base for BL99 (ppt)
       seaiceStabilityReferenceHeight     = 10._RKIND       ! reference height for stability (m)

  ! dynamics constants
  real(kind=RKIND), public :: &
       seaiceIceStrengthConstantHiblerP   = 2.75e4_RKIND  ,&! P* constant in Hibler strength formula
       seaiceIceStrengthConstantHiblerC   = 20._RKIND       ! C* constant in Hibler strength formula

  ! minimum sea ice area
  real(kind=RKIND), public :: &
       iceAreaMinimum, &
       iceThicknessMinimum, &
       snowThicknessMinimum

  ! biogeochemistry constants
  real(kind=RKIND), public :: &
       skeletalLayerThickness             = 0.03_RKIND    ,&! (m) skeletal layer thickness
       gramsCarbonPerMolCarbon   ! g carbon per mol carbon

   ! ocean biogeochemistry ISPOL values
   real(kind=RKIND), parameter, public :: &
         oceanAmmoniumISPOL        = 1.0_RKIND, & ! mmol N m-3
         oceanDMSISPOL             = 0.1_RKIND, & ! mmol S m-3
         oceanDMSPISPOL            = 0.1_RKIND, & ! mmol S m-3
         oceanDiatomsISPOL         = 1.0_RKIND, & ! mmol N m-3
         oceanSmallAlgaeISPOL      = 0.0057_RKIND, & ! mmol N m-3
         oceanPhaeocystisISPOL     = 0.0027_RKIND, & ! mmol N m-3
         oceanPolysaccharidsISPOL  = 16.2_RKIND, & ! mmol C m-3
         oceanLipidsISPOL          = 9.0_RKIND, & ! mmol C m-3
         oceanProteinsCarbonISPOL  = 9.0_RKIND, & ! mmol C m-3
         oceanDICISPOL             = 1.0_RKIND, & ! mmol C m-3
         oceanProteinsISPOL        = 12.9_RKIND, & ! mmol N m-3
         oceanDissolvedIronISPOL   = 0.4_RKIND, & ! mmol Fe m-3
         oceanParticulateIronISPOL = 2.0_RKIND,& ! mmol Fe m-3
         oceanHumicsISPOL          = 1.0_RKIND ! mmol C m-3

end module seaice_constants
