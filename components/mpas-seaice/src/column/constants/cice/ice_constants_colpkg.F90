!  SVN:$Id: ice_constants_colpkg.F90 1012 2015-06-26 12:34:09Z eclare $
!=======================================================================
!
! This module defines a variety of physical and numerical constants
! used in the column package
!
! author Elizabeth C. Hunke, LANL

      module ice_constants_colpkg

      use ice_kinds_mod

      implicit none
      save
      private

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         secday    = 86400.0_dbl_kind ,&! seconds in calendar day
         rhos      = 330.0_dbl_kind   ,&! density of snow (kg/m^3)
         rhoi      = 917.0_dbl_kind   ,&! density of ice (kg/m^3)
         rhow      = 1026.0_dbl_kind  ,&! density of seawater (kg/m^3)
         cp_air    = 1005.0_dbl_kind  ,&! specific heat of air (J/kg/K)
         ! (Briegleb JGR 97 11475-11485  July 1992)
         emissivity= 0.95_dbl_kind    ,&! emissivity of snow and ice
!echmod         emissivity= 0.985_dbl_kind    ,&! emissivity of snow and ice
         cp_ice    = 2106._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = 4218._dbl_kind   ,&! specific heat of ocn    (J/kg/K)
                                        ! freshwater value needed for enthalpy
         depressT  = 0.054_dbl_kind   ,&! Tf:brine salinity ratio (C/ppt)
         dragio    = 0.00536_dbl_kind ,&! ice-ocn drag coefficient
         albocn    = 0.06_dbl_kind    ,&! ocean albedo
         gravit    = 9.80616_dbl_kind ,&! gravitational acceleration (m/s^2)
         viscosity_dyn = 1.79e-3_dbl_kind, & ! dynamic viscosity of brine (kg/m/s)
         Tocnfrz   = -1.8_dbl_kind    ,&! freezing temp of seawater (C), used 
                                        ! as Tsfcn for open water only when 
                                        ! tfrz_option is 'minus1p8' or null
         rhofresh  = 1000.0_dbl_kind  ,&! density of fresh water (kg/m^3)
         zvir      = 0.606_dbl_kind   ,&! rh2o/rair - 1.0
         vonkar    = 0.4_dbl_kind     ,&! von Karman constant
         cp_wv     = 1.81e3_dbl_kind  ,&! specific heat of water vapor (J/kg/K)
         stefan_boltzmann = 567.0e-10_dbl_kind,&!  W/m^2/K^4
         Tffresh   = 273.15_dbl_kind  ,&! freezing temp of fresh ice (K)
         Lsub      = 2.835e6_dbl_kind ,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = 2.501e6_dbl_kind ,&! latent heat, vaporization freshwater (J/kg)
         Lfresh    = Lsub-Lvap        ,&! latent heat of melting of fresh ice (J/kg)
         Timelt    = 0.0_dbl_kind     ,&! melting temperature, ice top surface  (C)
         Tsmelt    = 0.0_dbl_kind     ,&! melting temperature, snow top surface (C)
         ice_ref_salinity = 4._dbl_kind ,&! (ppt)
         ! ocn_ref_salinity = 34.7_dbl_kind,&! (ppt)
         iceruf   = 0.0005_dbl_kind   ,&! ice surface roughness (m)
         cprho    = cp_ocn*rhow       ,&! for ocean mixed layer (J kg / K m^3)

         ! for ice strength
         Cp       = 0.5_dbl_kind*gravit*(rhow-rhoi)*rhoi/rhow ,&! proport const for PE 
         Pstar    = 2.75e4_dbl_kind   ,&! constant in Hibler strength formula 
                                        ! (kstrength = 0) 
         Cstar    = 20._dbl_kind      ,&! constant in Hibler strength formula 
                                        ! (kstrength = 0) 

         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
         kappav = 1.4_dbl_kind ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
         !kappan = 17.6_dbl_kind,&! vis extnctn coef in ice, wvlngth<700nm (1/m)

         ! kice is not used for mushy thermo
         kice   = 2.03_dbl_kind  ,&! thermal conductivity of fresh ice(W/m/deg)
         ! kseaice is used only for zero-layer thermo
         kseaice= 2.00_dbl_kind  ,&! thermal conductivity of sea ice (W/m/deg)
                                   ! (used in zero layer thermodynamics option)
         ksno   = 0.30_dbl_kind  ,&! thermal conductivity of snow  (W/m/deg)
         zref   = 10._dbl_kind   ,&! reference height for stability (m)
         hs_min = 1.e-4_dbl_kind ,&! min snow thickness for computing zTsn (m)
         snowpatch = 0.02_dbl_kind, & ! parameter for fractional snow area (m)

         ! biogeochemistry
         R_gC2molC = 12.0107   , & ! molar mass of carbon
         sk_l = 0.03_dbl_kind      ! (m) skeletal layer thickness

      integer (kind=int_kind), parameter, public :: & 
         nspint = 3              ,& ! number of solar spectral intervals
         nspint_5bd = 5            ! number of solar spectral intervals used in SNICAR

      ! weights for albedos 
      ! 4 Jan 2007 BPB  Following are appropriate for complete cloud
      ! in a summer polar atmosphere with 1.5m bare sea ice surface:
      ! .636/.364 vis/nir with only 0.5% direct for each band.
      real (kind=dbl_kind), parameter, public :: &           ! currently used only
         awtvdr = 0.00318_dbl_kind, &! visible, direct  ! for history and
         awtidr = 0.00182_dbl_kind, &! near IR, direct  ! diagnostics
         awtvdf = 0.63282_dbl_kind, &! visible, diffuse
         awtidf = 0.36218_dbl_kind   ! near IR, diffuse

      real (kind=dbl_kind), parameter, public :: &
         qqqice  = 11637800._dbl_kind   ,&! for qsat over ice
         TTTice  = 5897.8_dbl_kind      ,&! for qsat over ice
         qqqocn  = 627572.4_dbl_kind    ,&! for qsat over ocn
         TTTocn  = 5107.4_dbl_kind        ! for qsat over ocn

      ! orbital parameters
      integer (kind=int_kind), public :: iyear_AD  ! Year to calculate orbit for
      real(kind=dbl_kind),public :: eccen  ! Earth's orbital eccentricity
      real(kind=dbl_kind),public :: obliqr ! Earth's obliquity in radians
      real(kind=dbl_kind),public :: lambm0 ! Mean longitude of perihelion at the
                                           ! vernal equinox (radians)
      real(kind=dbl_kind),public :: mvelpp ! Earth's moving vernal equinox longitude
                                           ! of perihelion + pi (radians)
      real(kind=dbl_kind),public :: obliq  ! obliquity in degrees
      real(kind=dbl_kind),public :: mvelp  ! moving vernal equinox long
      real(kind=dbl_kind),public :: decln  ! solar declination angle in radians
      real(kind=dbl_kind),public :: eccf   ! earth orbit eccentricity factor
      logical(kind=log_kind),public :: log_print ! Flags print of status/error
    
      ! snow parameters
      real (kind=dbl_kind), parameter, public :: &
         snwlvlfac =   0.3_dbl_kind, & ! 30% rule: fractional increase in snow depth
                                       ! over ridged ice, compared with level ice
         rhosmin   = 100.0_dbl_kind    ! minimum snow density (kg/m^3)

      !-----------------------------------------------------------------
      ! numbers used in column package
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
        c0   = 0.0_dbl_kind, &
        c1   = 1.0_dbl_kind, &
        c1p5 = 1.5_dbl_kind, &
        c2   = 2.0_dbl_kind, &
        c3   = 3.0_dbl_kind, &
        c4   = 4.0_dbl_kind, &
        c5   = 5.0_dbl_kind, &
        c6   = 6.0_dbl_kind, &
        c8   = 8.0_dbl_kind, &
        c10  = 10.0_dbl_kind, &
        c15  = 15.0_dbl_kind, &
        c16  = 16.0_dbl_kind, &
        c20  = 20.0_dbl_kind, &
        c25  = 25.0_dbl_kind, &
        c100 = 100.0_dbl_kind, &
        c1000= 1000.0_dbl_kind, &
        p001 = 0.001_dbl_kind, &
        p01  = 0.01_dbl_kind, &
        p1   = 0.1_dbl_kind, &
        p2   = 0.2_dbl_kind, &
        p4   = 0.4_dbl_kind, &
        p5   = 0.5_dbl_kind, &
        p6   = 0.6_dbl_kind, &
        p05  = 0.05_dbl_kind, &
        p15  = 0.15_dbl_kind, &
        p25  = 0.25_dbl_kind, &
        p75  = 0.75_dbl_kind, &
        p333 = c1/c3, &
        p666 = c2/c3, &
        puny   = 1.0e-11_dbl_kind, &
        bignum = 1.0e+30_dbl_kind, &
        pi     = 3.14159265358979323846_dbl_kind, &
        pih    = p5*pi

!=======================================================================

      end module ice_constants_colpkg

!=======================================================================
