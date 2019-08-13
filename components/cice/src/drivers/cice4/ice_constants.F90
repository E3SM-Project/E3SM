!=======================================================================
!BOP
!
! !MODULE: ice_constants - sets physical constants
!
! !DESCRIPTION:
!
! This module defines a variety of physical and numerical constants
! used throughout the ice model \\
!
! Code originally based on constants.F in POP
!
! !REVISION HISTORY:
!  SVN:$Id: ice_constants.F90 143 2008-08-08 23:08:22Z eclare $
!
! author Elizabeth C. Hunke, LANL
!
! !INTERFACE:

      module ice_constants
!
! !USES:
!
      use ice_kinds_mod
!
!EOP
!
      implicit none
      save

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

#ifdef AOMIP
      real (kind=dbl_kind), parameter :: &
         rhos      = 300.0_dbl_kind   ,&! density of snow (kg/m^3)
         rhoi      = 900.0_dbl_kind   ,&! density of ice (kg/m^3)
         rhow      = 1025.0_dbl_kind  ,&! density of seawater (kg/m^3)
         cp_air    = 1004.0_dbl_kind  ,&! specific heat of air (J/kg/K)
         emissivity= 0.98_dbl_kind    ,&! emissivity of snow and ice
         cp_ice    = 2090._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = 4190._dbl_kind   ,&! specific heat of ocn    (J/kg/K)
         depressT  = 0.0575_dbl_kind  ,&! Tf:brine salinity ratio (C/ppt)
         dragio    = 0.0055_dbl_kind  ,&! ice-ocn drag coefficient
         albocn    = 0.10_dbl_kind      ! ocean albedo
#else
! CICE default parameters
      real (kind=dbl_kind), parameter :: &
         rhos      = 330.0_dbl_kind   ,&! density of snow (kg/m^3)
         rhoi      = 917.0_dbl_kind   ,&! density of ice (kg/m^3)
         rhow      = 1026.0_dbl_kind  ,&! density of seawater (kg/m^3)
         cp_air    = 1005.0_dbl_kind  ,&! specific heat of air (J/kg/K)
         ! (Briegleb JGR 97 11475-11485  July 1992)
         emissivity= 0.95_dbl_kind    ,&! emissivity of snow and ice
         cp_ice    = 2106._dbl_kind   ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = 4218._dbl_kind   ,&! specific heat of ocn    (J/kg/K)
                                        ! freshwater value needed for enthalpy
         depressT  = 0.054_dbl_kind   ,&! Tf:brine salinity ratio (C/ppt)
         dragio    = 0.00536_dbl_kind ,&! ice-ocn drag coefficient
         albocn    = 0.06_dbl_kind      ! ocean albedo
#endif

      real (kind=dbl_kind), parameter :: &
         gravit    = 9.80616_dbl_kind    ,&! gravitational acceleration (m/s^2)
         omega     = 7.292e-5_dbl_kind   ,&! angular velocity of earth (rad/sec)
         radius    = 6.37e6_dbl_kind       ! earth radius (m)

      real (kind=dbl_kind), parameter :: &
         pi = 3.14159265358979323846_dbl_kind,&! pi
         secday    = 86400.0_dbl_kind ,&! seconds in calendar day
         Tocnfrz   = -1.8_dbl_kind    ,&! freezing temp of seawater (C),
                                        ! used as Tsfcn for open water
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
         ice_ref_salinity = 4._dbl_kind ,&! (psu)
!        ocn_ref_salinity = 34.7_dbl_kind,&! (psu)
!        rho_air   = 1.2_dbl_kind     ,&! ambient air density (kg/m^3)
         spval_dbl = 1.0e30_dbl_kind    ! special value (double precision)

      real (kind=real_kind), parameter :: &
         spval     = 1.0e30_real_kind   ! special value for netCDF output

      real (kind=dbl_kind), parameter :: &
         iceruf   = 0.0005_dbl_kind   ,&! ice surface roughness (m)

         ! (Ebert, Schramm and Curry JGR 100 15965-15975 Aug 1995)
         kappav = 1.4_dbl_kind ,&! vis extnctn coef in ice, wvlngth<700nm (1/m)
         kappan = 17.6_dbl_kind,&! vis extnctn coef in ice, wvlngth<700nm (1/m)

         kice   = 2.03_dbl_kind  ,&! thermal conductivity of fresh ice(W/m/deg)
         kseaice= 2.00_dbl_kind  ,&! thermal conductivity of sea ice (W/m/deg)
                                   ! (used in zero layer thermodynamics option)
         ksno   = 0.30_dbl_kind  ,&! thermal conductivity of snow  (W/m/deg)
         zref   = 10._dbl_kind   ,&! reference height for stability (m)
         snowpatch = 0.02_dbl_kind ! parameter for fractional snow area (m)

      ! weights for albedos (match those for isccp shortwave forcing)
! 4 Jan 2007 BPB  Following are appropriate for complete cloud
! in a summer polar atmosphere with 1.5m bare sea ice surface:
! .636/.364 vis/nir with only 0.5% direct for each band.
      real (kind=dbl_kind), parameter :: &           ! currently used only
         awtvdr = 0.00318_dbl_kind, &! visible, direct  ! for history and
         awtidr = 0.63282_dbl_kind, &! near IR, direct  ! diagnostics
         awtvdf = 0.00182_dbl_kind, &! visible, diffuse
         awtidf = 0.36218_dbl_kind   ! near IR, diffuse

      real (kind=dbl_kind), parameter :: &
         qqqice  = 11637800._dbl_kind   ,&! for qsat over ice
         TTTice  = 5897.8_dbl_kind      ,&! for qsat over ice
         qqqocn  = 627572.4_dbl_kind    ,&! for qsat over ocn
         TTTocn  = 5107.4_dbl_kind        ! for qsat over ocn

      ! these are currently set so as to have no effect on the decomposition
      real (kind=dbl_kind), parameter :: &
         shlat  =  30.0_dbl_kind   ,&! artificial masking edge (deg)
         nhlat  = -30.0_dbl_kind     ! artificial masking edge (deg)
   
      !-----------------------------------------------------------------
      ! numbers
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter :: &
        c0   = 0.0_dbl_kind, &
        c1   = 1.0_dbl_kind, &
        c1p5 = 1.5_dbl_kind, &
        c2   = 2.0_dbl_kind, &
        c3   = 3.0_dbl_kind, &
        c4   = 4.0_dbl_kind, &
        c5   = 5.0_dbl_kind, &
        c6   = 6.0_dbl_kind, &
        c8   = 8.0_dbl_kind, &
        c9   = 9.0_dbl_kind, &
        c10  = 10.0_dbl_kind, &
        c12  = 12.0_dbl_kind, &
        c15  = 15.0_dbl_kind, &
        c16  = 16.0_dbl_kind, &
        c20  = 20.0_dbl_kind, &
        c25  = 25.0_dbl_kind, &
        c100 = 100.0_dbl_kind, &
        c180 = 180.0_dbl_kind, &
        c360 = 360.0_dbl_kind, &
        c365 = 365.0_dbl_kind, &
        c3600= 3600.0_dbl_kind, &
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
        p166 = c1/c6, &
        p333 = c1/c3, &
        p666 = c2/c3, &
        p111 = c1/c9, &
        p055 = p111*p5, &
        p027 = p055*p5, &
        p222 = c2/c9, &
        eps11  = 1.0e-11_dbl_kind, &
        eps13  = 1.0e-13_dbl_kind, &
        eps16  = 1.0e-16_dbl_kind, &
        puny   = eps11, &
        bignum = 1.0e+30_dbl_kind, &
        pih    = p5*pi, &
        pi2    = c2*pi

      !-----------------------------------------------------------------
      ! location of fields for staggered grids
      !-----------------------------------------------------------------

      integer (int_kind), parameter :: &   
        field_loc_unknown  =  0, & 
        field_loc_noupdate = -1, & 
        field_loc_center   =  1, & 
        field_loc_NEcorner =  2, & 
        field_loc_Nface    =  3, & 
        field_loc_Eface    =  4, &
        field_loc_Wface    =  5


      !-----------------------------------------------------------------
      ! field type attribute - necessary for handling
      ! changes of direction across tripole boundary
      !-----------------------------------------------------------------

      integer (int_kind), parameter :: &   
        field_type_unknown  =  0, & 
        field_type_noupdate = -1, & 
        field_type_scalar   =  1, & 
        field_type_vector   =  2, & 
        field_type_angle    =  3

      !-----------------------------------------------------------------
      ! conversion factors
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter :: &
        cm_to_m       = 0.01_dbl_kind   ,&! cm to meters
        m_to_cm       = 100._dbl_kind   ,&! meters to cm
        m2_to_km2     = 1.e-6_dbl_kind  ,&! m^2 to km^2
        kg_to_g       = 1000._dbl_kind  ,&! kilograms to grams
        mps_to_cmpdy  = 8.64e6_dbl_kind ,&! m per s to cm per day
        rad_to_deg    = 180._dbl_kind/pi  ! degree-radian conversion

!=======================================================================

      end module ice_constants

!=======================================================================
