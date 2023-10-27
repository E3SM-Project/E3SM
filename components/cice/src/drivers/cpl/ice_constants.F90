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
!  SVN:$Id: ice_constants.F90 37 2006-11-29 18:06:44Z eclare $
!
! author Elizabeth C. Hunke, LANL
!
! !INTERFACE:

      module ice_constants
!
! !USES:
!
      use shr_const_mod
      use ice_kinds_mod
!
!EOP
!
      implicit none
      save

      !-----------------------------------------------------------------
      ! physical constants
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter :: &
         pi        = SHR_CONST_PI    ,&! pi
         gravit    = SHR_CONST_G     ,&! gravitational acceleration (m/s^2)
         secday    = SHR_CONST_CDAY  ,&! seconds in calendar day
         omega     = SHR_CONST_OMEGA ,&! angular velocity of earth (rad/sec)
         radius    = SHR_CONST_REARTH,&! earth radius (m)
         rhos      = 330.0_dbl_kind  ,&! density of snow (kg/m^3)
         rhoi      = SHR_CONST_RHOICE,&! density of ice (kg/m^3)
         rhow      = SHR_CONST_RHOSW ,&! density of seawater (kg/m^3)
         rhofresh  = SHR_CONST_RHOFW ,&! density of fresh water (kg/m^3)
         zvir      = SHR_CONST_ZVIR  ,&! rh2o/rair - 1.0
         vonkar    = SHR_CONST_KARMAN,&! von Karman constant
         cp_air    = SHR_CONST_CPDAIR,&! specific heat of air (J/kg/K)
         cp_wv     = SHR_CONST_CPWV  ,&! specific heat of water vapor (J/kg/K)
         ! (Briegleb JGR 97 11475-11485  July 1992)
         emissivity = 0.95_dbl_kind  ,&! emissivity of snow and ice
         stefan_boltzmann = SHR_CONST_STEBOL,&!  W/m^2/K^4
         Tffresh   = SHR_CONST_TKFRZ ,&! freezing temp of fresh ice (K)
         cp_ice    = SHR_CONST_CPICE ,&! specific heat of fresh ice (J/kg/K)
         cp_ocn    = SHR_CONST_CPSW  ,&! specific heat of ocn    (J/kg/K)
         depressT  = 0.054_dbl_kind  ,&! Tf:brine salinity ratio (C/ppt)
         Lsub      = SHR_CONST_LATSUB,&! latent heat, sublimation freshwater (J/kg)
         Lvap      = SHR_CONST_LATVAP,&! latent heat, vaporization freshwater (J/kg)
         Lfresh    = SHR_CONST_LATICE,&! latent heat of melting of fresh ice (J/kg)
         Timelt    = SHR_CONST_TKFRZ-SHR_CONST_TKFRZ,&! melting temp. ice top surface  (C)
         Tsmelt    = SHR_CONST_TKFRZ-SHR_CONST_TKFRZ,&! melting temp. snow top surface (C)
         ice_ref_salinity = SHR_CONST_ICE_REF_SAL ,&! (psu)
!        ocn_ref_salinity = SHR_CONST_OCN_REF_SAL ,&! (psu)
         albocn = 0.06_dbl_kind                   ,&! ocean albedo
         dragio = 0.00536_dbl_kind                ,&! ice-ocn drag coefficient
!        rho_air   = SHR_CONST_RHODAIR,&! ambient air density (kg/m^3)
         spval_dbl = SHR_CONST_SPVAL  ,&! special value
         snowpatch = 0.005_dbl_kind     ! parameter for fractional snow area (m)

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
         Tocnfrz= -34.0_dbl_kind*depressT  ! freezing temp of seawater (C),
                                           ! used as Tsfcn for open water

      ! weights for albedos 
! 4 Jan 2007 BPB  Following are appropriate for complete cloud
! in a summer polar atmosphere with 1.5m bare sea ice surface:
! .636/.364 vis/nir with only 0.5% direct for each band.
      real (kind=dbl_kind), parameter :: &           ! currently used only
         awtvdr = 0.00318_dbl_kind, &! visible, direct  ! for history and
         awtidr = 0.00182_dbl_kind, &! near IR, direct
         awtvdf = 0.63282_dbl_kind, &! visible, diffuse  ! diagnostics
         awtidf = 0.36218_dbl_kind   ! near IR, diffuse

      real (kind=dbl_kind), parameter :: &
         hs0   = 0.03_dbl_kind,     &! parameter for delta-Eddington snow frac
         hsmin = 0.0001_dbl_kind     ! minimum snow thickness for dEdd

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
        c7   = 7.0_dbl_kind, &
        c8   = 8.0_dbl_kind, &
        c9   = 9.0_dbl_kind, &
        c10  = 10.0_dbl_kind, &
        c12  = 12.0_dbl_kind, &
        c15  = 15.0_dbl_kind, &
        c16  = 16.0_dbl_kind, &
        c20  = 20.0_dbl_kind, &
        c25  = 25.0_dbl_kind, &
        c90  = 90.0_dbl_kind, &
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
        eps04  = 1.0e-4_dbl_kind, &
        eps11  = 1.0e-11_dbl_kind, &
        eps12  = 1.0e-12_dbl_kind, &
        eps13  = 1.0e-13_dbl_kind, &
        eps15  = 1.0e-15_dbl_kind, &
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

#ifndef USE_ESMF
      integer (kind=int_kind), parameter :: &
         ESMF_SUCCESS = 0   ! otherwise ESMF defines this parameter
#endif

      ! useful for debugging
      integer (kind=int_kind), parameter :: &
         mtest = -999, itest = 1, jtest = 1, ntest = 1, btest = 1
!         mtest = 2, itest = 50, jtest = 53, ntest = 1, btest = 1

!=======================================================================

      end module ice_constants

!=======================================================================
