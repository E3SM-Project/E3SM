!  SVN:$Id: ice_colpkg_shared.F90 1012 2015-06-26 12:34:09Z eclare $
!=========================================================================
!
! flags for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module ice_colpkg_shared

      use ice_kinds_mod
      use ice_constants_colpkg, only: c3

      implicit none

      private

!-----------------------------------------------------------------------
! Parameters for thermodynamics
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: &
         ktherm          ! type of thermodynamics
                         ! 0 = 0-layer approximation
                         ! 1 = Bitz and Lipscomb 1999
                         ! 2 = mushy layer theory

      character (char_len), public :: &
         conduct, &      ! 'MU71' or 'bubbly'
         fbot_xfer_type  ! transfer coefficient type for ice-ocean heat flux

      logical (kind=log_kind), public :: &
         heat_capacity, &! if true, ice has nonzero heat capacity
                         ! if false, use zero-layer thermodynamics
         calc_Tsfc       ! if true, calculate surface temperature
                         ! if false, Tsfc is computed elsewhere and
                         ! atmos-ice fluxes are provided to CICE

      real (kind=dbl_kind), parameter, public :: &
         saltmax = 3.2_dbl_kind,   & ! max salinity at ice base for BL99 (ppt)
         ! phi_init and dSin0_frazil are used for mushy thermo, ktherm=2
         phi_init = 0.75_dbl_kind, & ! initial liquid fraction of frazil
         dSin0_frazil = c3 ! bulk salinity reduction of newly formed frazil

      real (kind=dbl_kind), public :: &
         ustar_min       ! minimum friction velocity for ice-ocean heat flux

      ! mushy thermo
      real(kind=dbl_kind), public :: &
         a_rapid_mode      , & ! channel radius for rapid drainage mode (m)
         Rac_rapid_mode    , & ! critical Rayleigh number for rapid drainage mode
         aspect_rapid_mode , & ! aspect ratio for rapid drainage mode (larger=wider)
         dSdt_slow_mode    , & ! slow mode drainage strength (m s-1 K-1)
         phi_c_slow_mode   , & ! liquid fraction porosity cutoff for slow mode
         phi_i_mushy           ! liquid fraction of congelation ice

!-----------------------------------------------------------------------
! Parameters for radiation
!-----------------------------------------------------------------------

      character (len=char_len), public :: &
         shortwave, & ! shortwave method, 'default' ('ccsm3') or 'dEdd'
         albedo_type  ! albedo parameterization, 'default' ('ccsm3') or 'constant'
                      ! shortwave='dEdd' overrides this parameter

      ! baseline albedos for ccsm3 shortwave, set in namelist
      real (kind=dbl_kind), public :: &
         albicev  , & ! visible ice albedo for h > ahmax
         albicei  , & ! near-ir ice albedo for h > ahmax
         albsnowv , & ! cold snow albedo, visible
         albsnowi , & ! cold snow albedo, near IR
         ahmax        ! thickness above which ice albedo is constant (m)

      ! dEdd tuning parameters, set in namelist
      real (kind=dbl_kind), public :: &
         R_ice    , & ! sea ice tuning parameter; +1 > 1sig increase in albedo
         R_pnd    , & ! ponded ice tuning parameter; +1 > 1sig increase in albedo
         R_snw    , & ! snow tuning parameter; +1 > ~.01 change in broadband albedo
         dT_mlt   , & ! change in temp for non-melt to melt snow grain 
                      ! radius change (C)
         rsnw_mlt , & ! maximum melting snow grain radius (10^-6 m)
         kalg         ! algae absorption coefficient for 0.5 m thick layer

!-----------------------------------------------------------------------
! Parameters for ridging and strength
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: & ! defined in namelist 
         kstrength  , & ! 0 for simple Hibler (1979) formulation 
                        ! 1 for Rothrock (1975) pressure formulation 
         krdg_partic, & ! 0 for Thorndike et al. (1975) formulation 
                        ! 1 for exponential participation function 
         krdg_redist    ! 0 for Hibler (1980) formulation 
                        ! 1 for exponential redistribution function 

      real (kind=dbl_kind), public :: &  
         mu_rdg, &      ! gives e-folding scale of ridged ice (m^.5) 
                        ! (krdg_redist = 1) 
         Cf             ! ratio of ridging work to PE change in ridging (kstrength = 1)

!-----------------------------------------------------------------------
! Parameters for atmosphere
!-----------------------------------------------------------------------

      character (len=char_len), public :: &
         atmbndy ! atmo boundary method, 'default' ('ccsm3') or 'constant'

      logical (kind=log_kind), public :: &
         calc_strair, &  ! if true, calculate wind stress components
         formdrag,    &  ! if true, calculate form drag
         highfreq        ! if true, use high frequency coupling

      integer (kind=int_kind), public :: &
         natmiter        ! number of iterations for boundary layer calculations

!-----------------------------------------------------------------------
! Parameters for ocean
!-----------------------------------------------------------------------

      logical (kind=log_kind), public :: &
         oceanmixed_ice           ! if true, use ocean mixed layer

      character(len=char_len), public :: &
         tfrz_option              ! form of ocean freezing temperature
                                  ! 'minus1p8' = -1.8 C
                                  ! 'linear_salt' = -depressT * sss
                                  ! 'mushy' conforms with ktherm=2

!-----------------------------------------------------------------------
! Parameters for the ice thickness distribution
!-----------------------------------------------------------------------

      integer (kind=int_kind), public :: &
         kitd        , & ! type of itd conversions
                         !   0 = delta function
                         !   1 = linear remap
         kcatbound       !   0 = old category boundary formula
                         !   1 = new formula giving round numbers
                         !   2 = WMO standard
                         !   3 = asymptotic formula

!-----------------------------------------------------------------------
! Parameters for melt ponds
!-----------------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         hs0             ! snow depth for transition to bare sea ice (m)

      ! level-ice ponds
      character (len=char_len), public :: &
         frzpnd          ! pond refreezing parameterization

      real (kind=dbl_kind), public :: &
         dpscale, &      ! alter e-folding time scale for flushing 
         rfracmin, &     ! minimum retained fraction of meltwater
         rfracmax, &     ! maximum retained fraction of meltwater
         pndaspect, &    ! ratio of pond depth to pond fraction
         hs1             ! tapering parameter for snow on pond ice

      ! topo ponds
      real (kind=dbl_kind), public :: &
         hp1             ! critical parameter for pond ice thickness

!=======================================================================

      end module ice_colpkg_shared

!=======================================================================
