!  SVN:$Id: ice_colpkg_shared.F90 1177 2017-03-08 18:17:21Z eclare $
!=========================================================================
!
! flags for the column package
!
! authors: Elizabeth C. Hunke, LANL

      module ice_colpkg_shared

      use ice_kinds_mod
      use ice_constants_colpkg, only: c3, c0, c1, p5, p1

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
         calc_Tsfc   ,  &! if true, calculate surface temperature
                         ! if false, Tsfc is computed elsewhere and
                         ! atmos-ice fluxes are provided to CICE
         solve_zsal  ,  &! if true, update salinity profile from solve_S_dt
         modal_aero      ! if true, use modal aerosal optical properties
                         ! only for use with tr_aero or tr_zaero

      real (kind=dbl_kind), parameter, public :: &
         saltmax = 3.2_dbl_kind,   & ! max salinity at ice base for BL99 (ppt)
         ! phi_init and dSin0_frazil are used for mushy thermo, ktherm=2
         phi_init = 0.75_dbl_kind, & ! initial liquid fraction of frazil
         min_salin = p1          , & ! threshold for brine pocket treatment 
         salt_loss =0.4_dbl_kind, &  ! fraction of salt retained in zsalinity 
         min_bgc        = 0.01_dbl_kind, & ! fraction of ocean bgc concentration in surface melt 
         dSin0_frazil = c3 ! bulk salinity reduction of newly formed frazil

      real (kind=dbl_kind), public :: &
         dts_b,   &      ! zsalinity timestep
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

      real (kind=dbl_kind), parameter, public :: &
         hi_ssl = 0.050_dbl_kind, & ! ice surface scattering layer thickness (m)
         hs_ssl = 0.040_dbl_kind    ! snow surface scattering layer thickness (m)

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

!-----------------------------------------------------------------------
! Parameters for biogeochemistry
!-----------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! dimensions
      !-----------------------------------------------------------------
      integer (kind=int_kind), parameter, public :: &
         max_algae  =   3       , & ! maximum number of algal types 
         max_dic    =   1       , & ! maximum number of dissolved inorganic carbon types 
         max_doc    =   3       , & ! maximum number of dissolved organic carbon types
         max_don    =   1       , & ! maximum number of dissolved organic nitrogen types
         max_fe     =   2       , & ! maximum number of iron types
         nmodal1    =   10      , & ! dimension for modal aerosol radiation parameters
         nmodal2    =   8       , & ! dimension for modal aerosol radiation parameters
         max_aero   =   6       , & ! maximum number of aerosols 
         max_nbtrcr = max_algae*2 & ! algal nitrogen and chlorophyll
                    + max_dic     & ! dissolved inorganic carbon
                    + max_doc     & ! dissolved organic carbon
                    + max_don     & ! dissolved organic nitrogen
                    + 5           & ! nitrate, ammonium, silicate, PON, and humics
                    + 3           & ! DMSPp, DMSPd, DMS
                    + max_fe*2    & ! dissolved Fe and  particulate Fe
                    + max_aero      ! aerosols

      !-----------------------------------------------------------------
      ! namelist
      !-----------------------------------------------------------------
      character(char_len_long), public :: & 
         bgc_data_dir   ! directory for biogeochemistry data

      character(char_len), public :: &          
         sil_data_type  , & ! 'default', 'clim'
         nit_data_type  , & ! 'default', 'clim'   
         fe_data_type   , & ! 'default', 'clim'      
         bgc_flux_type      ! type of ocean-ice piston velocity 
                            ! 'constant', 'Jin2006' 

      logical (kind=log_kind), public :: &
         z_tracers,      & ! if .true., bgc or aerosol tracers are vertically resolved
         scale_bgc,      & ! if .true., initialize bgc tracers proportionally with salinity
         solve_zbgc,     & ! if .true., solve vertical biochemistry portion of code
         dEdd_algae        ! if .true., algal absorption of Shortwave is computed in the
        
      logical (kind=log_kind), public :: & 
         skl_bgc            ! if true, solve skeletal biochemistry

      real (kind=dbl_kind), public :: & 
         grid_o      , & ! for bottom flux        
         l_sk        , & ! characteristic diffusive scale (zsalinity) (m)
         grid_o_t    , & ! top grid point length scale 
         phi_snow    , & ! porosity of snow
         initbio_frac, & ! fraction of ocean tracer concentration used to initialize tracer 
         frazil_scav     ! multiple of ocean tracer concentration due to frazil scavenging

      real (kind=dbl_kind), public :: & 
         grid_oS     , & ! for bottom flux (zsalinity)
         l_skS           ! 0.02 characteristic skeletal layer thickness (m) (zsalinity)

      logical (kind=log_kind), public :: & 
         restore_bgc      ! if true, restore nitrate

      !-----------------------------------------------------------------
      ! From ice_zbgc_shared.F90
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         ratio_Si2N_diatoms, &   ! algal Si to N (mol/mol)
         ratio_Si2N_sp     , &
         ratio_Si2N_phaeo  , &
         ratio_S2N_diatoms , &   ! algal S  to N (mol/mol)
         ratio_S2N_sp      , &
         ratio_S2N_phaeo   , &
         ratio_Fe2C_diatoms, &   ! algal Fe to C  (umol/mol)
         ratio_Fe2C_sp     , &
         ratio_Fe2C_phaeo  , &
         ratio_Fe2N_diatoms, &   ! algal Fe to N  (umol/mol)
         ratio_Fe2N_sp     , &
         ratio_Fe2N_phaeo  , &
         ratio_Fe2DON      , &   ! Fe to N of DON (nmol/umol)
         ratio_Fe2DOC_s    , &   ! Fe to C of DOC (nmol/umol) saccharids
         ratio_Fe2DOC_l    , &   ! Fe to C of DOC (nmol/umol) lipids
         fr_resp           , &   ! fraction of algal growth lost due to respiration        
         tau_min           , &   ! rapid mobile to stationary exchanges (s) = 1.5 hours
         tau_max           , &   ! long time mobile to stationary exchanges (s) = 2 days
         algal_vel         , &   ! 0.5 cm/d(m/s) Lavoie 2005  1.5 cm/day
         R_dFe2dust        , &   !  g/g (3.5% content) Tagliabue 2009
         dustFe_sol              ! solubility fraction

      !-----------------------------------------------------------------
      ! From algal_dyn in ice_algae.F90
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         chlabs_diatoms   , & ! chl absorption (1/m/(mg/m^3))
         chlabs_sp        , & !
         chlabs_phaeo     , & !
         alpha2max_low_diatoms , & ! light limitation (1/(W/m^2))  
         alpha2max_low_sp      , & 
         alpha2max_low_phaeo   , & 
         beta2max_diatoms , & ! light inhibition (1/(W/m^2))  
         beta2max_sp      , & 
         beta2max_phaeo   , & 
         mu_max_diatoms   , & ! maximum growth rate (1/day)       
         mu_max_sp        , & 
         mu_max_phaeo     , & 
         grow_Tdep_diatoms, & ! Temperature dependence of growth (1/C)
         grow_Tdep_sp     , & 
         grow_Tdep_phaeo  , & 
         fr_graze_diatoms , & ! Fraction grazed
         fr_graze_sp      , & 
         fr_graze_phaeo   , & 
         mort_pre_diatoms , & ! Mortality (1/day)
         mort_pre_sp      , & 
         mort_pre_phaeo   , & 
         mort_Tdep_diatoms, & ! T dependence of mortality (1/C)
         mort_Tdep_sp     , &  
         mort_Tdep_phaeo  , &  
         k_exude_diatoms  , & ! algal exudation (1/d)
         k_exude_sp       , &  
         k_exude_phaeo    , &  
         K_Nit_diatoms    , & ! nitrate half saturation (mmol/m^3)
         K_Nit_sp         , &  
         K_Nit_phaeo      , &  
         K_Am_diatoms     , & ! ammonium half saturation (mmol/m^3)
         K_Am_sp          , &   
         K_Am_phaeo       , &   
         K_Sil_diatoms    , & ! silicate half saturation (mmol/m^3)
         K_Sil_sp         , &   
         K_Sil_phaeo      , &   
         K_Fe_diatoms     , & ! iron half saturation (nM)
         K_Fe_sp          , &   
         K_Fe_phaeo       , &    
         f_don_protein    , & ! fraction of spilled grazing to proteins           
         kn_bac_protein   , & ! Bacterial degredation of DON (1/d)                
         f_don_Am_protein , & ! fraction of remineralized DON to ammonium         
         f_doc_s          , & ! fraction of mortality to DOC 
         f_doc_l          , &   
         f_exude_s        , & ! fraction of exudation to DOC
         f_exude_l        , & 
         k_bac_s          , & ! Bacterial degredation of DOC (1/d)
         k_bac_l          , & 
         T_max            , & ! maximum temperature (C)
         fsal             , & ! Salinity limitation (ppt)
         op_dep_min       , & ! Light attenuates for optical depths exceeding min
         fr_graze_s       , & ! fraction of grazing spilled or slopped
         fr_graze_e       , & ! fraction of assimilation excreted 
         fr_mort2min      , & ! fractionation of mortality to Am
         fr_dFe           , & ! fraction of remineralized nitrogen (in units of algal iron)
         k_nitrif         , & ! nitrification rate (1/day)           
         t_iron_conv      , & ! desorption loss pFe to dFe (day)
         max_loss         , & ! restrict uptake to % of remaining value 
         max_dfe_doc1     , & ! max ratio of dFe to saccharides in the ice (nM Fe/muM C)    
         fr_resp_s        , & ! DMSPd fraction of respiration loss as DMSPd
         y_sk_DMS         , & ! fraction conversion given high yield
         t_sk_conv        , & ! Stefels conversion time (d)
         t_sk_ox              ! DMS oxidation time (d)

      !-----------------------------------------------------------------
      ! former parameters now in namelist
      !-----------------------------------------------------------------

      real (kind=dbl_kind), public :: &
         algaltype_diatoms  , & ! mobility type
         algaltype_sp       , & !
         algaltype_phaeo    , & !
         nitratetype        , & !
         ammoniumtype       , & !
         silicatetype       , & !
         dmspptype          , & !
         dmspdtype          , & !
         humtype            , & !
         doctype_s          , & !
         doctype_l          , & !
         dontype_protein    , & !
         fedtype_1          , & !
         feptype_1          , & !
         zaerotype_bc1      , & !
         zaerotype_bc2      , & !
         zaerotype_dust1    , & !
         zaerotype_dust2    , & !
         zaerotype_dust3    , & !
         zaerotype_dust4    , & !
         ratio_C2N_diatoms  , & ! algal C to N ratio (mol/mol)
         ratio_C2N_sp       , & !
         ratio_C2N_phaeo    , & !
         ratio_chl2N_diatoms, & ! algal chlorophyll to N ratio (mg/mmol)
         ratio_chl2N_sp     , & !
         ratio_chl2N_phaeo  , & !
         F_abs_chl_diatoms  , & ! scales absorbed radiation for dEdd
         F_abs_chl_sp       , & !
         F_abs_chl_phaeo    , & !
         ratio_C2N_proteins     ! ratio of C to N in proteins (mol/mol)       

      !-----------------------------------------------------------------
      ! Transport type 
      !-----------------------------------------------------------------
      ! In delta Eddington, algal particles are assumed to cause no
      ! significant scattering (Brieglib and Light), only absorption
      ! in the visible spectral band (200-700 nm)
      ! Algal types: Diatoms, flagellates, Phaeocycstis
      ! DOC        : Proteins, EPS, Lipids
      !-----------------------------------------------------------------
      real (kind=dbl_kind), parameter, dimension(max_dic), public :: &
         dictype   = (/-c1/)  ! not in namelist

      real (kind=dbl_kind), dimension(max_algae), public :: &
         algaltype   ! tau_min for both retention and release

      real (kind=dbl_kind), dimension(max_doc), public :: &
         doctype 

      real (kind=dbl_kind), dimension(max_don), public :: &
         dontype  

      real (kind=dbl_kind), dimension(max_fe), public :: &
         fedtype 

      real (kind=dbl_kind), dimension(max_fe), public :: &
         feptype  

      !------------------------------------------------------------
      ! Aerosol order and type should be consistent with order/type 
      ! specified in delta Eddington:  1) hydrophobic black carbon;
      ! 2) hydrophilic black carbon; 3) dust (0.05-0.5 micron);
      ! 4) dust (0.5-1.25 micron); 5) dust (1.25-2.5 micron);
      ! 6) dust (2.5-5 micron) 
      !-------------------------------------------------------------
      real (kind=dbl_kind), dimension(max_aero), public :: &
         zaerotype  

      !-----------------------------------------------------------------
      ! Forcing input, history and diagnostic output
      !-----------------------------------------------------------------

      real (kind=dbl_kind), parameter, public :: &
         rhosi     = 940.0_dbl_kind, & ! average sea ice density
                                       ! Cox and Weeks, 1982: 919-974 kg/m^2
         sk_l      = 0.03_dbl_kind     ! skeletal layer thickness (m)

      real (kind=dbl_kind), dimension(max_algae), public :: &
         R_C2N     ,      & ! algal C to N (mole/mole) 
         R_chl2N   ,      & ! 3 algal chlorophyll to N (mg/mmol)
         F_abs_chl          ! to scale absorption in Dedd

      real (kind=dbl_kind), dimension(max_don), public :: &  ! increase compare to algal R_Fe2C
         R_C2N_DON 

!=======================================================================

      end module ice_colpkg_shared

!=======================================================================
