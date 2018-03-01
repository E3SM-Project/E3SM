module EDTypesMod

  use FatesConstantsMod,     only : r8 => fates_r8
  use FatesGlobals,          only : fates_log
  use shr_infnan_mod,        only : nan => shr_infnan_nan, assignment(=)

  use FatesHydraulicsMemMod, only : ed_cohort_hydr_type
  use FatesHydraulicsMemMod, only : ed_patch_hydr_type
  use FatesHydraulicsMemMod, only : ed_site_hydr_type

  implicit none
  save

  integer, parameter :: maxPatchesPerSite  = 10   ! maximum number of patches to live on a site
  integer, parameter :: maxCohortsPerPatch = 160  ! maximum number of cohorts per patch

  integer, parameter :: nclmax = 2                ! Maximum number of canopy layers
  integer, parameter :: ican_upper = 1            ! Nominal index for the upper canopy
  integer, parameter :: ican_ustory = 2           ! Nominal index for understory in two-canopy system

  integer, parameter :: nlevleaf = 40             ! number of leaf layers in canopy layer
  integer, parameter :: maxpft = 10               ! maximum number of PFTs allowed
                                                  ! the parameter file may determine that fewer
                                                  ! are used, but this helps allocate scratch
                                                  ! space and output arrays.
 

  ! TODO: we use this cp_maxSWb only because we have a static array q(size=2) of
  ! land-ice abledo for vis and nir.  This should be a parameter, which would
  ! get us on track to start using multi-spectral or hyper-spectral (RGK 02-2017)
  integer, parameter :: maxSWb = 2      ! maximum number of broad-bands in the
                                        ! shortwave spectrum cp_numSWb <= cp_maxSWb
                                        ! this is just for scratch-array purposes
                                        ! if cp_numSWb is larger than this value
                                        ! simply bump this number up as needed

  integer, parameter :: ivis = 1        ! This is the array index for short-wave
                                        ! radiation in the visible spectrum, as expected
                                        ! in boundary condition files and parameter
                                        ! files.  This will be compared with 
                                        ! the HLM's expectation in FatesInterfaceMod
  integer, parameter :: inir = 2        ! This is the array index for short-wave
                                        ! radiation in the near-infrared spectrum, as expected
                                        ! in boundary condition files and parameter
                                        ! files.  This will be compared with 
                                        ! the HLM's expectation in FatesInterfaceMod

  ! Switches that turn on/off ED dynamics process (names are self explanatory)
  ! IMPORTANT NOTE!!! THESE SWITCHES ARE EXPERIMENTAL.  
  ! THEY SHOULD CORRECTLY TURN OFF OR ON THE PROCESS, BUT.. THERE ARE VARIOUS 
  ! ASPECTS REGARDING DIAGNOSING RATES AND HOW THEY ARE REPORTED WHEN THESE 
  ! PROCESSES ARE OFF THAT NEED TO BE DISCUSSED AND CONSIDERED.
  ! TO-DO: THESE SHOULD BE PARAMETERS IN THE FILE OR NAMELIST - ADDING THESE
  ! WAS OUTSIDE THE SCOPE OF THE VERY LARGE CHANGESET WHERE THESE WERE FIRST
  ! INTRODUCED (RGK 03-2017)
  logical, parameter :: do_ed_phenology = .true.


  ! MODEL PARAMETERS
  real(r8), parameter :: AREA                 = 10000.0_r8 ! Notional area of simulated forest m2
  real(r8), parameter :: AREA_INV             = 1.0e-4_r8  ! Inverse of the notion area (faster math)

  integer, parameter :: numWaterMem           = 10         ! watermemory saved as site level var

  ! BIOLOGY/BIOGEOCHEMISTRY        
  integer , parameter :: external_recruitment = 0          ! external recruitment flag 1=yes  
  integer , parameter :: SENES                = 10         ! Window of time over which we track temp for cold sensecence (days)
  real(r8), parameter :: DINC_ED              = 1.0_r8     ! size of LAI bins. 
  integer , parameter :: N_DIST_TYPES         = 3          ! Disturbance Modes 1) tree-fall, 2) fire, 3) logging
  integer , parameter :: dtype_ifall          = 1          ! index for naturally occuring tree-fall generated event
  integer , parameter :: dtype_ifire          = 2          ! index for fire generated disturbance event
  integer , parameter :: dtype_ilog           = 3          ! index for logging generated disturbance event

  ! SPITFIRE     
  integer,  parameter :: NCWD                 = 4          ! number of coarse woody debris pools (twig,s branch,l branch, trunk)
  integer , parameter :: NFSC                 = NCWD+2     ! number fuel size classes  (4 cwd size classes, leaf litter, and grass)
  integer,  parameter :: lg_sf                = 6          ! array index of live grass pool for spitfire
  integer,  parameter :: dl_sf                = 1          ! array index of dead leaf pool for spitfire (dead grass and dead leaves)
  integer,  parameter :: tw_sf                = 2          ! array index of twig pool for spitfire
  integer,  parameter :: tr_sf                = 5          ! array index of dead trunk pool for spitfire
  integer,  parameter :: lb_sf                = 4          ! array index of large branch pool for spitfire 
  real(r8), parameter :: fire_threshold       = 50.0_r8    ! threshold for fires that spread or go out. KWm-2 (Pyne 1986)

  ! PATCH FUSION 
  real(r8), parameter :: NTOL                 = 0.05_r8    ! min plant density for hgt bin to be used in height profile comparisons 
  real(r8), parameter :: HITEMAX              = 30.0_r8    ! max dbh value used in hgt profile comparison 
  real(r8), parameter :: DBHMAX               = 150.0_r8   ! max dbh value used in hgt profile comparison 
  integer , parameter :: N_HITE_BINS          = 60         ! no. of hite bins used to distribute LAI
  integer , parameter :: N_DBH_BINS           = 5          ! no. of dbh bins used when comparing patches


  real(r8), parameter :: min_npm2       = 1.0E-8_r8  ! minimum cohort number density per m2 before termination
  real(r8), parameter :: min_patch_area = 0.001_r8   ! smallest allowable patch area before termination
  real(r8), parameter :: min_nppatch    = 1.0E-11_r8 ! minimum number of cohorts per patch (min_npm2*min_patch_area)
  real(r8), parameter :: min_n_safemath = 1.0E-15_r8 ! in some cases, we want to immediately remove super small
                                                     ! number densities of cohorts to prevent FPEs

  character*4 yearchar                    

  ! special mode to cause PFTs to create seed mass of all currently-existing PFTs
  logical, parameter :: homogenize_seed_pfts  = .false.

  integer, parameter :: nlevmclass_ed = 5      ! nlev "mortality" classes in ED
                                               ! Number of ways to die
                                               ! (background,hydraulic,carbon,impact,fire)

  character(len = 10), parameter,dimension(nlevmclass_ed) :: char_list = &
       (/"background","hydraulic ","carbon    ","impact    ","fire      "/)


  !************************************
  !** COHORT type structure          **
  !************************************
  type ed_cohort_type

     ! POINTERS
     type (ed_cohort_type) , pointer :: taller   => null()       ! pointer to next tallest cohort     
     type (ed_cohort_type) , pointer :: shorter  => null()       ! pointer to next shorter cohort     
     type (ed_patch_type)  , pointer :: patchptr => null()       ! pointer to patch that cohort is in
     type (ed_site_type)   , pointer :: siteptr  => null()       ! pointer to site that cohort is in

     ! VEGETATION STRUCTURE
     integer  ::  pft                                    ! pft number
     real(r8) ::  n                                      ! number of individuals in cohort per 'area' (10000m2 default)
     real(r8) ::  dbh                                    ! dbh: cm
     real(r8) ::  hite                                   ! height: meters
     integer  ::  indexnumber                            ! unique number for each cohort. (within clump?)
     real(r8) ::  balive                                 ! total living biomass: kGC per indiv
     real(r8) ::  bdead                                  ! dead biomass:  kGC per indiv
     real(r8) ::  bstore                                 ! stored carbon: kGC per indiv
     real(r8) ::  laimemory                              ! target leaf biomass- set from previous year: kGC per indiv
     integer  ::  canopy_layer                           ! canopy status of cohort (1 = canopy, 2 = understorey, etc.)
     real(r8) ::  canopy_layer_yesterday                 ! recent canopy status of cohort (1 = canopy, 2 = understorey, etc.)  real to be conservative during fusion
     real(r8) ::  b                                      ! total biomass: kGC per indiv
     real(r8) ::  bsw                                    ! sapwood in stem and roots: kGC per indiv
     real(r8) ::  bl                                     ! leaf biomass: kGC per indiv
     real(r8) ::  br                                     ! fine root biomass: kGC per indiv
     real(r8) ::  lai                                    ! leaf area index of cohort   m2/m2
     real(r8) ::  sai                                    ! stem area index of cohort   m2/m2
     real(r8) ::  gscan                                  ! Stomatal resistance of cohort. 
     real(r8) ::  canopy_trim                            ! What is the fraction of the maximum leaf biomass that we are targeting? :-
     real(r8) ::  leaf_cost                              ! How much does it cost to maintain leaves: kgC/m2/year-1
     real(r8) ::  excl_weight                            ! How much of this cohort is demoted each year, as a proportion of all cohorts:-
     real(r8) ::  prom_weight                            ! How much of this cohort is promoted each year, as a proportion of all cohorts:-
     integer  ::  nv                                     ! Number of leaf layers: -
     integer  ::  status_coh                             ! growth status of plant  (2 = leaves on , 1 = leaves off)
     real(r8) ::  c_area                                 ! areal extent of canopy (m2)
     real(r8) ::  treelai                                ! lai of tree (total leaf area (m2) / canopy area (m2)
     real(r8) ::  treesai                                ! stem area index of tree (total stem area (m2) / canopy area (m2)
     logical  ::  isnew                                  ! flag to signify a new cohort, new cohorts have not experienced
                                                         ! npp or mortality and should therefore not be fused or averaged
     integer  ::  size_class                             ! An index that indicates which diameter size bin the cohort currently resides in
                                                         ! this is used for history output. We maintain this in the main cohort memory
                                                         ! because we don't want to continually re-calculate the cohort's position when
                                                         ! performing size diagnostics at high-frequency calls
     integer  ::  size_by_pft_class                      ! An index that indicates the cohorts position of the joint size-class x functional
                                                         ! type classification. We also maintain this in the main cohort memory
                                                         ! because we don't want to continually re-calculate the cohort's position when
                                                         ! performing size diagnostics at high-frequency calls


     ! CARBON FLUXES 
     
     ! ----------------------------------------------------------------------------------
     ! NPP, GPP and RESP: Instantaneous, accumulated and accumulated-hold types.*
     ! 
     ! _tstep:    The instantaneous estimate that is calculated at each rapid plant biophysics
     !            time-step (ie photosynthesis, sub-hourly). (kgC/indiv/timestep)
     ! _acc:      The accumulation of the _tstep variable from the beginning to ending of
     !            the dynamics time-scale.  This variable is zero'd during initialization and
     !            after the dynamics call-sequence is completed.  (kgC/indiv/day)
     ! _acc_hold: While _acc is zero'd after the dynamics call sequence and then integrated, 
     !            _acc_hold "holds" the integrated value until the next time dynamics is 
     !            called. This is necessary for restarts. This variable also has units
     !            converted to a useful rate (kgC/indiv/yr)
     ! ----------------------------------------------------------------------------------

     real(r8) ::  gpp_tstep          ! Gross Primary Production (see above *)
     real(r8) ::  gpp_acc
     real(r8) ::  gpp_acc_hold

     real(r8) ::  npp_tstep          ! Net Primary Production (see above *)
     real(r8) ::  npp_acc
     real(r8) ::  npp_acc_hold

     real(r8) ::  resp_tstep         ! Autotrophic respiration (see above *)
     real(r8) ::  resp_acc
     real(r8) ::  resp_acc_hold

     ! Net Primary Production Partitions

     real(r8) ::  npp_leaf                               ! NPP into leaves (includes replacement of turnover):  KgC/indiv/year
     real(r8) ::  npp_froot                              ! NPP into fine roots (includes replacement of turnover):  KgC/indiv/year

     real(r8) ::  npp_bsw                                ! NPP into sapwood: KgC/indiv/year
     real(r8) ::  npp_bdead                              ! NPP into deadwood (structure):  KgC/indiv/year
     real(r8) ::  npp_bseed                              ! NPP into seeds: KgC/indiv/year
     real(r8) ::  npp_store                              ! NPP into storage: KgC/indiv/year

     real(r8) ::  ts_net_uptake(nlevleaf)              ! Net uptake of leaf layers: kgC/m2/s
     real(r8) ::  year_net_uptake(nlevleaf)            ! Net uptake of leaf layers: kgC/m2/year

     ! RESPIRATION COMPONENTS
     real(r8) ::  rdark                                  ! Dark respiration: kgC/indiv/s
     real(r8) ::  resp_g                                 ! Growth respiration:  kgC/indiv/timestep
     real(r8) ::  resp_m                                 ! Maintenance respiration:  kgC/indiv/timestep 
     real(r8) ::  livestem_mr                            ! Live stem        maintenance respiration: kgC/indiv/s
                                                         ! (Above ground)
     real(r8) ::  livecroot_mr                           ! Live stem        maintenance respiration: kgC/indiv/s
                                                         ! (below ground)
     real(r8) ::  froot_mr                               ! Live fine root   maintenance respiration: kgC/indiv/s

     ! ALLOCATION
     real(r8) ::  md                                     ! plant maintenance demand: kgC/indiv/year
     real(r8) ::  leaf_md                                ! leaf  maintenance demand: kgC/indiv/year
     real(r8) ::  root_md                                ! root  maintenance demand: kgC/indiv/year
     real(r8) ::  carbon_balance                         ! carbon remaining for growth and storage: kg/indiv/year
     real(r8) ::  seed_prod                              ! reproduction seed and clonal: KgC/indiv/year
     real(r8) ::  leaf_litter                            ! leaf litter from phenology: KgC/m2
     real(r8) ::  woody_turnover                         ! amount of wood lost each day: kgC/indiv/year. Currently set to zero.

     !MORTALITY
     real(r8) ::  dmort                                  ! proportional mortality rate. (year-1)

     ! Mortality Rate Partitions
     real(r8) ::  bmort                                  ! background mortality rate        n/year
     real(r8) ::  cmort                                  ! carbon starvation mortality rate n/year
     real(r8) ::  hmort                                  ! hydraulic failure mortality rate n/year
     real(r8) ::  imort                                  ! mortality from impacts by others n/year
     real(r8) ::  fmort                                  ! fire mortality                   n/year

      ! Logging Mortality Rate 
	 ! Yi Xu
     real(r8) ::  lmort_logging                          ! directly logging rate            %/per logging activity
     real(r8) ::  lmort_collateral                       ! collaterally damaged rate        %/per logging activity
     real(r8) ::  lmort_infra                            ! mechanically damaged rate        %/per logging activity
	      

     ! NITROGEN POOLS      
     ! ----------------------------------------------------------------------------------
     ! Nitrogen pools are not prognostic in the current implementation.
     ! They are diagnosed during photosynthesis using a simple C2N parameter. Local values
     ! used in that routine.
     ! ----------------------------------------------------------------------------------

     ! GROWTH DERIVIATIVES
     real(r8) ::  dndt                                   ! time derivative of cohort size  : n/year
     real(r8) ::  dhdt                                   ! time derivative of height       : m/year
     real(r8) ::  ddbhdt                                 ! time derivative of dbh          : cm/year
     real(r8) ::  dbalivedt                              ! time derivative of total living biomass : KgC/year
     real(r8) ::  dbdeaddt                               ! time derivative of dead biomass         : KgC/year
     real(r8) ::  dbstoredt                              ! time derivative of stored biomass       : KgC/year
     real(r8) ::  storage_flux                           ! flux from npp into bstore               : KgC/year

     ! FIRE
     real(r8) ::  cfa                                    ! proportion of crown affected by fire:-
     real(r8) ::  cambial_mort                           ! probability that trees dies due to cambial char:-
     real(r8) ::  crownfire_mort                         ! probability of tree post-fire mortality due to crown scorch:-
     real(r8) ::  fire_mort                              ! post-fire mortality from cambial and crown damage assuming two are independent:-

     ! Hydraulics
     type(ed_cohort_hydr_type), pointer :: co_hydr       ! All cohort hydraulics data, see FatesHydraulicsMemMod.F90


  end type ed_cohort_type

  !************************************
  !** Patch type structure           **
  !************************************

  type ed_patch_type

     ! POINTERS
     type (ed_cohort_type), pointer :: tallest => null()           ! pointer to patch's tallest cohort    
     type (ed_cohort_type), pointer :: shortest => null()          ! pointer to patch's shortest cohort
     type (ed_patch_type),  pointer :: older => null()             ! pointer to next older patch   
     type (ed_patch_type),  pointer :: younger => null()           ! pointer to next younger patch      
     type (ed_site_type),   pointer :: siteptr => null()           ! pointer to the site that the patch is in

     !INDICES
     integer  :: patchno                                           ! unique number given to each new patch created for tracking

     ! PATCH INFO
     real(r8) ::  age                                              ! average patch age: years                   
     integer  ::  age_class                                        ! age class of the patch for history binning purposes
     real(r8) ::  area                                             ! patch area: m2  
     integer  ::  countcohorts                                     ! Number of cohorts in patch
     integer  ::  ncl_p                                            ! Number of occupied canopy layers

     ! LEAF ORGANIZATION
     real(r8) ::  pft_agb_profile(maxpft,n_dbh_bins)            ! binned above ground biomass, for patch fusion: KgC/m2
     real(r8) ::  canopy_layer_lai(nclmax)                         ! lai that is shading this canopy layer: m2/m2 
     real(r8) ::  total_canopy_area                                ! area that is covered by vegetation : m2
     real(r8) ::  total_tree_area                                  ! area that is covered by woody vegetation : m2
     real(r8) ::  canopy_area                                      ! area that is covered by vegetation : m2 (is this different to total_canopy_area?
     real(r8) ::  bare_frac_area                                   ! bare soil in this patch expressed as a fraction of the total soil surface.
     real(r8) ::  lai                                              ! leaf area index of patch
     real(r8) ::  zstar                                            ! height of smallest canopy tree -- only meaningful in "strict PPA" mode

     real(r8) ::  tlai_profile(nclmax,maxpft,nlevleaf)        ! total   leaf area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  elai_profile(nclmax,maxpft,nlevleaf)        ! exposed leaf area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  tsai_profile(nclmax,maxpft,nlevleaf)        ! total   stem area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  esai_profile(nclmax,maxpft,nlevleaf)        ! exposed stem area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  layer_height_profile(nclmax,maxpft,nlevleaf)
     real(r8) ::  canopy_area_profile(nclmax,maxpft,nlevleaf) ! fraction of canopy in each canopy 
     ! layer, pft, and leaf layer:-
     integer  ::  present(nclmax,maxpft)                        ! is there any of this pft in this canopy layer?      
     integer  ::  nrad(nclmax,maxpft)                           ! number of exposed leaf layers for each canopy layer and pft
     integer  ::  ncan(nclmax,maxpft)                           ! number of total   leaf layers for each canopy layer and pft

     !RADIATION FLUXES      
     real(r8) ::  fabd_sun_z(nclmax,maxpft,nlevleaf)          ! sun fraction of direct light absorbed by each canopy 
     ! layer, pft, and leaf layer:-
     real(r8) ::  fabd_sha_z(nclmax,maxpft,nlevleaf)          ! shade fraction of direct light absorbed by each canopy 
     ! layer, pft, and leaf layer:-
     real(r8) ::  fabi_sun_z(nclmax,maxpft,nlevleaf)          ! sun fraction of indirect light absorbed by each canopy 
     ! layer, pft, and leaf layer:-
     real(r8) ::  fabi_sha_z(nclmax,maxpft,nlevleaf)          ! shade fraction of indirect light absorbed by each canopy 
     ! layer, pft, and leaf layer:-

     real(r8) ::  ed_laisun_z(nclmax,maxpft,nlevleaf)         ! amount of LAI in the sun   in each canopy layer, 
     ! pft, and leaf layer. m2/m2
     real(r8) ::  ed_laisha_z(nclmax,maxpft,nlevleaf)         ! amount of LAI in the shade in each canopy layer,
     real(r8) ::  ed_parsun_z(nclmax,maxpft,nlevleaf)         ! PAR absorbed  in the sun   in each canopy layer,
     real(r8) ::  ed_parsha_z(nclmax,maxpft,nlevleaf)         ! PAR absorbed  in the shade in each canopy layer,
     real(r8) ::  f_sun(nclmax,maxpft,nlevleaf)               ! fraction of leaves in the sun in each canopy layer, pft, 

     ! and leaf layer. m2/m2
     real(r8),allocatable ::  tr_soil_dir(:)                              ! fraction of incoming direct  radiation that (cm_numSWb)
     ! is transmitted to the soil as direct
     real(r8),allocatable ::  tr_soil_dif(:)                              ! fraction of incoming diffuse radiation that 
     ! is transmitted to the soil as diffuse
     real(r8),allocatable ::  tr_soil_dir_dif(:)                          ! fraction of incoming direct  radiation that 
     ! is transmitted to the soil as diffuse
     real(r8),allocatable ::  fab(:)                                      ! fraction of incoming total   radiation that is absorbed by the canopy
     real(r8),allocatable ::  fabd(:)                                     ! fraction of incoming direct  radiation that is absorbed by the canopy
     real(r8),allocatable ::  fabi(:)                                     ! fraction of incoming diffuse radiation that is absorbed by the canopy
     real(r8),allocatable ::  sabs_dir(:)                                 ! fraction of incoming direct  radiation that is absorbed by the canopy
     real(r8),allocatable ::  sabs_dif(:)                                 ! fraction of incoming diffuse radiation that is absorbed by the canopy


     !SEED BANK
     real(r8) :: seeds_in(maxpft)                               ! seed production KgC/m2/year
     real(r8) :: seed_decay(maxpft)                             ! seed decay in KgC/m2/year
     real(r8) :: seed_germination(maxpft)                       ! germination rate of seed pool in KgC/m2/year

     ! PHOTOSYNTHESIS       

     real(r8) ::  psn_z(nclmax,maxpft,nlevleaf)               ! carbon assimilation in each canopy layer, pft, and leaf layer. umolC/m2/s
!     real(r8) ::  gpp                                              ! total patch gpp: KgC/m2/year
!     real(r8) ::  npp                                              ! total patch npp: KgC/m2/year   

     ! ROOTS
     real(r8), allocatable ::  rootfr_ft(:,:)                      ! root fraction of each PFT in each soil layer:-
     real(r8), allocatable ::  rootr_ft(:,:)                       ! fraction of water taken from each PFT and soil layer:-
     real(r8) ::  btran_ft(maxpft)                              ! btran calculated seperately for each PFT:-   

     ! DISTURBANCE 
     real(r8) ::  disturbance_rates(n_dist_types)                  ! disturbance rate from 1) mortality 
                                                                   !                       2) fire: fraction/day 
                                                                   !                       3) logging mortatliy
     real(r8) ::  disturbance_rate                                 ! larger effective disturbance rate: fraction/day

     ! LITTER AND COARSE WOODY DEBRIS 
     ! Pools of litter (non respiring) 
     real(r8) ::  cwd_ag(ncwd)                                     ! above ground coarse wood debris litter that does not respire. KgC/m2
     real(r8) ::  cwd_bg(ncwd)                                     ! below ground coarse wood debris litter that does not respire. KgC/m2
     real(r8) ::  leaf_litter(maxpft)                           ! above ground leaf litter that does not respire. KgC/m2
     real(r8) ::  root_litter(maxpft)                           ! below ground fine root litter that does not respire. KgC/m2

     ! Fluxes of litter (non respiring) 
     real(r8) :: fragmentation_scaler                              ! Scale rate of litter fragmentation. 0 to 1.
     real(r8) :: cwd_ag_in(ncwd)                                   ! Flux into CWD_AG from turnover and mortality KgC/m2/y
     real(r8) :: cwd_bg_in(ncwd)                                   ! Flux into cwd_bg from root turnover and mortality KgC/m2/y
     real(r8) :: cwd_ag_out(ncwd)                                  ! Flux out of AG CWD into AG litter KgC/m2/y
     real(r8) :: cwd_bg_out(ncwd)                                  ! Flux out of BG CWD into BG litter KgC/m2/


     real(r8) :: leaf_litter_in(maxpft)                         ! Flux in  to AG leaf litter from leaf turnover and mortality KgC/m2/y
     real(r8) :: leaf_litter_out(maxpft)                        ! Flux out of AG leaf litter from fragmentation KgC/m2/y
     real(r8) :: root_litter_in(maxpft)                         ! Flux in  to BG root litter from leaf turnover and mortality KgC/m2/y
     real(r8) :: root_litter_out(maxpft)                        ! Flux out of BG root from fragmentation KgC/m2/y

     ! Derivatives of litter (non respiring) 
     real(r8) ::  dcwd_AG_dt(ncwd)                                 ! rate of change of above ground CWD in each size class: KgC/m2/year. 
     real(r8) ::  dcwd_BG_dt(ncwd)                                 ! rate of change of below ground CWD in each size class: KgC/m2/year. 
     real(r8) ::  dleaf_litter_dt(maxpft)                       ! rate of change of leaf litter in each size class: KgC/m2/year. 
     real(r8) ::  droot_litter_dt(maxpft)                       ! rate of change of root litter in each size class: KgC/m2/year. 

     real(r8) ::  repro(maxpft)                                 ! allocation to reproduction per PFT : KgC/m2

     !FUEL CHARECTERISTICS
     real(r8) ::  sum_fuel                                         ! total ground fuel related to ros (omits 1000hr fuels): KgC/m2
     real(r8) ::  fuel_frac(nfsc)                                  ! fraction of each litter class in the ros_fuel:-.  
     real(r8) ::  livegrass                                        ! total aboveground grass biomass in patch.  KgC/m2
     real(r8) ::  fuel_bulkd                                       ! average fuel bulk density of the ground fuel 
                                                                   ! (incl. live grasses. omits 1000hr fuels). KgC/m3
     real(r8) ::  fuel_sav                                         ! average surface area to volume ratio of the ground fuel 
                                                                   ! (incl. live grasses. omits 1000hr fuels).
     real(r8) ::  fuel_mef                                         ! average moisture of extinction factor 
                                                                   ! of the ground fuel (incl. live grasses. omits 1000hr fuels).
     real(r8) ::  fuel_eff_moist                                   ! effective avearage fuel moisture content of the ground fuel 
                                                                   ! (incl. live grasses. omits 1000hr fuels)
     real(r8) ::  litter_moisture(nfsc)

     ! FIRE SPREAD
     real(r8) ::  ros_front                                        ! rate of forward  spread of fire: m/min
     real(r8) ::  ros_back                                         ! rate of backward spread of fire: m/min
     real(r8) ::  effect_wspeed                                    ! windspeed modified by fraction of relative grass and tree cover: m/min
     real(r8) ::  tau_l                                            ! Duration of lethal heating: mins
     real(r8) ::  fi                                               ! average fire intensity of flaming front:  kj/m/s or kw/m
     integer  ::  fire                                             ! Is there a fire? 1=yes 0=no
     real(r8) ::  fd                                               ! fire duration: mins
     real(r8) ::  nf                                               ! number of fires initiated daily: n/gridcell/day
     real(r8) ::  sh                                               ! average scorch height: m 

     ! FIRE EFFECTS     
     real(r8) ::  ab                                               ! area burnt:  m2/day
     real(r8) ::  frac_burnt                                       ! fraction burnt: frac gridcell/day  
     real(r8) ::  tfc_ros                                          ! total fuel consumed - no trunks.  KgC/m2/day
     real(r8) ::  burnt_frac_litter(nfsc)                          ! fraction of each litter pool burned:-


     ! PLANT HYDRAULICS     
     type(ed_patch_hydr_type) , pointer :: pa_hydr                 ! All patch hydraulics data, see FatesHydraulicsMemMod.F90

   contains

  end type ed_patch_type

  
  !************************************
  !** Resources management type      **
  ! YX
  !************************************
  type ed_resources_management_type
    
     real(r8) ::  trunk_product_site                       ! Actual  trunk product at site level KgC/site

     !debug variables
     real(r8) ::  delta_litter_stock
     real(r8) ::  delta_biomass_stock
     real(r8) ::  delta_individual
  
  end type ed_resources_management_type



  !************************************
  !** Site type structure           **
  !************************************

  type ed_site_type
     
     ! POINTERS  
     type (ed_patch_type), pointer :: oldest_patch => null()   ! pointer to oldest patch at the site  
     type (ed_patch_type), pointer :: youngest_patch => null() ! pointer to yngest patch at the site
     
     ! Resource management
     type (ed_resources_management_type) :: resources_management ! resources_management at the site 



     ! INDICES 
     real(r8) ::  lat                                          ! latitude:  degrees 
     real(r8) ::  lon                                          ! longitude: degrees 

     ! CARBON BALANCE       
     real(r8) :: flux_in                                      ! for carbon balance purpose. C coming into biomass pool:  KgC/site
     real(r8) :: flux_out                                     ! for carbon balance purpose. C leaving ED pools  KgC/site
     real(r8) :: old_stock                                    ! for accounting purposes, remember biomass stock from last time:  KgC/site
     real(r8) :: npp                                          ! used for calculating NEP and NBP during BGC summarization phase
     real(r8) :: nep                                          ! Net ecosystem production, i.e. fast-timescale carbon balance that 
                                                              ! does not include disturbance [gC/m2/s]
     real(r8) :: nbp                                          ! Net biosphere production, i.e. slow-timescale carbon balance that 
                                                              ! integrates to total carbon change [gC/m2/s]
     real(r8) :: tot_seed_rain_flux                           ! [gC/m2/s] total flux of carbon from seed rain
     real(r8) :: fire_c_to_atm                                ! total fire carbon loss to atmosphere [gC/m2/s]
     real(r8) :: ed_litter_stock                              ! litter in [gC/m2]
     real(r8) :: cwd_stock                                    ! coarse woody debris [gC/m2]
     real(r8) :: biomass_stock                                ! total biomass at the column level in [gC / m2]
     real(r8) :: totfatesc                                    ! Total FATES carbon at the site, including vegetation, CWD, seeds, 
                                                              ! and FATES portion of litter [gC/m2] 
     real(r8) :: totbgcc                                      ! Total BGC carbon at the site, including litter, and soil pools [gC/m2] 
     real(r8) :: totecosysc                                   ! Total ecosystem C at the site, including vegetation, 
                                                              ! CWD, litter (from HLM and FATES), and soil pools [gC/m2]

     real(r8) :: totfatesc_old                                ! Total FATES C at the site from last call to balance check [gC/m2]
     real(r8) :: totbgcc_old                                  ! Total BGC C at the site from last call to balance check [gC/m2] 
     real(r8) :: totecosysc_old                               ! Total ecosystem C at the site from last call to balance check [gC/m2]
     
     real(r8) :: fates_to_bgc_this_ts                         ! total flux of carbon from FATES to BGC models on current timestep [gC/m2/s] 
     real(r8) :: fates_to_bgc_last_ts                         ! total flux of carbon from FATES to BGC models on previous timestep [gC/m2/s] 

     real(r8) :: cbal_err_fates                               ! [gC/m2/s]  total carbon balance error for FATES processes
     real(r8) :: cbal_err_bgc                                 ! [gC/m2/s]  total carbon balance error for BGC (HLM) processes
     real(r8) :: cbal_err_tot                                 ! [gC/m2/s]  total carbon balance error for all land processes

     real(r8) :: nep_timeintegrated                           ! Net ecosystem production accumulated over model time-steps [gC/m2]
     real(r8) :: hr_timeintegrated                            ! Heterotrophic respiration accumulated over model time-steps [gC/m2]
     real(r8) :: npp_timeintegrated                           ! Net primary production accumulated over model time-steps [gC/m2]
     real(r8) :: nbp_integrated                               ! Net biosphere production accumulated over model time-steps [gC/m2]


     ! PHENOLOGY 
     real(r8) ::  ED_GDD_site                                  ! ED Phenology growing degree days.
     integer  ::  status                                       ! are leaves in this pixel on or off for cold decid
     integer  ::  dstatus                                      ! are leaves in this pixel on or off for drought decid
     real(r8) ::  ncd                                          ! no chilling days:-
     real(r8) ::  last_n_days(senes)                           ! record of last 10 days temperature for senescence model. deg C
     integer  ::  leafondate                                   ! doy of leaf on:-
     integer  ::  leafoffdate                                  ! doy of leaf off:-
     integer  ::  dleafondate                                  ! doy of leaf on drought:-
     integer  ::  dleafoffdate                                 ! doy of leaf on drought:-
     real(r8) ::  water_memory(numWaterMem)                             ! last 10 days of soil moisture memory...

     !SEED BANK
     real(r8) :: seed_bank(maxpft)                              ! seed pool in KgC/m2/year
     real(r8) :: dseed_dt(maxpft)
     real(r8) :: seed_rain_flux(maxpft)                         ! flux of seeds from exterior KgC/m2/year (needed for C balance purposes)

     ! FIRE
     real(r8) ::  wind                                         ! daily wind in m/min for Spitfire units 
     real(r8) ::  acc_ni                                       ! daily nesterov index accumulating over time.
     real(r8) ::  fdi                                          ! daily probability an ignition event will start a fire
     real(r8) ::  frac_burnt                                   ! fraction of soil burnt in this day.
     real(r8) ::  total_burn_flux_to_atm                       ! total carbon burnt to the atmosphere in this day. KgC/site
     real(r8) ::  cwd_ag_burned(ncwd)
     real(r8) ::  leaf_litter_burned(maxpft)

     ! PLANT HYDRAULICS
     type(ed_site_hydr_type), pointer :: si_hydr
        
     ! TERMINATION, RECRUITMENT, DEMOTION, and DISTURBANCE

     real(r8), allocatable :: terminated_nindivs(:,:,:) ! number of individuals that were in cohorts which were terminated this timestep, on size x pft x canopy array. 
     real(r8) :: termination_carbonflux(2)                     ! carbon flux from live to dead pools associated with termination mortality, per canopy level
     real(r8) :: recruitment_rate(1:maxpft)                     ! number of individuals that were recruited into new cohorts
     real(r8), allocatable :: demotion_rate(:)                ! rate of individuals demoted from canopy to understory per FATES timestep
     real(r8) :: demotion_carbonflux                           ! biomass of demoted individuals from canopy to understory [kgC/ha/day]
     real(r8), allocatable :: promotion_rate(:)               ! rate of individuals promoted from understory to canopy per FATES timestep
     real(r8) :: promotion_carbonflux                          ! biomass of promoted individuals from understory to canopy [kgC/ha/day]

     ! some diagnostic-only (i.e. not resolved by ODE solver) flux of carbon to CWD and litter pools from termination and canopy mortality
     real(r8) :: CWD_AG_diagnostic_input_carbonflux(1:ncwd)       ! diagnostic flux to AG CWD [kg C / m2 / yr]
     real(r8) :: CWD_BG_diagnostic_input_carbonflux(1:ncwd)       ! diagnostic flux to BG CWD [kg C / m2 / yr]
     real(r8) :: leaf_litter_diagnostic_input_carbonflux(1:maxpft) ! diagnostic flux to AG litter [kg C / m2 / yr]
     real(r8) :: root_litter_diagnostic_input_carbonflux(1:maxpft) ! diagnostic flux to BG litter [kg C / m2 / yr]

     ! Canopy Spread
     real(r8) ::  spread                                          ! dynamic canopy allometric term [unitless]
     
  end type ed_site_type

contains

  ! =====================================================================================

  subroutine val_check_ed_vars(currentPatch,var_aliases,return_code)

     ! ----------------------------------------------------------------------------------
     ! Perform numerical checks on variables of interest.
     ! The input string is of the form:  'VAR1_NAME:VAR2_NAME:VAR3_NAME'
     ! ----------------------------------------------------------------------------------


     use FatesUtilsMod,only : check_hlm_list
     use FatesUtilsMod,only : check_var_real

     ! Arguments
     type(ed_patch_type),intent(in), target :: currentPatch
     character(len=*),intent(in)            :: var_aliases
     integer,intent(out)                    :: return_code ! return 0 for all fine
                                                           ! return 1 if a nan detected
                                                           ! return 10+ if an overflow
                                                           ! return 100% if an underflow
     ! Locals
     type(ed_cohort_type), pointer          :: currentCohort

     
     ! Check through a registry of variables to check
     
     if ( check_hlm_list(trim(var_aliases),'co_n') ) then

        currentCohort => currentPatch%shortest
        do while(associated(currentCohort))
           call check_var_real(currentCohort%n,'cohort%n',return_code)
           if(.not.(return_code.eq.0)) then
              call dump_site(currentPatch%siteptr)
              call dump_patch(currentPatch)
              call dump_cohort(currentCohort)
              return
           end if
           currentCohort => currentCohort%taller
        end do
     end if
     
     if ( check_hlm_list(trim(var_aliases),'co_dbh') ) then

        currentCohort => currentPatch%shortest
        do while(associated(currentCohort))        
           call check_var_real(currentCohort%dbh,'cohort%dbh',return_code)
           if(.not.(return_code.eq.0)) then
              call dump_site(currentPatch%siteptr)
              call dump_patch(currentPatch)
              call dump_cohort(currentCohort)
              return
           end if
           currentCohort => currentCohort%taller
        end do
     end if

     if ( check_hlm_list(trim(var_aliases),'pa_area') ) then

        call check_var_real(currentPatch%area,'patch%area',return_code)
        if(.not.(return_code.eq.0)) then
           call dump_site(currentPatch%siteptr)
           call dump_patch(currentPatch)
           return
        end if
     end if
     


     return
  end subroutine val_check_ed_vars

  ! =====================================================================================

  subroutine dump_site(csite) 

     type(ed_site_type),intent(in),target :: csite


     ! EDTypes is 

     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) ' Site Coordinates                       '
     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) 'latitude                    = ', csite%lat
     write(fates_log(),*) 'longitude                   = ', csite%lon
     write(fates_log(),*) '----------------------------------------'
     return

  end subroutine dump_site

  ! =====================================================================================


  subroutine dump_patch(cpatch)

     type(ed_patch_type),intent(in),target :: cpatch

     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) ' Dumping Patch Information              '
     write(fates_log(),*) ' (omitting arrays)                      '
     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) 'pa%patchno            = ',cpatch%patchno
     write(fates_log(),*) 'pa%age                = ',cpatch%age
     write(fates_log(),*) 'pa%age_class          = ',cpatch%age_class
     write(fates_log(),*) 'pa%area               = ',cpatch%area
     write(fates_log(),*) 'pa%countcohorts       = ',cpatch%countcohorts
     write(fates_log(),*) 'pa%ncl_p              = ',cpatch%ncl_p
     write(fates_log(),*) 'pa%total_canopy_area  = ',cpatch%total_canopy_area
     write(fates_log(),*) 'pa%total_tree_area    = ',cpatch%total_tree_area
     write(fates_log(),*) 'pa%canopy_area        = ',cpatch%canopy_area
     write(fates_log(),*) 'pa%bare_frac_area     = ',cpatch%bare_frac_area
     write(fates_log(),*) 'pa%lai                = ',cpatch%lai
     write(fates_log(),*) 'pa%zstar              = ',cpatch%zstar
     write(fates_log(),*) 'pa%disturbance_rate   = ',cpatch%disturbance_rate
     write(fates_log(),*) '----------------------------------------'
     return

  end subroutine dump_patch

  ! =====================================================================================
  
  subroutine dump_cohort(ccohort)


     type(ed_cohort_type),intent(in),target :: ccohort
     
     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) ' Dumping Cohort Information             '
     write(fates_log(),*) '----------------------------------------'
     write(fates_log(),*) 'co%pft                    = ', ccohort%pft
     write(fates_log(),*) 'co%n                      = ', ccohort%n                         
     write(fates_log(),*) 'co%dbh                    = ', ccohort%dbh                                        
     write(fates_log(),*) 'co%hite                   = ', ccohort%hite                                
     write(fates_log(),*) 'co%b                      = ', ccohort%b                            
     write(fates_log(),*) 'co%balive                 = ', ccohort%balive
     write(fates_log(),*) 'co%bdead                  = ', ccohort%bdead                          
     write(fates_log(),*) 'co%bstore                 = ', ccohort%bstore
     write(fates_log(),*) 'co%laimemory              = ', ccohort%laimemory
     write(fates_log(),*) 'co%bsw                    = ', ccohort%bsw                  
     write(fates_log(),*) 'co%bl                     = ', ccohort%bl
     write(fates_log(),*) 'co%br                     = ', ccohort%br
     write(fates_log(),*) 'co%lai                    = ', ccohort%lai                         
     write(fates_log(),*) 'co%sai                    = ', ccohort%sai  
     write(fates_log(),*) 'co%gscan                  = ', ccohort%gscan
     write(fates_log(),*) 'co%leaf_cost              = ', ccohort%leaf_cost
     write(fates_log(),*) 'co%canopy_layer           = ', ccohort%canopy_layer
     write(fates_log(),*) 'co%canopy_layer_yesterday = ', ccohort%canopy_layer_yesterday
     write(fates_log(),*) 'co%nv                     = ', ccohort%nv
     write(fates_log(),*) 'co%status_coh             = ', ccohort%status_coh
     write(fates_log(),*) 'co%canopy_trim            = ', ccohort%canopy_trim
     write(fates_log(),*) 'co%status_coh             = ', ccohort%status_coh               
     write(fates_log(),*) 'co%excl_weight            = ', ccohort%excl_weight               
     write(fates_log(),*) 'co%prom_weight            = ', ccohort%prom_weight               
     write(fates_log(),*) 'co%size_class             = ', ccohort%size_class
     write(fates_log(),*) 'co%size_by_pft_class      = ', ccohort%size_by_pft_class
     write(fates_log(),*) 'co%gpp_acc_hold           = ', ccohort%gpp_acc_hold
     write(fates_log(),*) 'co%gpp_acc                = ', ccohort%gpp_acc
     write(fates_log(),*) 'co%gpp_tstep              = ', ccohort%gpp_tstep
     write(fates_log(),*) 'co%npp_acc_hold           = ', ccohort%npp_acc_hold
     write(fates_log(),*) 'co%npp_tstep              = ', ccohort%npp_tstep
     write(fates_log(),*) 'co%npp_acc                = ', ccohort%npp_acc
     write(fates_log(),*) 'co%resp_tstep             = ', ccohort%resp_tstep
     write(fates_log(),*) 'co%resp_acc               = ', ccohort%resp_acc
     write(fates_log(),*) 'co%resp_acc_hold          = ', ccohort%resp_acc_hold
     write(fates_log(),*) 'co%npp_leaf               = ', ccohort%npp_leaf
     write(fates_log(),*) 'co%npp_froot              = ', ccohort%npp_froot
     write(fates_log(),*) 'co%npp_bsw                = ', ccohort%npp_bsw
     write(fates_log(),*) 'co%npp_bdead              = ', ccohort%npp_bdead
     write(fates_log(),*) 'co%npp_bseed              = ', ccohort%npp_bseed
     write(fates_log(),*) 'co%npp_store              = ', ccohort%npp_store
     write(fates_log(),*) 'co%rdark                  = ', ccohort%rdark
     write(fates_log(),*) 'co%resp_m                 = ', ccohort%resp_m
     write(fates_log(),*) 'co%resp_g                 = ', ccohort%resp_g
     write(fates_log(),*) 'co%livestem_mr            = ', ccohort%livestem_mr
     write(fates_log(),*) 'co%livecroot_mr           = ', ccohort%livecroot_mr
     write(fates_log(),*) 'co%froot_mr               = ', ccohort%froot_mr
     write(fates_log(),*) 'co%md                     = ', ccohort%md
     write(fates_log(),*) 'co%leaf_md                = ', ccohort%leaf_md
     write(fates_log(),*) 'co%root_md                = ', ccohort%root_md
     write(fates_log(),*) 'co%carbon_balance         = ', ccohort%carbon_balance
     write(fates_log(),*) 'co%dmort                  = ', ccohort%dmort
     write(fates_log(),*) 'co%seed_prod              = ', ccohort%seed_prod
     write(fates_log(),*) 'co%treelai                = ', ccohort%treelai
     write(fates_log(),*) 'co%treesai                = ', ccohort%treesai
     write(fates_log(),*) 'co%leaf_litter            = ', ccohort%leaf_litter
     write(fates_log(),*) 'co%c_area                 = ', ccohort%c_area
     write(fates_log(),*) 'co%woody_turnover         = ', ccohort%woody_turnover
     write(fates_log(),*) 'co%cmort                  = ', ccohort%cmort
     write(fates_log(),*) 'co%bmort                  = ', ccohort%bmort
     write(fates_log(),*) 'co%imort                  = ', ccohort%imort
     write(fates_log(),*) 'co%fmort                  = ', ccohort%fmort
     write(fates_log(),*) 'co%hmort                  = ', ccohort%hmort
     write(fates_log(),*) 'co%isnew                  = ', ccohort%isnew
     write(fates_log(),*) 'co%dndt                   = ', ccohort%dndt
     write(fates_log(),*) 'co%dhdt                   = ', ccohort%dhdt
     write(fates_log(),*) 'co%ddbhdt                 = ', ccohort%ddbhdt
     write(fates_log(),*) 'co%dbalivedt              = ', ccohort%dbalivedt
     write(fates_log(),*) 'co%dbdeaddt               = ', ccohort%dbdeaddt
     write(fates_log(),*) 'co%dbstoredt              = ', ccohort%dbstoredt
     write(fates_log(),*) 'co%storage_flux           = ', ccohort%storage_flux
     write(fates_log(),*) 'co%cfa                    = ', ccohort%cfa
     write(fates_log(),*) 'co%fire_mort              = ', ccohort%fire_mort
     write(fates_log(),*) 'co%crownfire_mort         = ', ccohort%crownfire_mort
     write(fates_log(),*) 'co%cambial_mort           = ', ccohort%cambial_mort
     write(fates_log(),*) 'co%size_class             = ', ccohort%size_class
     write(fates_log(),*) 'co%size_by_pft_class      = ', ccohort%size_by_pft_class
     write(fates_log(),*) '----------------------------------------'
     return
  end subroutine dump_cohort

end module EDTypesMod
