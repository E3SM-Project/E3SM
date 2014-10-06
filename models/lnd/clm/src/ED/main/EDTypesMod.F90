  module EDtypesMod

     use shr_kind_mod , only : r8 => shr_kind_r8;
     use clm_varpar   , only : nlevcan_ed, nclmax, numrad, nlevgrnd
     use domainMod    , only : domain_type

     implicit none
     save

                                                           !SWITCHES THAT ARE READ IN
  integer         RESTART                                  ! restart flag, 1= read initial system state 0 = bare ground

                                                           ! MODEL PARAMETERS
  real(r8)            :: timestep_secs                     ! subdaily timestep in seconds (e.g. 1800 or 3600) 
  integer             :: n_sub                             ! num of substeps in year 
  real(r8), parameter :: AREA                 = 10000.0_r8 ! Notional area of simulated forest m2
  integer  doy

  integer, parameter  :: invalidValue         = -9999      ! invalid value for gcells,
                                                           ! cohorts, and patches

                                                           ! for setting number of patches per gridcell and number of cohorts per patch
                                                           ! for I/O and converting to a vector
  integer, parameter :: numPatchesPerGridCell = 4          !
  integer, parameter :: numCohortsPerPatch    = 20         !
  integer, parameter :: cohorts_per_gcell     = 80         ! should be numPatchesPerGridCell*numCohortsPerPatch
                                                           ! portions of the pft 
  integer, parameter :: numWaterMem           = 10         ! watermemory saved as site
                                                           ! level var

                                                           ! BIOLOGY/BIOGEOCHEMISTRY        
  integer , parameter :: INTERNAL_RECRUITMENT = 1          ! internal recruitment fla  1=yes  
  integer , parameter :: EXTERNAL_RECRUITMENT = 0          ! external recruitment flag 1=yes  
  integer , parameter :: SENES                = 10         ! Window of time over which we track temp for cold sensecence (days)
  real(r8), parameter :: DINC_ED              = 1.0_r8     ! size of LAI bins. 
  integer , parameter :: N_DIST_TYPES         = 2          ! number of disturbance types (mortality, fire)
  integer , parameter :: numpft_ed            = 2          ! number of PFTs used in ED. 
  
                                                           ! SPITFIRE     
  integer , parameter :: NLSC                 = 5          ! number carbon compartments in above ground litter array 
  integer , parameter :: NFSC                 = 6          ! number fuel size classes  
  integer , parameter :: N_EF                 = 7          ! number of emission factors. One per trace gas or aerosol species.
  integer,  parameter :: NCWD                 = 4          ! number of coarse woody debris pools
  integer,  parameter :: lg_sf                = 6          ! array index of live grass pool for spitfire
  integer,  parameter :: dg_sf                = 1          ! array index of dead grass pool for spitfire
  integer,  parameter :: tr_sf                = 5          ! array index of dead trunk pool for spitfire
  integer,  parameter :: lb_sf                = 4          ! array index of lrge branch pool for spitfire 
  real(r8), parameter :: fire_threshold       = 35.0_r8    ! threshold for fires that spread or go out. KWm-2
 
                                                           ! COHORT FUSION          
  real(r8), parameter :: FUSETOL              = 0.6_r8     ! min fractional difference in dbh between cohorts

                                                           ! PATCH FUSION 
  real(r8), parameter :: NTOL                 = 0.05_r8    ! min plant density for hgt bin to be used in height profile comparisons 
  real(r8), parameter :: HITEMAX              = 30.0_r8    ! max dbh value used in hgt profile comparison 
  real(r8), parameter :: DBHMAX               = 150.0_r8   ! max dbh value used in hgt profile comparison 
  integer , parameter :: N_HITE_BINS          = 60         ! no. of hite bins used to distribute LAI
  integer , parameter :: N_DBH_BINS           = 5          ! no. of dbh bins used when comparing patches

  character*4 yearchar                     

                                                         !************************************
                                                         !** COHORT type structure          **
                                                         !************************************
  type cohort

     ! POINTERS
     type (cohort) , pointer :: taller   => null()       ! pointer to next tallest cohort     
     type (cohort) , pointer :: shorter  => null()       ! pointer to next shorter cohort     
     type (patch)  , pointer :: patchptr => null()       ! pointer to patch that cohort is in
     type (site)   , pointer :: siteptr  => null()       ! pointer to site that cohort is in

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

     ! CARBON FLUXES 
     real(r8) ::  gpp                                    ! GPP:  kgC/indiv/year
     real(r8) ::  gpp_acc                                ! GPP:  kgC/indiv/day   
     real(r8) ::  gpp_clm                                ! GPP:  kgC/indiv/timestep
     real(r8) ::  npp                                    ! NPP:  kgC/indiv/year
     real(r8) ::  npp_acc                                ! NPP:  kgC/indiv/day   
     real(r8) ::  npp_clm                                ! NPP:  kgC/indiv/timestep
     real(r8) ::  resp                                   ! Resp: kgC/indiv/year
     real(r8) ::  resp_acc                               ! Resp: kgC/indiv/day
     real(r8) ::  resp_clm                               ! Resp: kgC/indiv/timestep

     real(r8) ::  ts_net_uptake(nlevcan_ed)              ! Net uptake of leaf layers: kgC/m2/s
     real(r8) ::  year_net_uptake(nlevcan_ed)            ! Net uptake of leaf layers: kgC/m2/year

     ! RESPIRATION COMPONENTS
     real(r8) ::  rd                                     ! Dark respiration: umol/indiv/s
     real(r8) ::  resp_g                                 ! Growth respiration:  kgC/indiv/timestep
     real(r8) ::  resp_m                                 ! Maintenance respiration:  kgC/indiv/timestep 
     real(r8) ::  livestem_mr                            ! Live stem        maintenance respiration: kgC/indiv/s
     real(r8) ::  livecroot_mr                           ! Live coarse root maintenance respiration: kgC/indiv/s
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

     ! NITROGEN POOLS      
     real(r8) ::  livestemn                              ! live stem nitrogen       : KgN/invid
     real(r8) ::  livecrootn                             ! live coarse root nitrogen: KgN/invid
     real(r8) ::  frootn                                 ! fine root  nitrogen      : KgN/invid

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

  end type cohort

  !************************************
  !** Patch type structure           **
  !************************************
  type patch

     ! POINTERS
     type (cohort), pointer :: tallest => null()                   ! pointer to patch's tallest cohort    
     type (cohort), pointer :: shortest => null()                  ! pointer to patch's shortest cohort
     type (patch),  pointer :: older => null()                     ! pointer to next older patch   
     type (patch),  pointer :: younger => null()                   ! pointer to next younger patch      
     type (site),   pointer :: siteptr => null()                   ! pointer to the site that the patch is in

     !INDICES
     integer  :: patchno                                           ! unique number given to each new patch created for tracking
     integer  :: clm_pno                                           ! clm patch number (index of p vector)

     ! PATCH INFO
     real(r8) ::  age                                              ! average patch age: years                   
     real(r8) ::  area                                             ! patch area: m2  
     integer  ::  countcohorts                                     ! Number of cohorts in patch
     integer  ::  ncl_p                                            ! Number of occupied canopy layers

     ! LEAF ORGANIZATION
     real(r8) ::  spread(nclmax)                                   ! dynamic ratio of dbh to canopy area: cm/m2
     real(r8) ::  pft_agb_profile(numpft_ed,n_dbh_bins)            ! binned above ground biomass, for patch fusion: KgC/m2
     real(r8) ::  canopy_layer_lai(nclmax)                         ! lai that is shading this canopy layer: m2/m2 
     real(r8) ::  total_canopy_area                                ! area that is covered by vegetation : m2
     real(r8) ::  total_tree_area                                  ! area that is covered by woody vegetation : m2
     real(r8) ::  canopy_area                                      ! area that is covered by vegetation : m2 (is this different to total_canopy_area?
     real(r8) ::  bare_frac_area                                   ! bare soil in this patch expressed as a fraction of the total soil surface.
     real(r8) ::  lai                                              ! leaf area index of patch

     real(r8) ::  tlai_profile(nclmax,numpft_ed,nlevcan_ed)        ! total   leaf area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  elai_profile(nclmax,numpft_ed,nlevcan_ed)        ! exposed leaf area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  tsai_profile(nclmax,numpft_ed,nlevcan_ed)        ! total   stem area in each canopy layer, pft, and leaf layer. m2/m2
     real(r8) ::  esai_profile(nclmax,numpft_ed,nlevcan_ed)        ! exposed stem area in each canopy layer, pft, and leaf layer. m2/m2

     real(r8) ::  canopy_area_profile(nclmax,numpft_ed,nlevcan_ed) ! fraction of canopy in each canopy 
                                                                   ! layer, pft, and leaf layer:-
     integer  ::  present(nclmax,numpft_ed)                        ! is there any of this pft in this canopy layer?      
     integer  ::  nrad(nclmax,numpft_ed)                           ! number of exposed leaf layers for each canopy layer and pft
     integer  ::  ncan(nclmax,numpft_ed)                           ! number of total   leaf layers for each canopy layer and pft

     !RADIATION FLUXES      
     real(r8) ::  fabd_sun_z(nclmax,numpft_ed,nlevcan_ed)          ! sun fraction of direct light absorbed by each canopy 
                                                                   ! layer, pft, and leaf layer:-
     real(r8) ::  fabd_sha_z(nclmax,numpft_ed,nlevcan_ed)          ! shade fraction of direct light absorbed by each canopy 
                                                                   ! layer, pft, and leaf layer:-
     real(r8) ::  fabi_sun_z(nclmax,numpft_ed,nlevcan_ed)          ! sun fraction of indirect light absorbed by each canopy 
                                                                   ! layer, pft, and leaf layer:-
     real(r8) ::  fabi_sha_z(nclmax,numpft_ed,nlevcan_ed)          ! shade fraction of indirect light absorbed by each canopy 
                                                                   ! layer, pft, and leaf layer:-

     real(r8) ::  ed_laisun_z(nclmax,numpft_ed,nlevcan_ed)         ! amount of LAI in the sun   in each canopy layer, 
                                                                   ! pft, and leaf layer. m2/m2
     real(r8) ::  ed_laisha_z(nclmax,numpft_ed,nlevcan_ed)         ! amount of LAI in the shade in each canopy layer,
     real(r8) ::  ed_parsun_z(nclmax,numpft_ed,nlevcan_ed)         ! PAR absorbed  in the sun   in each canopy layer,
     real(r8) ::  ed_parsha_z(nclmax,numpft_ed,nlevcan_ed)         ! PAR absorbed  in the shade in each canopy layer,
     real(r8) ::  f_sun(nclmax,numpft_ed,nlevcan_ed)               ! fraction of leaves in the sun in each canopy layer, pft, 
                                                                   ! and leaf layer. m2/m2
     real(r8) ::  tr_soil_dir(numrad)                              ! fraction of incoming direct  radiation that 
                                                                   ! is transmitted to the soil as direct
     real(r8) ::  tr_soil_dif(numrad)                              ! fraction of incoming diffuse radiation that 
                                                                   ! is transmitted to the soil as diffuse
     real(r8) ::  tr_soil_dir_dif(numrad)                          ! fraction of incoming direct  radiation that 
                                                                   ! is transmitted to the soil as diffuse
     real(r8) ::  fab(numrad)                                      ! fraction of incoming total   radiation that is absorbed by the canopy
     real(r8) ::  fabd(numrad)                                     ! fraction of incoming direct  radiation that is absorbed by the canopy
     real(r8) ::  fabi(numrad)                                     ! fraction of incoming diffuse radiation that is absorbed by the canopy
     real(r8) ::  sabs_dir(numrad)                                 ! fraction of incoming direct  radiation that is absorbed by the canopy
     real(r8) ::  sabs_dif(numrad)                                 ! fraction of incoming diffuse radiation that is absorbed by the canopy


     !SEED BANK
     real(r8) :: seed_bank(numpft_ed)                              ! seed pool in KgC/m2/year
     real(r8) :: seeds_in(numpft_ed)                               ! seed production KgC/m2/year
     real(r8) :: seed_decay(numpft_ed)                             ! seed decay in KgC/m2/year
     real(r8) :: seed_germination(numpft_ed)                       ! germination rate of seed pool in KgC/m2/year
     real(r8) :: dseed_dt(numpft_ed)

     ! PHOTOSYNTHESIS       
     real(r8) ::  psn_z(nclmax,numpft_ed,nlevcan_ed)               ! carbon assimilation in each canopy layer, pft, and leaf layer. umolC/m2/s
     real(r8) ::  gpp                                              ! total patch gpp: KgC/m2/year
     real(r8) ::  npp                                              ! total patch npp: KgC/m2/year   

     ! ROOTS
     real(r8), allocatable ::  rootfr_ft(:,:)                      ! root fraction of each PFT in each soil layer:-
     real(r8), allocatable ::  rootr_ft(:,:)                       ! fraction of water taken from each PFT and soil layer:-
     real(r8) ::  btran_ft(numpft_ed)                              ! btran calculated seperately for each PFT:-   

     ! DISTURBANCE 
     real(r8) ::  disturbance_rates(n_dist_types)                  ! disturbance rate from 1) mortality and 2) fire: fraction/day
     real(r8) ::  disturbance_rate                                 ! larger effective disturbance rate: fraction/day

     ! LITTER AND COARSE WOODY DEBRIS 
                                                                   ! Pools of litter (non respiring) 
     real(r8) ::  cwd_ag(ncwd)                                     ! above ground coarse wood debris litter that does not respire. KgC/m2
     real(r8) ::  cwd_bg(ncwd)                                     ! below ground coarse wood debris litter that does not respire. KgC/m2
     real(r8) ::  leaf_litter(numpft_ed)                           ! above ground leaf litter that does not respire. KgC/m2
     real(r8) ::  root_litter(numpft_ed)                           ! below ground fine root litter that does not respire. KgC/m2

     ! Fluxes of litter (non respiring) 
     real(r8) :: fragmentation_scaler    ! Scale rate of litter fragmentation. 0 to 1.
     real(r8) :: cwd_ag_in(ncwd)                                   ! Flux into CWD_AG from turnover and mortality KgC/m2/y
     real(r8) :: cwd_bg_in(ncwd)                                   ! Flux into cwd_bg from root turnover and mortality KgC/m2/y
     real(r8) :: cwd_ag_out(ncwd)                                  ! Flux out of AG CWD into AG litter KgC/m2/y
     real(r8) :: cwd_bg_out(ncwd)                                  ! Flux out of BG CWD into BG litter KgC/m2/


     real(r8) :: leaf_litter_in(numpft_ed)                         ! Flux in  to AG leaf litter from leaf turnover and mortality KgC/m2/y
     real(r8) :: leaf_litter_out(numpft_ed)                        ! Flux out of AG leaf litter from fragmentation KgC/m2/y
     real(r8) :: root_litter_in(numpft_ed)                         ! Flux in  to BG root litter from leaf turnover and mortality KgC/m2/y
     real(r8) :: root_litter_out(numpft_ed)                        ! Flux out of BG root from fragmentation KgC/m2/y

     ! Derivatives of litter (non respiring) 
     real(r8) ::  dcwd_AG_dt(ncwd)                                 ! rate of change of above ground CWD in each size class: KgC/m2/year. 
     real(r8) ::  dcwd_BG_dt(ncwd)                                 ! rate of change of below ground CWD in each size class: KgC/m2/year. 
     real(r8) ::  dleaf_litter_dt(numpft_ed)                       ! rate of change of leaf litter in each size class: KgC/m2/year. 
     real(r8) ::  droot_litter_dt(numpft_ed)                       ! rate of change of root litter in each size class: KgC/m2/year. 

     real(r8) ::  canopy_mortality_woody_litter                    ! flux of wood litter in to litter pool: KgC/m2/year
     real(r8) ::  canopy_mortality_leaf_litter(numpft_ed)          ! flux in to  leaf litter from tree death: KgC/m2/year
     real(r8) ::  canopy_mortality_root_litter(numpft_ed)          ! flux in to froot litter  from tree death: KgC/m2/year

     real(r8) ::  repro(numpft_ed)                                 ! allocation to reproduction per PFT : KgC/m2

     !FUEL CHARECTERISTICS
     real(r8) ::  sum_fuel                                         ! total ground fuel related to ros (omits 1000hr fuels): KgC/m2
     real(r8) ::  fuel_frac(ncwd+2)                                ! fraction of each litter class in the ros_fuel:-.  
     real(r8) ::  livegrass                                        ! total aboveground grass biomass in patch.  KgC/m2
     real(r8) ::  fuel_bulkd                                       ! average fuel bulk density of the ground fuel 
                                                                   ! (incl. live grasses. omits 1000hr fuels). KgC/m3
     real(r8) ::  fuel_sav                                         ! average surface area to volume ratio of the ground fuel 
                                                                   ! (incl. live grasses. omits 1000hr fuels).
     real(r8) ::  fuel_mef                                         ! average moisture of extinction factor 
                                                                   ! of the ground fuel (incl. live grasses. omits 1000hr fuels).
     real(r8) ::  fuel_eff_moist                                   ! effective avearage fuel moisture content of the ground fuel 
                                                                   ! (incl. live grasses. omits 1000hr fuels)
     real(r8) ::  litter_moisture(ncwd+2)

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

  end type patch

  !************************************
  !** Site type structure           **
  !************************************
  type site
     ! POINTERS  
     type (patch), pointer :: oldest_patch => null()        ! pointer to oldest patch at the site  
     type (patch), pointer :: youngest_patch => null()      ! pointer to yngest patch at the site

     ! INDICES 
     real(r8) ::  lat                 ! latitude:  degrees 
     real(r8) ::  lon                 ! longitude: degrees 
     integer  ::  clmgcell            ! gridcell index
     integer  ::  clmcolumn           ! column index (assuming there is only one soil column in each gcell.
     integer  ::  istheresoil         ! are there any soil columns, or is this all ice/rocks/lakes?

     ! CARBON BALANCE       
     real(r8) ::  flux_in             ! for carbon balance purpose. C coming into biomass pool:  KgC/site
     real(r8) ::  flux_out            ! for carbon balance purpose. C leaving ED pools  KgC/site
     real(r8) ::  old_stock           ! for accounting purposes, remember biomass stock from last time:  KgC/site
 
     ! DISTURBANCE
     real(r8) ::  disturbance_mortality     ! site level disturbance rates from mortality.
     real(r8) ::  disturbance_fire          ! site level disturbance rates from fire.  
     integer  ::  dist_type                 ! disturbance dist_type id.
     real(r8) ::  disturbance_rate          ! site total dist rate

     ! PHENOLOGY 
     integer  ::  status              ! are leaves in this pixel on or off for cold decid
     integer  ::  dstatus             ! are leaves in this pixel on or off for drought decid
     real(r8) ::  gdd                 ! growing degree days: deg C. 
     real(r8) ::  ncd                 ! no chilling days:-
     real(r8) ::  last_n_days(senes)  ! record of last 10 days temperature for senescence model. deg C
     integer  ::  leafondate                ! doy of leaf on:-
     integer  ::  leafoffdate               ! doy of leaf off:-
     integer  ::  dleafondate               ! doy of leaf on drought:-
     integer  ::  dleafoffdate              ! doy of leaf on drought:-
     real(r8) ::  water_memory(10)          ! last 10 days of soil moisture memory...
     real(r8) ::  cwd_ag_burned(ncwd)
     real(r8) :: leaf_litter_burned(numpft_ed)
     
     ! FIRE 
     real(r8) ::  acc_ni              ! daily nesterov index accumulating over time.
     real(r8) ::  ab                  ! daily burnt area: m2
     real(r8) ::  frac_burnt          ! fraction of soil burnt in this day.
 
  end type site

  !************************************
  !** Userdata type structure       **
  !************************************

  type userdata
     ! POINTERS 
     type (site)  , pointer :: firstsite_pnt => null()         ! pointer to the first site in the system
     type (cohort), pointer :: storesmallcohort => null()      ! storage of the smallest cohort for insertion routine
     type (cohort), pointer :: storebigcohort => null()        ! storage of the largest cohort for insertion routine 

     ! TIMING PARAMETERS
     integer  ::   time_period              ! Within year timestep (1:N_SUB) day of year
     real(r8) ::   deltat                   ! fraction of year used for each timestep (1/N_SUB)
     integer  ::   currentindex             ! Counts up the number of cohorts which have been made. 
     integer  ::   restart_year             ! Which year of simulation are we starting in? 

  end type userdata

  !****************************************
  !** Pointer arrays for linking to CLM  **
  !****************************************               
  type gridcell_edstate_type
     type(site),pointer :: spnt  
  end type gridcell_edstate_type

  !
  ! cohortype.  mimics the types in clmtype.  For now this holds ED data that is
  ! necessary in the rest of CLM
  !

  type cohort_type
     integer , pointer :: gridcell(:) !index into gridcell level quantities
  end type cohort_type

  type(cohort_type)     , public :: coh

  !----------------------------------------------------
  ! Declare single instance of arrays
  !----------------------------------------------------
  type(userdata)             , public, target         :: udata
  type(gridcell_edstate_type),            allocatable :: gridCellEdState(:) !patch ED structure

end module EDtypesMod
