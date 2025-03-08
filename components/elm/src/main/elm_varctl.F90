module elm_varctl

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing run control variables
  !
  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8, SHR_KIND_CL, SHR_KIND_CS
  use shr_sys_mod , only: shr_sys_abort ! cannot use endrun here due to circular dependency

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  implicit none
  public :: elm_varctl_set    ! Set variables
  public :: cnallocate_carbon_only_set
  public :: cnallocate_carbon_only
  public :: cnallocate_carbonnitrogen_only_set
  public :: cnallocate_carbonnitrogen_only
  public :: cnallocate_carbonphosphorus_only_set
  public :: cnallocate_carbonphosphorus_only
  public :: get_carbontag ! get the tag for carbon simulations
  !
  private
  save
  !
  ! !PUBLIC TYPES:
  !
  integer , parameter, public ::  iundef = -9999999
  real(r8), parameter, public ::  rundef = -9999999._r8
  integer , parameter, public ::  fname_len = SHR_KIND_CL   ! max length of file names in this module
  !----------------------------------------------------------
  !
  ! Run control variables
  !
  ! case id
  character(len=256), public :: caseid  = ' '

  ! case title
  character(len=256), public :: ctitle  = ' '

  ! Type of run
  integer, public :: nsrest             = iundef

  ! Startup from initial conditions
  integer, public, parameter :: nsrStartup  = 0

  ! Continue from restart files
  integer, public, parameter :: nsrContinue = 1

  ! Branch from restart files
  integer, public, parameter :: nsrBranch   = 2

  ! true => allow case name to remain the same for branch run
  ! by default this is not allowed
  logical, public :: brnch_retain_casename = .false.

  !true => no valid land points -- do NOT run
  logical, public :: noland = .false.

  ! Hostname of machine running on
  character(len=256), public :: hostname = ' '

  ! username of user running program
  character(len=256), public :: username = ' '

  ! description of this source
  character(len=256), public :: source   = "E3SM Land Model"

  ! version of program
  character(len=256), public :: version  = " "

  ! dataset conventions
  character(len=256), public :: conventions = "CF-1.7"

  !----------------------------------------------------------
  ! Unit Numbers
  !----------------------------------------------------------
  !
  integer, public :: iulog = 6        ! "stdout" log file unit number, default is 6

  !----------------------------------------------------------
  ! Output NetCDF files
  !----------------------------------------------------------

  logical, public :: outnc_large_files = .true.         ! large file support for output NetCDF files

  !----------------------------------------------------------
  ! Run input files
  !----------------------------------------------------------

  character(len=fname_len), public :: finidat    = ' '        ! initial conditions file name
  character(len=fname_len), public :: fsurdat    = ' '        ! surface data file name
  character(len=fname_len), public :: fatmgrid   = ' '        ! atm grid file name
  character(len=fname_len), public :: fatmlndfrc = ' '        ! lnd frac file on atm grid
  character(len=fname_len), public :: fatmtopo   = ' '        ! topography on atm grid
  character(len=fname_len), public :: flndtopo   = ' '        ! topography on lnd grid
  character(len=fname_len), public :: paramfile  = ' '        ! ASCII data file with PFT physiological constants
  character(len=fname_len), public :: nrevsn     = ' '        ! restart data file name for branch run
  character(len=fname_len), public :: fsnowoptics  = ' '      ! snow optical properties file name
  character(len=fname_len), public :: fsnowaging   = ' '      ! snow aging parameters file name
  character(len=fname_len), public :: fsoilordercon    = ' '  ! ASCII data file with soil order dependent  constants

  !----------------------------------------------------------
  ! Flag to turn on MEGAN VOC's
  !----------------------------------------------------------

  logical, public :: use_voc = .true.

  !----------------------------------------------------------
  ! Interpolation of finidat if requested
  !----------------------------------------------------------

  logical, public :: bound_h2osoi = .true. ! for debugging

  ! If finidat_interp_source is non-blank and finidat is blank then interpolation will be done from
  ! finidat_interp_source to finidat_interp_dest

  character(len=fname_len), public :: finidat_interp_source = ' '
  character(len=fname_len), public :: finidat_interp_dest   = 'finidat_interp_dest.nc'

  !----------------------------------------------------------
  ! Irrigate logic
  !----------------------------------------------------------

  ! do not irrigate by default
  logical, public :: irrigate = .false.

  !----------------------------------------------------------
  ! Two-way coupled irrigation with MOSART
  !----------------------------------------------------------

  ! True is 2way, false is 1way
  logical, public :: tw_irr = .false.

  !----------------------------------------------------------
  ! Extra groundwater pumping for irrigation
  !----------------------------------------------------------

  ! True is extra pumping, false is stick with the gw fraction
  logical, public :: extra_gw_irr = .false.

  !----------------------------------------------------------
  ! FIRRIG data
  !----------------------------------------------------------

  ! True is read from surface data, false is constant
  logical, public :: firrig_data = .false.

  !----------------------------------------------------------
  ! Landunit logic
  !----------------------------------------------------------

  ! true => separate crop landunit is not created by default
  logical, public :: create_crop_landunit = .false.

  !----------------------------------------------------------
  ! Other subgrid logic
  !----------------------------------------------------------

  ! true => make ALL patches, cols & landunits active (even if weight is 0)
  logical, public :: all_active = .false.
  !$acc declare copyin(all_active)
  !----------------------------------------------------------
  ! BGC logic and datasets
  !----------------------------------------------------------

  ! values of 'prognostic','diagnostic','constant'
  character(len=16), public :: co2_type = 'constant'

  ! State of the model for the accelerated decomposition (AD) spinup.
  ! 0 (default) = normal model; 1 = AD SPINUP
  integer, public  :: spinup_state            = 0
  integer, public  :: nyears_ad_carbon_only   = 0
  real(r8), public :: spinup_mortality_factor = 1._r8
  !$acc declare copyin(spinup_state           )
  !$acc declare copyin(nyears_ad_carbon_only  )
  !$acc declare copyin(spinup_mortality_factor)

  ! true => anoxia is applied to heterotrophic respiration also considered in CH4 model
  ! default value reset in controlMod
  logical, public :: anoxia  = .true.
  !$acc declare copyin(anoxia)
  ! used to override an error check on reading in restart files
  logical, public :: override_bgc_restart_mismatch_dump = .false.

  ! Set in AllocationInit (TO DO - had to move it here to avoid circular dependency)
  logical, public  :: carbon_only
  logical, public  :: carbonnitrogen_only
  logical, public  :: carbonphosphorus_only
  !$acc declare create(carbon_only          )
  !$acc declare create(carbonnitrogen_only  )
  !$acc declare create(carbonphosphorus_only)
  !----------------------------------------------------------
  ! Physics
  !----------------------------------------------------------

  ! use subgrid fluxes
  integer,  public :: subgridflag = 1

  ! true => write global average diagnostics to std out
  logical,  public :: wrtdia       = .false.

  ! atmospheric CO2 molar ratio (by volume) (umol/mol)
  real(r8), public :: co2_ppmv     = 355._r8
  !$acc declare copyin(co2_ppmv)

  ! Use constant climate during transient run (CPL_BYPASS only)
  logical, public :: const_climate_hist  = .false.

  !$acc declare copyin(const_climate_hist)
  !----------------------------------------------------------
  ! C isotopes
  !----------------------------------------------------------

  logical, public :: use_c13 = .false.                  ! true => use C-13 model
  logical, public :: use_c14 = .false.                  ! true => use C-14 model

  !----------------------------------------------------------
  !  FATES switches
  !----------------------------------------------------------

  logical, public            :: use_fates = .false.                     ! true => use  ED
  integer, public            :: fates_spitfire_mode = 0                 ! 0 for no fire; 1 for constant ignitions
  character(len=256), public :: fates_harvest_mode = ''                 ! five different harvest modes; see namelist_definitions
  logical, public            :: use_fates_fixed_biogeog = .false.       ! true => use fixed biogeography mode
  logical, public            :: use_fates_planthydro = .false.          ! true => turn on fates hydro
  logical, public            :: use_fates_cohort_age_tracking = .false. ! true => turn on cohort age tracking
  logical, public            :: use_fates_tree_damage = .false.         ! true => turn on tree damage module
  logical, public            :: use_fates_ed_st3   = .false.            ! true => static stand structure
  logical, public            :: use_fates_ed_prescribed_phys = .false.  ! true => prescribed physiology
  logical, public            :: use_fates_inventory_init = .false.      ! true => initialize fates from inventory
  logical, public            :: use_fates_nocomp = .false.              ! true => no competition mode
  logical, public            :: use_fates_sp = .false.                  ! true => FATES satellite phenology mode
  logical, public            :: use_fates_luh = .false.                 ! true => FATES land use transitions mode
  logical, public            :: use_fates_lupft = .false.               ! true => FATES land use x pft mode
  logical, public            :: use_fates_potentialveg = .false.        ! true => FATES potential veg only
  character(len=256), public :: fluh_timeseries = ''                    ! filename for land use harmonization data
  character(len=256), public :: flandusepftdat = ''                     ! filename for fates landuse x pft data
  character(len=256), public :: fates_inventory_ctrl_filename = ''      ! filename for inventory control
  integer, public            :: fates_parteh_mode = -9                  ! 1 => carbon only
                                                                        ! 2 => C+N+P (not enabled yet)
                                                                        ! no others enabled
  integer, public            :: fates_seeddisp_cadence = iundef         ! 0 => no seed dispersal across gridcells
                                                                        ! 1, 2, 3  => daily, monthly, or yearly seed dispersal

  ! FATES history dimension level
  ! fates can produce history at either the daily timescale (dynamics)
  ! and the model step timescale. It can also generate output on the extra dimension
  ! Performing this output can be expensive, so we allow different history dimension
  ! levels.
  ! The first index is output at the model timescale
  ! The second index is output at the dynamics (daily) timescale      
  ! 0 - no output
  ! 1 - include only column level means (3D)
  ! 2 - include output that includes the 4th dimension

  integer, dimension(2), public   :: fates_history_dimlevel = (/2,2/)

  
  !----------------------------------------------------------
  !  BeTR switches
  !----------------------------------------------------------
  logical, public :: use_betr = .false.          ! true=> use BeTR
  !$acc declare create(use_betr)

  !----------------------------------------------------------
  ! lai streams switch for Sat. Phenology
  !----------------------------------------------------------

  logical, public :: use_lai_streams = .false. ! true => use lai streams in SatellitePhenologyMod.F90
  !----------------------------------------------------------
  ! plant hydraulic stress switch
  !----------------------------------------------------------

  logical, public :: use_hydrstress = .false. ! true => use plant hydraulic stress calculation

  !$acc declare create(use_hydrstress)
  !$acc declare create(use_lai_streams)
  !----------------------------------------------------------
  ! dynamic root switch
  !----------------------------------------------------------

  logical, public :: use_dynroot = .false. ! true => use dynamic root module
  !$acc declare create(use_dynroot)
  !----------------------------------------------------------
  ! glacier_mec control variables: default values (may be overwritten by namelist)
  ! NOTE: glc_smb must have the same values for CLM and GLC
  !----------------------------------------------------------

  ! glacier_mec landunit is not created (set in controlMod)
  logical , public :: create_glacier_mec_landunit = .false.
  !$acc declare create(create_glacier_mec_landunit)
  ! if true, pass surface mass balance info to GLC
  logical , public :: glc_smb = .true.
  !$acc declare create(glc_smb)
  ! if false, pass positive-degree-day info to GLC

  ! true => CLM glacier area & topography changes dynamically
  logical , public :: glc_do_dynglacier = .false.
  !$acc declare copyin(glc_do_dynglacier)
  ! true => downscale precip division into rain & snow
  logical , public :: glcmec_downscale_rain_snow_convert = .false.

  ! true => downscale longwave radiation
  logical , public :: glcmec_downscale_longwave = .true.

  ! number of days before one considers the perennially snow-covered point 'land ice'
  integer , public :: glc_snow_persistence_max_days = 7300
  !$acc declare copyin(glc_snow_persistence_max_days)

  ! glc_grid used to determine fglcmask
  character(len=256), public :: glc_grid = ' '

  ! glacier mask file name (based on glc_grid)
  character(len=fname_len), public :: fglcmask = ' '
  !
  !----------------------------------------------------------
  ! single column control variables
  !----------------------------------------------------------

  logical,  public :: single_column = .false. ! true => single column mode
  logical,  public :: scm_multcols  = .false. ! true => SCM domain mode
  real(r8), public :: scmlat        = rundef  ! single column lat
  real(r8), public :: scmlon        = rundef  ! single column lon
  integer,  public :: scm_nx        = -1      ! doubly periodic points x direction
  integer,  public :: scm_ny        = -1      ! doubly periodic points y direction

  !----------------------------------------------------------
  ! instance control
  !----------------------------------------------------------

  integer, public :: inst_index
  character(len=16), public :: inst_name
  character(len=16), public :: inst_suffix

  !----------------------------------------------------------
  ! Decomp control variables
  !----------------------------------------------------------

  ! number of segments per clump for decomp
  integer, public :: nsegspc = 20

  !----------------------------------------------------------
  ! Derived variables (run, history and restart file)
  !----------------------------------------------------------

  ! directory name for local restart pointer file
  character(len=256), public :: rpntdir = '.'

  ! file name for local restart pointer file
  character(len=256), public :: rpntfil = 'rpointer.lnd'

  ! moved hist_wrtch4diag from histFileMod.F90 to here - caused compiler error with intel
  ! namelist: write CH4 extra diagnostic output
  logical, public :: hist_wrtch4diag = .false.

  !----------------------------------------------------------
  ! ED/FATES
  !----------------------------------------------------------
  character(len=fname_len), public :: fates_paramfile  = ' '

  !----------------------------------------------------------
  ! Migration of CPP variables
  !----------------------------------------------------------

  logical, public :: use_nofire          = .false.
  logical, public :: use_lch4            = .false.
  logical, public :: use_vertsoilc       = .false.
  logical, public :: use_extralakelayers = .false.
  logical, public :: use_vichydro        = .false.
  logical, public :: use_century_decomp  = .false.
  logical, public :: use_cn              = .false.
  logical, public :: use_cndv            = .false.
  logical, public :: use_crop            = .false.
  logical, public :: use_snicar_frc      = .false.
  logical, public :: use_snicar_ad       = .false.
  logical, public :: use_extrasnowlayers = .false.
  logical, public :: use_firn_percolation_and_compaction  = .false.
  logical, public :: use_vancouver       = .false.
  logical, public :: use_mexicocity      = .false.
  logical, public :: use_noio            = .false.
  logical, public :: use_var_soil_thick  = .false.
  logical, public :: use_T_rho_dependent_snowthk     = .false.
  logical, public :: use_atm_downscaling_to_topunit  = .false.
  character(len = SHR_KIND_CS), public :: precip_downscaling_method  = 'ERMM' ! Precip downscaling method values can be ERMM or FNM
  logical, public :: use_lake_wat_storage = .false.
  logical, public :: use_top_solar_rad   = .false.  ! TOP : sub-grid topographic effect on surface solar radiation

  !----------------------------------------------------------
  ! Fan controls (use_fan)
  !----------------------------------------------------------
  logical, public :: use_fan             = .false.
  character(len=32), public :: fan_mode  = 'none'
  logical, public :: fan_nh3_to_atm      = .false.
  logical, public :: fan_to_bgc_crop     = .false.
  logical, public :: fan_to_bgc_veg      = .false.
 

  !----------------------------------------------------------
  ! VSFM switches
  !----------------------------------------------------------
  logical          , public :: use_vsfm                    = .false.
  logical          , public :: vsfm_use_dynamic_linesearch = .false.
  logical          , public :: vsfm_include_seepage_bc     = .false.
  character(len=32), public :: vsfm_satfunc_type           = 'smooth_brooks_corey_bz3'
  character(len=32), public :: vsfm_lateral_model_type     = 'none'

  !----------------------------------------------------------
  ! PETSc-based thermal model switches
  !----------------------------------------------------------
  logical, public :: use_petsc_thermal_model = .false.

  !----------------------------------------------------------
  ! Stub EM switches
  !----------------------------------------------------------
  logical          , public :: use_em_stub = .false.

  !----------------------------------------------------------
  ! To retrieve namelist
  !----------------------------------------------------------
  character(len=SHR_KIND_CL), public :: NLFilename_in ! Namelist filename
  !
  logical, private :: elmvarctl_isset = .false.
  !-----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  ! nutrient competition (nu_com), default is relative demand approach (RD)
  character(len=15), public :: nu_com = 'RD'
  !$acc declare copyin(nu_com)

  !-----------------------------------------------------------------------
  ! forest N/P fertilization
  logical, public :: forest_fert_exp = .false.
  !$acc declare copyin(forest_fert_exp)
  !-----------------------------------------------------------------------
  ! ECA regular spinup with P on, keep labile, secondary, occluded, parent
  ! material P being constant or not
  logical, public :: ECA_Pconst_RGspin = .false.
  !$acc declare copyin(ECA_Pconst_RGspin)
  !-----------------------------------------------------------------------
  ! Priority of plant to get symbiotic N fixation, phosphatase
  logical, public :: NFIX_PTASE_plant = .false.
  !$acc declare create(NFIX_PTASE_plant)
  !-----------------------------------------------------------------------
  !CO2 and warming experiments
  character(len=8), public :: startdate_add_temperature ='99991231'
  character(len=8), public :: startdate_add_co2         ='99991231'
  real(r8), public         :: add_co2 = 0d0
  real(r8), public         :: add_temperature = 0d0

  !-----------------------------------------------------------------------
  ! Lateral grid connectivity
  !-----------------------------------------------------------------------
  logical, public            :: lateral_connectivity  = .false.
  character(len=256), public :: domain_decomp_type    = 'round_robin'

  !-----------------------------------------------------------------------
  ! Subgrid hillslope hydrologic connectivity (through topounits)
  !-----------------------------------------------------------------------
  logical, public            :: use_IM2_hillslope_hydrology  = .false.
 
  !-----------------------------------------------------------------------
  ! flux limiter for phenology flux calculation
  logical, public :: use_pheno_flux_limiter = .false.

  ! Soil erosion
  !-----------------------------------------------------------------------
  logical, public :: use_erosion    = .false.   ! switch for turning on the soil erosion model
  logical, public :: ero_ccycle     = .false.   ! switch for turning on soil C, N and P loss by erosion (only valid when user_erosion = .true.)

  !$acc declare copyin(use_pheno_flux_limiter)
  !$acc declare copyin(use_erosion)
  !$acc declare copyin(ero_ccycle )
  !-----------------------------------------------------------------------
  ! bgc & pflotran interface
  !
  logical, public :: use_elm_interface  = .false.
  logical, public :: use_elm_bgc        = .false.
  logical, public :: use_pflotran       = .false.
  logical, public :: pf_surfaceflow     = .false.

  !$acc declare copyin(use_elm_interface)
  !$acc declare copyin(use_elm_bgc      )
  !$acc declare copyin(use_pflotran     )
  !$acc declare copyin(pf_surfaceflow   )

  ! the following switches will allow flexibility of coupling CLM with PFLOTRAN (which in fact runs in 3 modes individually or coupled)
  logical, public :: pf_cmode      = .false.                 ! switch for 'C' mode coupling (will be updated in interface)
  logical, public :: pf_hmode      = .false.                 ! switch for 'H' mode coupling (will be updated in interface)
  logical, public :: pf_tmode      = .false.                 ! switch for 'T' mode coupling (will be updated in interface)
  logical, public :: pf_frzmode    = .false.                 ! switch for 'freezing' mode availablity in PF-thmode (will be updated in interface)
  logical, public :: initth_pf2clm = .false.                 ! switch for initializing CLM TH states from pflotran
  integer, public :: pf_elmnstep0  = 0                       ! the ELM timestep of start/restart
  !$acc declare copyin(pf_cmode     )
  !$acc declare copyin(pf_hmode     )
  !$acc declare copyin(pf_tmode     )
  !$acc declare copyin(pf_frzmode   )
  !$acc declare copyin(initth_pf2clm)
  !$acc declare copyin(pf_clmnstep0 )

  ! cpl_bypass
   character(len=fname_len), public :: metdata_type   = ' '    ! metdata type for CPL_BYPASS mode
   character(len=fname_len), public :: metdata_bypass = ' '    ! met data directory for CPL_BYPASS mode (site, qian, cru_ncep)
   character(len=fname_len), public :: metdata_biases = ' '    ! met biases files for CPL_BYPASS mode
   character(len=fname_len), public :: co2_file       = ' '    ! co2 file for CPL_BYPASS mode
   character(len=fname_len), public :: aero_file      = ' '    ! aerosol deposition file for CPL_BYPASS mode

  !$acc declare create(use_fates)

  !$acc declare copyin(use_c13, use_cn, use_lch4, glcmec_downscale_rain_snow_convert)
  !$acc declare copyin(use_c14)
  !$acc declare copyin(glcmec_downscale_longwave, subgridflag)
  !$acc declare copyin(use_nofire         )
  !$acc declare copyin(use_lch4           )
  !$acc declare copyin(use_nitrif_denitrif)
  !$acc declare copyin(use_vertsoilc      )
  !$acc declare copyin(use_extralakelayers)
  !$acc declare copyin(use_vichydro       )
  !$acc declare copyin(use_century_decomp )
  !$acc declare copyin(use_cn             )
  !$acc declare copyin(use_crop           )
  !$acc declare copyin(use_snicar_frc     )
  !$acc declare copyin(use_snicar_ad      )
  !$acc declare copyin(use_vancouver      )
  !$acc declare copyin(use_mexicocity     )
  !$acc declare copyin(use_noio           )
  !$acc declare copyin(use_var_soil_thick )
  !$acc declare copyin(tw_irr)
  !$acc declare copyin(use_vsfm                   )
  !$acc declare copyin(vsfm_use_dynamic_linesearch)
  !$acc declare copyin(vsfm_include_seepage_bc    )
  !$acc declare copyin(vsfm_satfunc_type          )
  !$acc declare copyin(vsfm_lateral_model_type    )


  !----------------------------------------------------------
  ! Budgets
  !----------------------------------------------------------
   logical, public :: do_budgets   = .false.
   integer, public :: budget_inst  = 0
   integer, public :: budget_daily = 0
   integer, public :: budget_month = 1
   integer, public :: budget_ann   = 1
   integer, public :: budget_ltann = 1
   integer, public :: budget_ltend = 0

   !----------------------------------------------------------
   ! land river two way coupling
   !----------------------------------------------------------
   logical, public :: use_lnd_rof_two_way = .false.
   integer, public :: lnd_rof_coupling_nstep = 0
   
   
   !----------------------------------------------------------
   ! SNICAR-AD
   !----------------------------------------------------------
   character(len=256), public :: snow_shape = 'sphere'
   character(len=256), public :: snicar_atm_type = 'default'
   logical, public :: use_dust_snow_internal_mixing = .false.

   !----------------------------------------------------------
   ! MPI syncing
   !----------------------------------------------------------
   integer, public :: mpi_sync_nstep_freq = 0
   
   !----------------------------------------------------------
   ! Modified infiltration scheme in surface water storage
   !----------------------------------------------------------
   logical, public :: use_modified_infil = .false.

contains

  !---------------------------------------------------------------------------
  subroutine elm_varctl_set( caseid_in, ctitle_in, brnch_retain_casename_in,    &
       single_column_in, scm_multcols_in, scmlat_in, scmlon_in, scm_nx_in, scm_ny_in, &
       nsrest_in, version_in, hostname_in, username_in)
    !
    ! !DESCRIPTION:
    ! Set input control variables.
    !
    ! !ARGUMENTS:
    character(len=256), optional, intent(IN) :: caseid_in                ! case id
    character(len=256), optional, intent(IN) :: ctitle_in                ! case title
    logical,            optional, intent(IN) :: brnch_retain_casename_in ! true => allow case name to remain the
                                                                         ! same for branch run
    logical,            optional, intent(IN) :: single_column_in         ! true => single column mode
    logical,            optional, intent(IN) :: scm_multcols_in          ! SCM for multiple columns
    real(r8),           optional, intent(IN) :: scmlat_in                ! single column lat
    real(r8),           optional, intent(IN) :: scmlon_in                ! single column lon
    integer,            optional, intent(IN) :: scm_nx_in                ! x direction points
    integer,            optional, intent(IN) :: scm_ny_in                ! y direction points
    integer,            optional, intent(IN) :: nsrest_in                ! 0: initial run. 1: restart: 3: branch
    character(len=256), optional, intent(IN) :: version_in               ! model version
    character(len=256), optional, intent(IN) :: hostname_in              ! hostname running on
    character(len=256), optional, intent(IN) :: username_in              ! username running job
    !-----------------------------------------------------------------------

    if ( elmvarctl_isset )then
       call shr_sys_abort(' ERROR:: control variables already set, cannot call this routine')
    end if

    if ( present(caseid_in       ) ) caseid        = caseid_in
    if ( present(ctitle_in       ) ) ctitle        = ctitle_in
    if ( present(single_column_in) ) single_column = single_column_in
    if ( present(scm_multcols_in ) ) scm_multcols  = scm_multcols_in
    if ( present(scmlat_in       ) ) scmlat        = scmlat_in
    if ( present(scmlon_in       ) ) scmlon        = scmlon_in
    if ( present(scm_nx_in       ) ) scm_nx        = scm_nx_in
    if ( present(scm_ny_in       ) ) scm_ny        = scm_ny_in
    if ( present(nsrest_in       ) ) nsrest        = nsrest_in
    if ( present(brnch_retain_casename_in) ) brnch_retain_casename = brnch_retain_casename_in
    if ( present(version_in      ) ) version       = version_in
    if ( present(username_in     ) ) username      = username_in
    if ( present(hostname_in     ) ) hostname      = hostname_in

  end subroutine elm_varctl_set

  ! Set module carbon_only flag
  subroutine cnallocate_carbon_only_set(carbon_only_in)
    logical, intent(in) :: carbon_only_in
    carbon_only = carbon_only_in
  end subroutine cnallocate_carbon_only_set

  ! Get module carbon_only flag
  logical function CNAllocate_Carbon_only()
    cnallocate_carbon_only = carbon_only
  end function CNAllocate_Carbon_only

  ! Set module carbonnitrogen_only flag
  subroutine cnallocate_carbonnitrogen_only_set(carbonnitrogen_only_in)
    logical, intent(in) :: carbonnitrogen_only_in
    carbonnitrogen_only = carbonnitrogen_only_in
  end subroutine cnallocate_carbonnitrogen_only_set

  ! Get module carbonnitrogen_only flag
  logical function CNAllocate_CarbonNitrogen_only()
    cnallocate_carbonnitrogen_only = carbonnitrogen_only
  end function CNAllocate_CarbonNitrogen_only


  ! Set module carbonphosphorus_only flag
  subroutine cnallocate_carbonphosphorus_only_set(carbonphosphorus_only_in)
    logical, intent(in) :: carbonphosphorus_only_in
    carbonphosphorus_only = carbonphosphorus_only_in
  end subroutine cnallocate_carbonphosphorus_only_set

  ! Get module carbonphosphorus_only flag
  logical function CNAllocate_CarbonPhosphorus_only()
    cnallocate_carbonphosphorus_only = carbonphosphorus_only
  end function CNAllocate_CarbonPhosphorus_only

  function get_carbontag(carbon_type)result(ctag)
    !$acc routine seq
    implicit none
    character(len=*) :: carbon_type

    character(len=3) :: ctag

    if(carbon_type=='c12')then
       ctag = 'C'
    elseif(carbon_type=='c13')then
       ctag = 'C13'
    elseif(carbon_type=='c14')then
       ctag = 'C14'
    endif
  end function get_carbontag

end module elm_varctl
