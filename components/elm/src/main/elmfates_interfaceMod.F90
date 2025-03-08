module ELMFatesInterfaceMod

   ! -------------------------------------------------------------------------------------
   ! This module contains various functions and definitions to aid in the
   ! coupling of the FATES library/API with the CLM/ALM/ATS/etc model driver.
   ! All connections between the two models should occur in this file alone.
   !
   ! This is also the only location where CLM code is allowed to see FATES memory
   ! structures.
   ! The routines here, that call FATES library routines, cannot pass most types defined
   ! by the driving land model (HLM), only native type arrays (int,real,log, etc), implementations
   ! of fates abstract classes, and references into fates boundary condition structures.
   !
   ! Note that CLM/ALM does use Shared Memory Parallelism (SMP), where processes such as
   ! the update of state variables are forked.  However, IO is not assumed to be
   ! threadsafe and therefore memory spaces reserved for IO must be continuous vectors,
   ! and moreover they must be pushed/pulled from history IO for each individual
   ! bounds_proc memory space as a unit.
   !
   ! Therefore, the state variables in the alm_fates communicator is vectorized by
   ! threadcount, and the IO communication arrays are not.
   !
   !
   ! Conventions:
   ! keep line widths within 90 spaces
   ! HLM acronym = Host Land Model
   !
   ! -------------------------------------------------------------------------------------

   !  use ed_driver_interface, only:

   ! Used CLM Modules
   use VegetationType    , only : veg_pp
   use shr_kind_mod      , only : r8 => shr_kind_r8
   use decompMod         , only : bounds_type
   use WaterStateType    , only : waterstate_type
   use WaterFluxType     , only : waterflux_type
   use CanopyStateType   , only : canopystate_type
   use TemperatureType   , only : temperature_type
   use EnergyFluxType    , only : energyflux_type

   use SoilStateType     , only : soilstate_type
   use elm_varctl        , only : iulog
   use elm_varctl        , only : use_fates
   use elm_varctl        , only : use_vertsoilc
   use elm_varctl        , only : fates_spitfire_mode
   use elm_varctl        , only : fates_harvest_mode
   use elm_varctl        , only : fates_parteh_mode
   use elm_varctl        , only : fates_seeddisp_cadence
   use elm_varctl        , only : use_fates_planthydro
   use elm_varctl        , only : use_fates_cohort_age_tracking
   use elm_varctl        , only : use_fates_ed_st3
   use elm_varctl        , only : use_fates_ed_prescribed_phys
   use elm_varctl        , only : use_fates_inventory_init
   use elm_varctl        , only : use_fates_fixed_biogeog
   use elm_varctl        , only : use_fates_nocomp
   use elm_varctl        , only : use_fates_sp
   use elm_varctl        , only : use_fates_luh
   use elm_varctl        , only : use_fates_lupft
   use elm_varctl        , only : use_fates_potentialveg
   use elm_varctl        , only : flandusepftdat
   use elm_varctl        , only : use_fates_tree_damage
   use elm_varctl        , only : nsrest, nsrBranch
   use elm_varctl        , only : fates_inventory_ctrl_filename
   use elm_varctl        , only : fates_history_dimlevel
   use elm_varctl        , only : use_lch4
   use elm_varctl        , only : use_century_decomp
   use elm_varcon        , only : tfrz
   use elm_varcon        , only : spval
   use elm_varcon        , only : denice
   use elm_varcon        , only : ispval
   use elm_varctl        , only : nu_com
   use elm_varpar        , only : numrad
   use elm_varpar        , only : ivis
   use elm_varpar        , only : inir
   use elm_varpar        , only : nlevsoi
   use elm_varpar        , only : nlevdecomp
   use elm_varpar        , only : nlevdecomp_full
   use elm_varpar        , only : i_met_lit, i_cel_lit, i_lig_lit
   use elm_varpar        , only : surfpft_lb, surfpft_ub
   use elm_varpar        , only : natpft_size, max_patch_per_col, maxpatch_urb
   use elm_varpar        , only : numcft, maxpatch_urb
   use PhotosynthesisType , only : photosyns_type
   Use TopounitDataType  , only : topounit_atmospheric_flux, topounit_atmospheric_state
   use atm2lndType       , only : atm2lnd_type
   use SurfaceAlbedoType , only : surfalb_type
   use SolarAbsorbedType , only : solarabs_type
   use FrictionVelocityType , only : frictionvel_type
   use elm_time_manager  , only : is_restart
   use ncdio_pio         , only : file_desc_t, ncd_int, ncd_double
   use restUtilMod,        only : restartvar
   use elm_time_manager  , only : get_days_per_year, &
                                  get_curr_date,     &
                                  get_ref_date,      &
                                  timemgr_datediff,  &
                                  is_beg_curr_day,   &
                                  is_end_curr_month, &
                                  get_step_size,     &
                                  get_nstep
   use perf_mod          , only : t_startf, t_stopf

   use spmdMod           , only : masterproc
   use decompMod         , only : get_proc_bounds,   &
                                  get_proc_clumps,   &
                                  get_proc_global,   &
                                  get_clump_bounds

   use GridcellType      , only : grc_pp
   use TopounitDataType  , only : top_as
   use ColumnType        , only : col_pp
   use ColumnDataType    , only : col_es, col_ws, col_wf, col_cs, col_cf
   use ColumnDataType    , only : col_nf, col_pf
   use VegetationDataType, only : veg_es, veg_wf, veg_ws
   use LandunitType      , only : lun_pp

   use landunit_varcon   , only : istsoil
   use abortutils        , only : endrun
   use shr_log_mod       , only : errMsg => shr_log_errMsg
   use elm_varcon        , only : dzsoi_decomp
   use FuncPedotransferMod, only: get_ipedof
   use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type

   ! Used FATES Modules
   use FatesConstantsMod     , only : ifalse
   use FatesConstantsMod     , only : fates_check_param_set
   use FatesInterfaceMod     , only : fates_interface_type
   use FatesInterfaceMod     , only : allocate_bcin
   use FatesInterfaceMod     , only : allocate_bcpconst
   use FatesInterfaceMod     , only : allocate_bcout
   use FatesInterfaceMod     , only : set_bcpconst
   use FatesInterfaceMod     , only : SetFatesTime
   use FatesInterfaceMod     , only : SetFatesGlobalElements1
   use FatesInterfaceMod     , only : SetFatesGlobalElements2
   use FatesInterfaceMod     , only : set_fates_ctrlparms
   use FatesInterfaceMod     , only : zero_bcs
   use FatesInterfaceMod     , only : FatesInterfaceInit
   use FatesInterfaceMod     , only : UpdateFatesRMeansTStep
   use FatesInterfaceMod     , only : InitTimeAveragingGlobals
   use FatesInterfaceTypesMod, only : fates_maxPatchesPerSite
   use FatesInterfaceMod     , only : DetermineGridCellNeighbors
   use FatesHistoryInterfaceMod, only : fates_hist
   use FatesRestartInterfaceMod, only : fates_restart_interface_type
   use FatesInterfaceTypesMod,   only : hlm_num_luh2_states
   use FatesIOVariableKindMod, only : group_dyna_simple, group_dyna_complx
   use PRTGenericMod         , only : num_elements
   use FatesPatchMod         , only : fates_patch_type
   use FatesDispersalMod     , only : lneighbors, dispersal_type, IsItDispersalTime
   use FatesInterfaceTypesMod, only : hlm_stepsize, hlm_current_day
   use EDMainMod             , only : ed_ecosystem_dynamics
   use EDMainMod             , only : ed_update_site
   use EDInitMod             , only : zero_site
   use EDInitMod             , only : init_site_vars
   use EDInitMod             , only : init_patches
   use EDInitMod             , only : set_site_properties
   use EDPftVarcon           , only : EDpftvarcon_inst
   use FatesRadiationDriveMod, only : FatesSunShadeFracs
   use FatesRadiationDriveMod, only : FatesNormalizedCanopyRadiation
   use EDBtranMod            , only : btran_ed, &
                                      get_active_suction_layers
   use EDCanopyStructureMod  , only : canopy_summarization, update_hlm_dynamics
   use EDCanopyStructureMod  , only : UpdateFatesAvgSnowDepth
   use FatesPlantRespPhotosynthMod, only : FatesPlantRespPhotosynthDrive
   use EDAccumulateFluxesMod , only : AccumulateFluxes_ED
   use FatesSoilBGCFluxMod   , only : UnPackNutrientAquisitionBCs
   use FatesSoilBGCFluxMod   , only : FluxIntoLitterPools
   use PRTGenericMod, only : prt_cnp_flex_allom_hyp
   use FatesPlantHydraulicsMod, only : hydraulics_drive
   use FatesPlantHydraulicsMod, only : HydrSiteColdStart
   use FatesPlantHydraulicsMod, only : InitHydrSites
   use FatesPlantHydraulicsMod, only : RestartHydrStates
   
   use FATESFireBase          , only : fates_fire_base_type
   use FATESFireFactoryMod    , only : no_fire, scalar_lightning, successful_ignitions,&
                                       anthro_ignitions, anthro_suppression

   use dynHarvestMod          , only : num_harvest_vars, harvest_varnames, wood_harvest_units
   use dynHarvestMod          , only : harvest_rates  ! these are dynamic in space and time
   use FatesConstantsMod      , only : hlm_harvest_area_fraction
   use FatesConstantsMod      , only : hlm_harvest_carbon

   use dynFATESLandUseChangeMod, only : num_landuse_transition_vars
   use dynFATESLandUseChangeMod, only : num_landuse_state_vars
   use dynFATESLandUseChangeMod, only : num_landuse_harvest_vars
   use dynFATESLandUseChangeMod, only : landuse_transitions
   use dynFATESLandUseChangeMod, only : landuse_states
   use dynFATESLandUseChangeMod, only : landuse_harvest
   use dynFATESLandUseChangeMod, only : landuse_transition_varnames
   use dynFATESLandUseChangeMod, only : landuse_state_varnames
   use dynFATESLandUseChangeMod, only : landuse_harvest_varnames
   use dynFATESLandUseChangeMod, only : landuse_harvest_units
   use dynFATESLandUseChangeMod, only : fates_harvest_no_logging
   use dynFATESLandUseChangeMod, only : fates_harvest_hlmlanduse
   use dynFATESLandUseChangeMod, only : fates_harvest_luh_area
   use dynFATESLandUseChangeMod, only : fates_harvest_luh_mass

   use FatesInterfaceTypesMod       , only : bc_in_type, bc_out_type

   use ELMFatesParamInterfaceMod, only : fates_param_reader_ctsm_impl
   use FatesParametersInterface, only : fates_param_reader_type
   use FatesParametersInterface, only : fates_parameters_type
   
   use perf_mod          , only : t_startf, t_stopf

   use FatesInterfaceTypesMod, only : fates_dispersal_cadence_none
   
   use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)

   implicit none


   type, public :: f2hmap_type

      ! This is the associated column index of each FATES site
      integer, allocatable :: fcolumn (:)

      ! This is the associated site index of any HLM columns
      ! This vector may be sparse, and non-sites have index 0
      integer, allocatable :: hsites  (:)

   end type f2hmap_type


   type, public :: hlm_fates_interface_type

      ! See above for descriptions of the sub-types populated
      ! by thread.  This type is somewhat self-explanatory, in that it simply
      ! breaks up memory and process by thread.  Each thread will have its
      ! own list of sites, and boundary conditions for those sites

      type(fates_interface_type), allocatable :: fates (:)


      ! This memory structure is used to map fates sites
      ! into the host model.  Currently, the FATES site
      ! and its column number matching are its only members

      type(f2hmap_type), allocatable  :: f2hmap(:)

      ! fates_restart is the inteface calss for restarting the model
      type(fates_restart_interface_type) :: fates_restart
      
      ! fates_fire_data_method determines the fire data passed from HLM to FATES
      class(fates_fire_base_type), allocatable :: fates_fire_data_method

      ! Type structure that holds allocatable arrays for mpi-based seed dispersal
      type(dispersal_type) :: fates_seed
      
   contains

      procedure, public :: init
      procedure, public :: check_hlm_active
      procedure, public :: restart
      procedure, public :: init_coldstart
      procedure, public :: dynamics_driv
      procedure, public :: wrap_sunfrac
      procedure, public :: wrap_btran
      procedure, public :: wrap_photosynthesis
      procedure, public :: wrap_accumulatefluxes
      procedure, public :: prep_canopyfluxes
      procedure, public :: wrap_canopy_radiation
      procedure, public :: wrap_WoodProducts
      procedure, public :: wrap_update_hifrq_hist
      procedure, public :: TransferZ0mDisp
      procedure, public :: InterpFileInputs  ! Interpolate inputs from files
      procedure, public :: Init2  ! Initialization after determining subgrid weights
      procedure, public :: InitAccBuffer ! Initialize any accumulation buffers
      procedure, public :: InitAccVars   ! Initialize any accumulation variables
      procedure, public :: UpdateAccVars ! Update any accumulation variables
      procedure, public :: UpdateLitterFluxes
      procedure, private :: init_history_io
      procedure, private :: wrap_update_hlmfates_dyn
      procedure, private :: init_soil_depths
      procedure, public  :: ComputeRootSoilFlux
      procedure, public  :: wrap_hydraulics_drive
      procedure, public  :: WrapUpdateFatesRmean
      procedure, public  :: WrapGlobalSeedDispersal
      procedure, public  :: WrapUpdateFatesSeedInOut

   end type hlm_fates_interface_type
   
   ! hlm_bounds_to_fates_bounds is not currently called outside the interface.
   ! Although there may be good reasons to, I privatized it so that the next
   ! developer will at least question its usage (RGK)
   private :: hlm_bounds_to_fates_bounds

   ! The GetAndSetTime function is used to get the current time from the CLM
   ! time procedures and then set to the fates global time variables during restart,
   ! init_coldstart, and dynamics_driv function calls
   private :: GetAndSetTime


   logical :: debug  = .false.

   character(len=*), parameter, private :: sourcefile = &
        __FILE__

   integer, parameter :: num_landuse_pft_vars = 4

   public  :: ELMFatesGlobals1
   public  :: ELMFatesGlobals2
   public  :: ELMFatesTimesteps
   public  :: CrossRefHistoryFields
   
contains

  subroutine ELMFatesGlobals1()

    ! --------------------------------------------------------------------------------
    ! This is the first call to fates
    ! We open the fates parameter file. And use that and some info on
    ! namelist variables to determine how many patches need to be allocated
    ! in ELM
    ! --------------------------------------------------------------------------------
    integer                                        :: pass_biogeog
    integer                                        :: pass_nocomp
    integer                                        :: pass_sp
    integer                                        :: pass_masterproc
    logical                                        :: verbose_output
    type(fates_param_reader_ctsm_impl)             :: var_reader

    if (use_fates) then

       verbose_output = .false.
       call FatesInterfaceInit(iulog, verbose_output)

       ! Force FATES parameters that are recieve type, to the unset value
       call set_fates_ctrlparms('flush_to_unset')

       ! Send parameters individually
       if(use_fates_fixed_biogeog)then
          pass_biogeog = 1
       else
          pass_biogeog = 0
       end if
       call set_fates_ctrlparms('use_fixed_biogeog',ival=pass_biogeog)

       if(use_fates_nocomp)then
          pass_nocomp = 1
       else
          pass_nocomp = 0
       end if
       call set_fates_ctrlparms('use_nocomp',ival=pass_nocomp)

       if(use_fates_sp)then
          pass_sp = 1
       else
          pass_sp = 0
       end if
       call set_fates_ctrlparms('use_sp',ival=pass_sp)

       if(masterproc)then
          pass_masterproc = 1
       else
          pass_masterproc = 0
       end if
       call set_fates_ctrlparms('masterproc',ival=pass_masterproc)

    end if

    ! The following call reads in the parameter file
    ! and then uses that to determine the number of patches
    ! FATES requires. We pass that to CLM here
    ! so that it can perform some of its allocations.
    ! During init 2, we will perform more size checks
    ! and allocations on the FATES side, which require
    ! some allocations from CLM (like soil layering)

    ! The third argument (num crop pfts) should always be zero.
    ! If use_crop is off then the natpft_size includes the
    ! crop pfts tacked on. If use_crop is true, then we don't
    ! want fates to handle crops, so again, it should be ignored.
    ! (RGK 07-2022)
    
    call SetFatesGlobalElements1(use_fates,natpft_size,0,var_reader)

    natpft_size = fates_maxPatchesPerSite
    max_patch_per_col= max(natpft_size, numcft, maxpatch_urb)

    return
  end subroutine ELMFatesGlobals1
  
  ! ====================================================================================
   
  subroutine ELMFatesGlobals2()

     ! --------------------------------------------------------------------------------
     ! This is one of the first calls to fates
     ! Used for setting dimensions.  This MUST
     ! be called after NL variables are specified and
     ! after the FATES parameter file has been read in
     ! Aside from setting global dimension info, which
     ! is used in the history file, we also transfer
     ! over the NL variables to FATES global settings.
     ! --------------------------------------------------------------------------------

     integer                                        :: pass_vertsoilc
     integer                                        :: pass_ed_st3
     integer                                        :: pass_logging
     integer                                        :: pass_ch4
     integer                                        :: pass_ed_prescribed_phys
     integer                                        :: pass_planthydro
     integer                                        :: pass_inventory_init
     integer                                        :: pass_is_restart
     integer                                        :: pass_cohort_age_tracking
     integer                                        :: pass_num_lu_harvest_cats
     integer                                        :: pass_lu_harvest
     integer                                        :: pass_tree_damage
     integer                                        :: pass_use_luh
     integer                                        :: pass_use_potentialveg     
     integer                                        :: pass_num_luh_states
     integer                                        :: pass_num_luh_transitions
     integer                                        :: pass_lupftdat
     ! ----------------------------------------------------------------------------------
     ! FATES lightning definitions
     ! 1 : use a global constant lightning rate found in fates_params.
     ! 2 : use an external lightning dataset.
     ! 3 : use an external confirmed ignitions dataset.
     ! 4 : use external lightning and population datasets to simulate
     !     both natural and anthropogenic
     ! Special note* external_lightning is not passed as a definition, because
     ! this definition is not used in FATES (but will be added in later API)
     ! as it is the resulting logic if not 0,1,3 or 4.
     ! ----------------------------------------------------------------------------------

     integer, parameter :: no_fire = 0
     integer, parameter :: scalar_lightning = 1
     integer, parameter :: external_lightning = 2
     integer, parameter :: successful_ignitions = 3
     integer, parameter :: anthro_ignitions= 4
  
     ! We will use this switch temporarily, until  we complete
     ! the ELM-FATES harvest integration
     logical, parameter :: do_elm_fates_harvest = .false.

     if (use_fates) then

        ! Send parameters individually
        call set_fates_ctrlparms('num_sw_bbands',ival=numrad)
        call set_fates_ctrlparms('vis_sw_index',ival=ivis)
        call set_fates_ctrlparms('nir_sw_index',ival=inir)
        call set_fates_ctrlparms('num_lev_soil',ival=nlevsoi)
        call set_fates_ctrlparms('hlm_name',cval='ELM')
        call set_fates_ctrlparms('hio_ignore_val',rval=spval)
        call set_fates_ctrlparms('soilwater_ipedof',ival=get_ipedof(0))
        call set_fates_ctrlparms('parteh_mode',ival=fates_parteh_mode)
        call set_fates_ctrlparms('seeddisp_cadence',ival=fates_seeddisp_cadence)

        call set_fates_ctrlparms('hist_hifrq_dimlevel',ival=fates_history_dimlevel(1))
        call set_fates_ctrlparms('hist_dynam_dimlevel',ival=fates_history_dimlevel(2))
        
        if(use_fates_tree_damage)then
           pass_tree_damage = 1
        else
           pass_tree_damage = 0
        end if
        call set_fates_ctrlparms('use_tree_damage',ival=pass_tree_damage)
        
        if((trim(nu_com)=='ECA') .or. (trim(nu_com)=='MIC')) then
           call set_fates_ctrlparms('nu_com',cval='ECA')
        else
           call set_fates_ctrlparms('nu_com',cval='RD')
        end if
        
        ! Tell fates which decomposition method is being used by the host land model
        if (use_century_decomp) then 
           call set_fates_ctrlparms('decomp_method',cval='CENTURY')
        else
           call set_fates_ctrlparms('decomp_method',cval='CTC')
        end if

        ! ELM ALWAYS has nitrogen and phosphorus "on"
        ! These may be in a non-limiting status (ie when supplements)
        ! are added, but they are always allocated and cycled non-the less
        ! FATES may want to interact differently with other models
        ! that don't even have these arrays allocated.
        ! FATES also checks that if NO3 is cycled in ELM, then
        ! any plant affinity parameters are checked.

        call set_fates_ctrlparms('nitrogen_spec',ival=1)
        call set_fates_ctrlparms('phosphorus_spec',ival=1)
        
        if(is_restart() .or. nsrest == nsrBranch) then
           pass_is_restart = 1
        else
           pass_is_restart = 0
        end if
        call set_fates_ctrlparms('is_restart',ival=pass_is_restart)

        if(use_lch4) then
           pass_ch4 = 1
        else
           pass_ch4 = 0
        end if
        call set_fates_ctrlparms('use_ch4',ival=pass_ch4)

        if(use_vertsoilc) then
           pass_vertsoilc = 1
        else
           pass_vertsoilc = 0
        end if
        call set_fates_ctrlparms('use_vertsoilc',ival=pass_vertsoilc)


        ! A note about spitfire: The CTSM team has added dataset
        ! read capabilities for lightning and ignition sources.
        ! THerefore since spitfire can be run in these different
        ! ignition "modes", the "use_spitfire" on the FATES side
        ! is no longer a binary switch, and is a "mode" flag where
        ! 0 = no fire
        ! 1 = spitfire on, with a scalar lighting rate
        ! 2 = spitfire on, lightning strike rate from data
        ! 3 = spitfire on, successful ignitions from data
        ! 4 = spitfire on, anthro ignitions from data
        ! Until we add this data stream connectivity in E3SM-Fates
        ! we will always pass 0 or 1.
        ! But.. we do have to define these modes, for now.


        call set_fates_ctrlparms('spitfire_mode',ival=fates_spitfire_mode)
        call set_fates_ctrlparms('sf_nofire_def',ival=no_fire)
        call set_fates_ctrlparms('sf_scalar_lightning_def',ival=scalar_lightning)
        call set_fates_ctrlparms('sf_successful_ignitions_def',ival=successful_ignitions)
        call set_fates_ctrlparms('sf_anthro_ignitions_def',ival=anthro_ignitions)

        ! FATES logging and harvest modes
        pass_logging = 0
        pass_lu_harvest = 0
        pass_num_lu_harvest_cats = 0
        if (trim(fates_harvest_mode) /= fates_harvest_no_logging) then
           pass_logging = 1 ! Time driven logging, without landuse harvest
           ! CLM landuse timeseries driven harvest rates
           if (trim(fates_harvest_mode) == fates_harvest_hlmlanduse) then
              pass_num_lu_harvest_cats = num_harvest_vars
              pass_lu_harvest = 1

           ! LUH2 landuse timeseries driven  harvest rates
           else if (trim(fates_harvest_mode) == fates_harvest_luh_area .or. &
                    trim(fates_harvest_mode) == fates_harvest_luh_mass) then
              pass_lu_harvest = 1
              pass_num_lu_harvest_cats = num_landuse_harvest_vars
           end if
        end if

        call set_fates_ctrlparms('use_lu_harvest',ival=pass_lu_harvest)
        call set_fates_ctrlparms('num_lu_harvest_cats',ival=pass_num_lu_harvest_cats)
        call set_fates_ctrlparms('use_logging',ival=pass_logging)

        if(use_fates_luh) then
           pass_use_luh = 1
           pass_num_luh_states = num_landuse_state_vars
           pass_num_luh_transitions = num_landuse_transition_vars
        else
           pass_use_luh = 0
           pass_num_luh_states = 0
           pass_num_luh_transitions = 0
        end if

        call set_fates_ctrlparms('use_luh2',ival=pass_use_luh)
        call set_fates_ctrlparms('num_luh2_states',ival=pass_num_luh_states)
        call set_fates_ctrlparms('num_luh2_transitions',ival=pass_num_luh_transitions)

        if ( use_fates_potentialveg ) then
           pass_use_potentialveg = 1
        else
           pass_use_potentialveg = 0
        end if
        call set_fates_ctrlparms('use_fates_potentialveg',ival=pass_use_potentialveg)

        if(use_fates_ed_st3) then
           pass_ed_st3 = 1
        else
           pass_ed_st3 = 0
        end if
        call set_fates_ctrlparms('use_ed_st3',ival=pass_ed_st3)

        if(use_fates_ed_prescribed_phys) then
           pass_ed_prescribed_phys = 1
        else
           pass_ed_prescribed_phys = 0
        end if
        call set_fates_ctrlparms('use_ed_prescribed_phys',ival=pass_ed_prescribed_phys)

        if(use_fates_planthydro) then
           pass_planthydro = 1
        else
           pass_planthydro = 0
        end if
        call set_fates_ctrlparms('use_planthydro',ival=pass_planthydro)

        if(use_fates_cohort_age_tracking) then
           pass_cohort_age_tracking = 1
        else
           pass_cohort_age_tracking = 0
        end if
        call set_fates_ctrlparms('use_cohort_age_tracking',ival=pass_cohort_age_tracking)

        if(use_fates_inventory_init) then
           pass_inventory_init = 1
        else
           pass_inventory_init = 0
        end if
        call set_fates_ctrlparms('use_inventory_init',ival=pass_inventory_init)

        call set_fates_ctrlparms('inventory_ctrl_file',cval=fates_inventory_ctrl_filename)

        ! Check through FATES parameters to see if all have been set
        call set_fates_ctrlparms('check_allset')

     end if

     ! This determines the total amount of space it requires in its largest
     ! dimension.  We are currently calling that the "cohort" dimension, but
     ! it is really a utility dimension that captures the models largest
     ! size need.
     ! Sets:
     ! fates_maxElementsPerPatch
     ! num_elements
     ! fates_maxElementsPerSite (where a site is roughly equivalent to a column)
     ! (Note: this needs to be called when use_fates=.false. as well, becuase
     ! it will return some nominal dimension sizes of 1

     call SetFatesGlobalElements2(use_fates)
     
     return
   end subroutine ELMFatesGlobals2

   ! ===================================================================================

   subroutine CrossRefHistoryFields

     ! This routine only needs to be called on the masterproc.
     ! Here we cross reference the ELM history master
     ! list and make sure that all fields that start
     ! with fates have been allocated. If it has
     ! not, then we give a more constructive error
     ! message than what is possible in PIO. The user
     ! most likely needs to increase the history density
     ! level

     use histFileMod, only: getname
     use histFileMod, only: hist_fincl1,hist_fincl2,hist_fincl3,hist_fincl4
     use histFileMod, only: hist_fincl5,hist_fincl6
     use histFileMod, only: max_tapes, max_flds, max_namlen

     integer :: t     ! iterator index for history tapes
     integer :: f     ! iterator index for registered history field names
     integer :: nh    ! iterator index for fates registered history
     logical :: is_fates_field ! Does this start with FATES_ ?
     logical :: found ! if true, than the history field is either
                      ! not part of the fates set, or was found in
                      ! the fates set
     character(len=64) :: fincl_name
     ! This is a copy of the public in histFileMod, copied
     ! here because it isn't filled at the time of this call
     character(len=max_namlen+2) :: fincl(max_flds,max_tapes)

     fincl(:,1)  = hist_fincl1(:)
     fincl(:,2)  = hist_fincl2(:)
     fincl(:,3)  = hist_fincl3(:)
     fincl(:,4)  = hist_fincl4(:)
     fincl(:,5)  = hist_fincl5(:)
     fincl(:,6)  = hist_fincl6(:)

     do t = 1,max_tapes

        f = 1
        search_fields: do while (f < max_flds .and. fincl(f,t) /= ' ')

           fincl_name = getname(fincl(f,t))
           is_fates_field = fincl_name(1:6)=='FATES_'

           if(is_fates_field) then
              found = .false.
              do_fates_hist: do nh = 1,fates_hist%num_history_vars()
                 if(trim(fates_hist%hvars(nh)%vname) == &
                      trim(fincl_name)) then
                    found=.true.
                    exit do_fates_hist
                 end if
              end do do_fates_hist

              if(.not.found)then
                 write(iulog,*) 'the history field: ',trim(fincl_name)
                 write(iulog,*) 'was requested in the namelist, but was'
                 write(iulog,*) 'not found in the list of fates_hist%hvars.'
                 write(iulog,*) 'Most likely, this is because this history variable'
                 write(iulog,*) 'was specified in the user namelist, but the user'
                 write(iulog,*) 'specified a FATES history output dimension level'
                 write(iulog,*) 'that does not contain that variable in its valid set.'
                 write(iulog,*) 'You may have to increase the namelist setting: fates_history_dimlevel'
                 write(iulog,*) 'current fates_history_dimlevel: ',fates_history_dimlevel(:)
                 if (debug) then
                    !if you want to list all fates history variables in registry turn on debug
                    do_fates_hist2: do nh = 1,fates_hist%num_history_vars()
                       write(iulog,*) trim(fates_hist%hvars(nh)%vname)
                    end do do_fates_hist2
                 end if
                 call endrun(msg=errMsg(sourcefile, __LINE__))
              end if
           end if
           f = f + 1
        end do search_fields

     end do
   end subroutine CrossRefHistoryFields
   
   ! ====================================================================================
   
   subroutine ELMFatesTimesteps()
     
     hlm_stepsize = real(get_step_size(),r8)
     
     call InitTimeAveragingGlobals()
     
     return
   end subroutine ELMFatesTimesteps
   
   ! ====================================================================================

   subroutine init(this, bounds_proc, flandusepftdat)

      ! ---------------------------------------------------------------------------------
      ! This initializes the hlm_fates_interface_type
      !
      ! sites is the root of the ED state hierarchy (instantaneous info on
      ! the state of the ecosystem).  As such, it governs the connection points between
      ! the host (which also dictates its allocation) and its patch structures.
      !
      ! sites may associate with different scales in different models. In
      ! CLM, it is being designed to relate to column scale.
      !
      ! This global may become relegated to this module.
      !
      ! Note: CLM/ALM currently wants sites to be allocated even if ed
      ! is not turned on
      ! ---------------------------------------------------------------------------------

      use spmdMod,                  only : npes
      use decompMod,                only : procinfo
      use FatesInterfaceMod,        only : FatesReportParameters
      use FatesParameterDerivedMod, only : param_derived
      use FatesInterfaceTypesMod,   only : numpft_fates => numpft
      use elm_varsur,               only : wt_nat_patch
      use topounit_varcon,          only : max_topounits, has_topounit
      use FATESFireFactoryMod,      only : create_fates_fire_data_method

      implicit none

      ! Input Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                   :: bounds_proc
      character(len=*), intent(in)                   :: flandusepftdat

      ! local variables
      integer                                        :: nclumps   ! Number of threads
      integer                                        :: nc        ! thread index
      integer                                        :: s         ! FATES site index
      integer                                        :: t         ! topounit index (HLM)
      integer                                        :: c         ! HLM column index
      integer                                        :: l         ! HLM LU index
      integer                                        :: g         ! HLM grid index
      integer                                        :: m         ! HLM PFT index
      integer                                        :: ft        ! FATES PFT index
      integer                                        :: pi,pf
      integer, allocatable                           :: collist (:)
      type(bounds_type)                              :: bounds_clump
      integer                                        :: nmaxcol
      integer                                        :: ndecomp
      integer                                        :: numg

      real(r8), allocatable :: landuse_pft_map(:,:,:)
      real(r8), allocatable :: landuse_bareground(:)

      ! Initialize the FATES communicators with the HLM
      ! This involves to stages
      ! 1) allocate the vectors
      ! 2) add the history variables defined in clm_inst to the history machinery
      call param_derived%Init( numpft_fates )

      ! Initialize dispersal
      if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then
         ! Initialize fates global seed dispersal array for all nodes
         call get_proc_global(ng=numg)
         call this%fates_seed%init(npes,numg,procinfo%ncells,numpft_fates)

         ! Initialize the array of nearest neighbors for fates-driven grid cell communications
         ! This must be called after surfrd_get_data and decompInit_lnd
         call DetermineGridCellNeighbors(lneighbors,this%fates_seed,numg)
      end if

      nclumps = get_proc_clumps()
      allocate(this%fates(nclumps))
      allocate(this%f2hmap(nclumps))

      if(debug)then
         write(iulog,*) 'alm_fates%init():  allocating for ',nclumps,' threads'
      end if

      ! Retrieve the landuse x pft static data if the optional switch has been set
      if (use_fates_fixed_biogeog .and. use_fates_luh) then
         call GetLandusePFTData(bounds_proc, flandusepftdat, landuse_pft_map, landuse_bareground)
      end if

      nclumps = get_proc_clumps()

      !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,nmaxcol,s,c,l,g,collist,pi,pf,ft)
      do nc = 1,nclumps

         call get_clump_bounds(nc, bounds_clump)
         nmaxcol = bounds_clump%endc - bounds_clump%begc + 1

         allocate(collist(1:nmaxcol))

         ! Allocate the mapping that points columns to FATES sites, 0 is NA
         allocate(this%f2hmap(nc)%hsites(bounds_clump%begc:bounds_clump%endc))

         ! Initialize all columns with a zero index, which indicates no FATES site
         this%f2hmap(nc)%hsites(:) = 0

         s = 0
         do c = bounds_clump%begc,bounds_clump%endc
            l = col_pp%landunit(c)

            ! These are the key constraints that determine if this column
            ! will have a FATES site associated with it

            ! INTERF-TODO: WE HAVE NOT FILTERED OUT FATES SITES ON INACTIVE COLUMNS.. YET
            ! NEED A RUN-TIME ROUTINE THAT CLEARS AND REWRITES THE SITE LIST

            if ( (lun_pp%itype(l) == istsoil) .and. (col_pp%active(c)) ) then
               s = s + 1
               collist(s) = c
               this%f2hmap(nc)%hsites(c) = s
               col_pp%is_fates(c) = .true.

               if(debug)then
                  write(iulog,*) 'alm_fates%init(): thread',nc,': found column',c,'with lu',l
                  write(iulog,*) 'LU type:', lun_pp%itype(l)
               end if
            endif

         enddo

         if(debug)then
           write(iulog,*) 'alm_fates%init(): thread',nc,': allocated ',s,' sites'
         end if

         ! Allocate vectors that match FATES sites with HLM columns
         ! RGK: Sites and fcolumns are forced as args during clm_driv() as of 6/4/2016
         ! We may have to give these a dummy allocation of 1, which should
         ! not be a problem since we always iterate on nsites.

         allocate(this%f2hmap(nc)%fcolumn(s))

         ! Assign the h2hmap indexing
         this%f2hmap(nc)%fcolumn(1:s)         =  collist(1:s)

         ! Deallocate the temporary arrays
         deallocate(collist)

         ! Set the number of FATES sites
         this%fates(nc)%nsites = s

         ! Allocate the FATES sites
         allocate (this%fates(nc)%sites(this%fates(nc)%nsites))

         ! Allocate the FATES boundary arrays (in)
         allocate(this%fates(nc)%bc_in(this%fates(nc)%nsites))

         ! Allocate the FATES boundary arrays (out)
         allocate(this%fates(nc)%bc_out(this%fates(nc)%nsites))

         ! Allocate and Initialize the Boundary Condition Arrays
         ! These are staticaly allocated at maximums, so
         ! No information about the patch or cohort structure is needed at this step

         ! Parameter Constants defined by FATES, but used in ELM
         ! Note that FATES has its parameters defined, so we can also set the values
         call allocate_bcpconst(this%fates(nc)%bc_pconst,nlevdecomp)

         ! This also needs
         call set_bcpconst(this%fates(nc)%bc_pconst,nlevdecomp)

         do s = 1, this%fates(nc)%nsites

            c = this%f2hmap(nc)%fcolumn(s)
            this%fates(nc)%sites(s)%h_gid = c

            if (use_vertsoilc) then
               ndecomp = col_pp%nlevbed(c)
            else
               ndecomp = 1
            end if

            call allocate_bcin(this%fates(nc)%bc_in(s), col_pp%nlevbed(c), ndecomp, &
                               num_harvest_vars, num_landuse_state_vars, num_landuse_transition_vars, &
                               surfpft_lb, surfpft_ub)
            call allocate_bcout(this%fates(nc)%bc_out(s),col_pp%nlevbed(c),ndecomp)
            call zero_bcs(this%fates(nc),s)

            ! Pass any grid-cell derived attributes to the site
            ! ---------------------------------------------------------------------------

            g = col_pp%gridcell(c)

            this%fates(nc)%sites(s)%lat = grc_pp%latdeg(g)
            this%fates(nc)%sites(s)%lon = grc_pp%londeg(g)

            ! Check whether or not the surface dataset has topounits.  If it doesn't set the
            ! index t to max_topounits, which should be 1.  Otherwise, determine the index
            ! from the columntype
            if (has_topounit) then
               t = col_pp%topounit(c)
            else
               if (max_topounits .ne. 1) then
                  write(iulog,*) 'max_topounits should only be one when has_topounit is false'
                  write(iulog,*) 'max_topounits, has_topounit: ', max_topounits, has_topounit
                  call endrun(msg=errMsg(sourcefile, __LINE__))
               else
                  t = max_topounits
               endif
            endif

            if (use_fates_fixed_biogeog) then
               ! Transfer the landuse x pft data to fates via bc_in if landuse mode engaged
               if (use_fates_luh) then
                  this%fates(nc)%bc_in(s)%pft_areafrac_lu(:,1:num_landuse_pft_vars) = landuse_pft_map(g,:,1:num_landuse_pft_vars)
                  this%fates(nc)%bc_in(s)%baregroundfrac = landuse_bareground(g)
               else
                  ! initialize static layers for reduced complexity FATES versions from HLM
                  this%fates(nc)%bc_in(s)%pft_areafrac(:)=0._r8
                  do m = surfpft_lb,surfpft_ub
                     ft = m-surfpft_lb
                     this%fates(nc)%bc_in(s)%pft_areafrac(ft)=wt_nat_patch(g,t,m)
                  end do

                  if(abs(sum(this%fates(nc)%bc_in(s)%pft_areafrac(surfpft_lb:surfpft_ub))-1.0_r8).gt.1.0e-9)then
                     write(iulog,*) 'pft_area error in interfc ',s, sum(this%fates(nc)%bc_in(s)%pft_areafrac(:))-1.0_r8
                     call endrun(msg=errMsg(sourcefile, __LINE__))
                  end if
               end if
            end if

         end do

         ! Initialize site-level static quantities dictated by the HLM
         ! currently ground layering depth

         call this%init_soil_depths(nc)

         if (use_fates_planthydro) then
            call InitHydrSites(this%fates(nc)%sites,this%fates(nc)%bc_in)
         end if

         if( this%fates(nc)%nsites == 0 ) then
            write(iulog,*) 'Clump ',nc,' had no valid FATES sites'
            write(iulog,*) 'This will likely cause problems until code is improved'
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if

         ! Set patch itypes on natural veg columns to nonsense
         ! This will force a crash if the model outside of FATES tries to think
         ! of the patch as a PFT.

         do s = 1, this%fates(nc)%nsites
            c = this%f2hmap(nc)%fcolumn(s)
            pi = col_pp%pfti(c)+1
            pf = col_pp%pftf(c)
            veg_pp%is_fates(pi:pf) = .true.
         end do

      end do
      !$OMP END PARALLEL DO

      call this%init_history_io(bounds_proc)

      ! Report Fates Parameters (debug flag in lower level routines)
      call FatesReportParameters(masterproc)
      
      ! Fire data to send to FATES
      call create_fates_fire_data_method( this%fates_fire_data_method )

      ! deallocate the local landuse x pft array
      if (use_fates_fixed_biogeog .and. use_fates_luh) then
         deallocate(landuse_pft_map)
         deallocate(landuse_bareground)
      end if

    end subroutine init

    ! ===================================================================================

    subroutine check_hlm_active(this, nc, bounds_clump)

      ! ---------------------------------------------------------------------------------
      ! This subroutine is not currently used.  It is just a utility that may come
      ! in handy when we have dynamic sites in FATES
      ! ---------------------------------------------------------------------------------

      implicit none
      class(hlm_fates_interface_type), intent(inout) :: this
      integer                                        :: nc
      type(bounds_type),intent(in)                   :: bounds_clump

      ! local variables
      integer :: c

      do c = bounds_clump%begc,bounds_clump%endc

         ! FATES ACTIVE BUT HLM IS NOT
         if(this%f2hmap(nc)%hsites(c)>0 .and. .not.col_pp%active(c)) then

            write(iulog,*) 'INACTIVE COLUMN WITH ACTIVE FATES SITE'
            write(iulog,*) 'c = ',c
            call endrun(msg=errMsg(sourcefile, __LINE__))

         elseif (this%f2hmap(nc)%hsites(c)==0 .and. col_pp%active(c)) then

            write(iulog,*) 'ACTIVE COLUMN WITH INACTIVE FATES SITE'
            write(iulog,*) 'c = ',c
            call endrun(msg=errMsg(sourcefile, __LINE__))
         end if
      end do

   end subroutine check_hlm_active

   ! ------------------------------------------------------------------------------------

   subroutine dynamics_driv(this, bounds_clump, top_as_inst,          &
         top_af_inst, atm2lnd_inst, soilstate_inst, &
         canopystate_inst, frictionvel_inst, soil_water_retention_curve )

      use FatesConstantsMod       , only : m2_per_km2

      ! This wrapper is called daily from clm_driver
      ! This wrapper calls ed_driver, which is the daily dynamics component of FATES
      ! ed_driver is not a hlm_fates_inst_type procedure because we need an extra step
      ! to process array bounding information

      implicit none
      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type),intent(in)                   :: bounds_clump
      type(topounit_atmospheric_state), intent(in)   :: top_as_inst
      type(topounit_atmospheric_flux),  intent(in)   :: top_af_inst
      type(atm2lnd_type)      , intent(in)           :: atm2lnd_inst
      type(soilstate_type)    , intent(in)           :: soilstate_inst
      type(canopystate_type)  , intent(inout)        :: canopystate_inst
      type(frictionvel_type)  , intent(inout)        :: frictionvel_inst
      class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
      
      ! !LOCAL VARIABLES:
      integer  :: s, i                     ! site index
      integer  :: c                        ! column index (HLM)
      integer  :: j                        ! Soil layer index
      integer  :: t                        ! topounit index (HLM)
      integer  :: ifp                      ! patch index
      integer  :: ft                       ! patch functional type index
      integer  :: p                        ! HLM patch index
      integer  :: g                        ! HLM grid index
      integer  :: nc                       ! clump index
      integer  :: nlevsoil                 ! number of soil layers at the site
      integer  :: ier                      ! allocate status code
      real(r8) :: s_node, smp_node         ! local for relative water content and potential
      
      real(r8), pointer :: lnfm24(:)       ! 24-hour averaged lightning data
      real(r8), pointer :: gdp_lf_col(:)          ! gdp data

      logical  :: do_landuse_update        ! local flag to pass transitions update to fates

      !-----------------------------------------------------------------------


      if (masterproc) then
         write(iulog, *) 'FATES dynamics start'
      end if


      nc = bounds_clump%clump_index

      ! ---------------------------------------------------------------------------------
      ! Part I.
      ! Prepare input boundary conditions for FATES dynamics
      ! Note that timing information is the same across all sites, this may
      ! seem redundant, but it is possible that we may have asynchronous site simulations
      ! one day.  The cost of holding site level boundary conditions is minimal
      ! and it keeps all the boundaries in one location
      ! ---------------------------------------------------------------------------------

      ! Set the FATES global time and date variables
      call GetAndSetTime
      
      if (fates_spitfire_mode > scalar_lightning) then
         allocate(lnfm24(bounds_clump%begg:bounds_clump%endg), stat=ier)
         if (ier /= 0) then
            call endrun(msg="allocation error for lnfm24"//&
                 errmsg(sourcefile, __LINE__))
         endif
         lnfm24 = this%fates_fire_data_method%GetLight24()
      end if
      
      if (fates_spitfire_mode == anthro_suppression) then
         allocate(gdp_lf_col(bounds_clump%begc:bounds_clump%endc), stat=ier)
         if (ier /= 0) then
            call endrun(msg="allocation error for gdp"//&
                 errmsg(sourcefile, __LINE__))
         endif
         gdp_lf_col = this%fates_fire_data_method%GetGDP()
      end if

      do s=1,this%fates(nc)%nsites

         c = this%f2hmap(nc)%fcolumn(s)
         g = col_pp%gridcell(c)
         t = col_pp%topounit(c)
         
         if (fates_spitfire_mode > scalar_lightning) then
            do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
               
               this%fates(nc)%bc_in(s)%lightning24(ifp) = lnfm24(g) * 24._r8  ! #/km2/hr to #/km2/day
               
               if (fates_spitfire_mode >= anthro_ignitions) then
                  this%fates(nc)%bc_in(s)%pop_density(ifp) = this%fates_fire_data_method%forc_hdm(g)
               end if

            end do ! ifp

            if (fates_spitfire_mode == anthro_suppression) then
               ! Placeholder for future fates use of gdp - comment out before integration
               !this%fates(nc)%bc_in(s)%gdp = gdp_lf_col(c) ! k US$/capita(g)
            end if
         end if

         nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

         ! Decomposition fluxes
         this%fates(nc)%bc_in(s)%w_scalar_sisl(1:nlevsoil) = col_cf%w_scalar(c,1:nlevsoil)
         this%fates(nc)%bc_in(s)%t_scalar_sisl(1:nlevsoil) = col_cf%t_scalar(c,1:nlevsoil)

         ! Soil water
         this%fates(nc)%bc_in(s)%h2o_liqvol_sl(1:nlevsoil)  = &
               col_ws%h2osoi_vol(c,1:nlevsoil)

         this%fates(nc)%bc_in(s)%max_rooting_depth_index_col = &
              min(nlevsoil, canopystate_inst%altmax_lastyear_indx_col(c))

         do j = 1,nlevsoil
            this%fates(nc)%bc_in(s)%tempk_sl(j) = col_es%t_soisno(c,j)
         end do

         call get_active_suction_layers(this%fates(nc)%nsites, &
             this%fates(nc)%sites,  &
             this%fates(nc)%bc_in,  &
             this%fates(nc)%bc_out)

         do j = 1,nlevsoil
            if(this%fates(nc)%bc_out(s)%active_suction_sl(j)) then
               s_node = max(col_ws%h2osoi_liqvol(c,j)/soilstate_inst%eff_porosity_col(c,j) ,0.01_r8)
               call soil_water_retention_curve%soil_suction(soilstate_inst%sucsat_col(c,j), &
                       s_node, &
                       soilstate_inst%bsw_col(c,j), &
                       smp_node)
               this%fates(nc)%bc_in(s)%smp_sl(j) = smp_node
            end if
         end do
         
         do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
            p = ifp+col_pp%pfti(c)

            this%fates(nc)%bc_in(s)%precip24_pa(ifp) = &
                  top_af_inst%prec24h(t)

            this%fates(nc)%bc_in(s)%relhumid24_pa(ifp) = &
                  top_as_inst%rh24h(t)

            this%fates(nc)%bc_in(s)%wind24_pa(ifp) = &
                  top_as_inst%wind24h(t)

         end do

         if(use_fates_sp)then
            do ft = surfpft_lb,surfpft_ub !set of pfts in HLM
                ! here we are mapping from P space in the HLM to FT space in the sp_input arrays.
                p = ft + col_pp%pfti(c) ! for an FT of 1 we want to use
                this%fates(nc)%bc_in(s)%hlm_sp_tlai(ft) = canopystate_inst%tlai_patch(p)
                this%fates(nc)%bc_in(s)%hlm_sp_tsai(ft) = canopystate_inst%tsai_patch(p)
                this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = canopystate_inst%htop_patch(p)
                if(canopystate_inst%htop_patch(p).lt.1.0e-20)then ! zero htop causes inifinite/nans.
                  this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = 0.01_r8
                endif
            end do ! p
          end if ! SP

         if(use_fates_planthydro)then
            this%fates(nc)%bc_in(s)%hksat_sisl(1:nlevsoil)  = soilstate_inst%hksat_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil) = soilstate_inst%watsat_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%watres_sisl(1:nlevsoil) = soilstate_inst%watmin_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil) = soilstate_inst%sucsat_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil)    = soilstate_inst%bsw_col(c,1:nlevsoil)
            this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil) =  col_ws%h2osoi_liq(c,1:nlevsoil)
         end if

         ! get the harvest data, which is by gridcell
         ! for now there is one veg column per gridcell, so store all harvest data in each site
         ! this will eventually change
         ! the harvest data are zero if today is before the start of the harvest time series
         if (trim(fates_harvest_mode) == fates_harvest_hlmlanduse) then
            this%fates(nc)%bc_in(s)%hlm_harvest_rates = harvest_rates(:,g)
            this%fates(nc)%bc_in(s)%hlm_harvest_catnames = harvest_varnames
            this%fates(nc)%bc_in(s)%hlm_harvest_units = wood_harvest_units
         else if (trim(fates_harvest_mode) == fates_harvest_luh_area .or. &
                  trim(fates_harvest_mode) == fates_harvest_luh_mass) then
            this%fates(nc)%bc_in(s)%hlm_harvest_rates = landuse_harvest(:,g)
            this%fates(nc)%bc_in(s)%hlm_harvest_catnames = landuse_harvest_varnames
            this%fates(nc)%bc_in(s)%hlm_harvest_units = landuse_harvest_units
         end if
         this%fates(nc)%bc_in(s)%site_area=col_pp%wtgcell(c)*grc_pp%area(g)*m2_per_km2

         if (use_fates_luh) then
            this%fates(nc)%bc_in(s)%hlm_luh_states = landuse_states(:,g)
            this%fates(nc)%bc_in(s)%hlm_luh_state_names = landuse_state_varnames
            this%fates(nc)%bc_in(s)%hlm_luh_transitions = landuse_transitions(:,g)
            this%fates(nc)%bc_in(s)%hlm_luh_transition_names = landuse_transition_varnames
         end if

      end do

      ! Nutrient uptake fluxes have been accumulating with each short
      ! timestep, here, we unload them from the boundary condition
      ! structures into the cohort structures.
      call UnPackNutrientAquisitionBCs(this%fates(nc)%sites, this%fates(nc)%bc_in)
      
      ! Distribute any seeds from neighboring gridcells into the current gridcell
      ! Global seed availability array populated by WrapGlobalSeedDispersal call
      if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then
         call this%WrapUpdateFatesSeedInOut(bounds_clump)
      end if

      ! ---------------------------------------------------------------------------------
      ! Flush arrays to values defined by %flushval (see registry entry in
      ! subroutine define_history_vars()
      ! ---------------------------------------------------------------------------------
      call fates_hist%flush_hvars(nc,upfreq_in=group_dyna_simple)
      call fates_hist%flush_hvars(nc,upfreq_in=group_dyna_complx)

      ! ---------------------------------------------------------------------------------
      ! Part II: Call the FATES model now that input boundary conditions have been
      ! provided.
      ! ---------------------------------------------------------------------------------

      do s = 1,this%fates(nc)%nsites
           
            call ed_ecosystem_dynamics(this%fates(nc)%sites(s),    &
                  this%fates(nc)%bc_in(s), &
                  this%fates(nc)%bc_out(s))

            call ed_update_site(this%fates(nc)%sites(s), &
                  this%fates(nc)%bc_in(s), &
                  this%fates(nc)%bc_out(s), &
                  is_restarting = .false.)
      enddo

      ! ---------------------------------------------------------------------------------
      ! Part III
      ! Update diagnostics of the FATES ecosystem structure that are used in the HLM.
      ! ---------------------------------------------------------------------------------
      call this%wrap_update_hlmfates_dyn(nc,               &
                                         bounds_clump,     &
                                         canopystate_inst, &
                                         frictionvel_inst, .false.)

      ! ---------------------------------------------------------------------------------
      ! Part IV:
      ! Update history IO fields that depend on ecosystem dynamics
      ! ---------------------------------------------------------------------------------
      call fates_hist%update_history_dyn( nc,                    &
           this%fates(nc)%nsites, &
           this%fates(nc)%sites,  &
           this%fates(nc)%bc_in)

      if (masterproc) then
         write(iulog, *) 'FATES dynamics complete'
      end if


      return
   end subroutine dynamics_driv


   ! ====================================================================================

   !subroutine WrapFatesNutrientSources(this,bounds_clump)


     ! Add in any N source from FATES
     !sminn_vr(c,:)    = sminn_vr(c,:) + alm_fates(ci)%bc_out(s)%source_nh4(:)
     !smin_nh4_vr(c,:) = smin_nh4_vr(c,:) + alm_fates(ci)%bc_out(s)%source_nh4(:)


    ! return
   !end subroutine WrapFatesNutrientSources


   ! ====================================================================================

   subroutine UpdateLitterFluxes(this,bounds_clump)

      implicit none
      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type)              , intent(in)    :: bounds_clump

      ! !LOCAL VARIABLES:
      integer  :: s                        ! site index
      integer  :: c                        ! column index (HLM)
      integer  :: nc                       ! clump index
      integer  :: nld_si
      real(r8) :: dtime

      dtime = real(get_step_size(),r8)
      nc = bounds_clump%clump_index

      do s = 1, this%fates(nc)%nsites
         c = this%f2hmap(nc)%fcolumn(s)

         col_cf%decomp_cpools_sourcesink(c,1:nlevdecomp,i_met_lit) = &
              col_cf%decomp_cpools_sourcesink(c,1:nlevdecomp,i_met_lit) + &
              this%fates(nc)%bc_out(s)%litt_flux_lab_c_si(1:nlevdecomp) * dtime
         col_cf%decomp_cpools_sourcesink(c,1:nlevdecomp,i_cel_lit) = &
              col_cf%decomp_cpools_sourcesink(c,1:nlevdecomp,i_cel_lit) + &
              this%fates(nc)%bc_out(s)%litt_flux_cel_c_si(1:nlevdecomp)* dtime
         col_cf%decomp_cpools_sourcesink(c,1:nlevdecomp,i_lig_lit) = &
              col_cf%decomp_cpools_sourcesink(c,1:nlevdecomp,i_lig_lit) + &
              this%fates(nc)%bc_out(s)%litt_flux_lig_c_si(1:nlevdecomp) * dtime

         col_cf%litfall(c) = &
               sum(this%fates(nc)%bc_out(s)%litt_flux_lab_c_si(1:nlevdecomp) * this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp)) + &
               sum(this%fates(nc)%bc_out(s)%litt_flux_cel_c_si(1:nlevdecomp) * this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp)) + &
               sum(this%fates(nc)%bc_out(s)%litt_flux_lig_c_si(1:nlevdecomp) * this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp))


         ! Since N and P are always allocated in ELM, AND, since on the FATES
         ! side we have prepped these arrays, which may be zero fluxes in the case of
         ! prescribed FATES nutrient mode, we can send the fluxes into the source pools

         select case(fates_parteh_mode)
         case (prt_cnp_flex_allom_hyp )

            col_pf%decomp_ppools_sourcesink(c,1:nlevdecomp,i_met_lit) = &
                 col_pf%decomp_ppools_sourcesink(c,1:nlevdecomp,i_met_lit) + &
                 this%fates(nc)%bc_out(s)%litt_flux_lab_p_si(1:nlevdecomp) * dtime

            col_pf%decomp_ppools_sourcesink(c,1:nlevdecomp,i_cel_lit) = &
                 col_pf%decomp_ppools_sourcesink(c,1:nlevdecomp,i_cel_lit) + &
                 this%fates(nc)%bc_out(s)%litt_flux_cel_p_si(1:nlevdecomp)* dtime

            col_pf%decomp_ppools_sourcesink(c,1:nlevdecomp,i_lig_lit) = &
                 col_pf%decomp_ppools_sourcesink(c,1:nlevdecomp,i_lig_lit) + &
                 this%fates(nc)%bc_out(s)%litt_flux_lig_p_si(1:nlevdecomp) * dtime

            ! Diagnostic for mass balancing (gP/m2/s)
            col_pf%plant_to_litter_pflux(c) = &
                 sum(this%fates(nc)%bc_out(s)%litt_flux_lab_p_si(1:nlevdecomp)*this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp)) + &
                 sum(this%fates(nc)%bc_out(s)%litt_flux_cel_p_si(1:nlevdecomp)*this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp)) + &
                 sum(this%fates(nc)%bc_out(s)%litt_flux_lig_p_si(1:nlevdecomp)*this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp))

            ! Transfer Nitrogen
            col_nf%decomp_npools_sourcesink(c,1:nlevdecomp,i_met_lit) = &
                 col_nf%decomp_npools_sourcesink(c,1:nlevdecomp,i_met_lit) + &
                 this%fates(nc)%bc_out(s)%litt_flux_lab_n_si(1:nlevdecomp) * dtime

            col_nf%decomp_npools_sourcesink(c,1:nlevdecomp,i_cel_lit) = &
                 col_nf%decomp_npools_sourcesink(c,1:nlevdecomp,i_cel_lit) + &
                 this%fates(nc)%bc_out(s)%litt_flux_cel_n_si(1:nlevdecomp)* dtime

            col_nf%decomp_npools_sourcesink(c,1:nlevdecomp,i_lig_lit) = &
                 col_nf%decomp_npools_sourcesink(c,1:nlevdecomp,i_lig_lit) + &
                 this%fates(nc)%bc_out(s)%litt_flux_lig_n_si(1:nlevdecomp) * dtime

            ! Diagnostic for mass balancing  (gN/m2/s)
            col_nf%plant_to_litter_nflux(c) = &
                 sum(this%fates(nc)%bc_out(s)%litt_flux_lab_n_si(1:nlevdecomp)*this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp)) + &
                 sum(this%fates(nc)%bc_out(s)%litt_flux_cel_n_si(1:nlevdecomp)*this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp)) + &
                 sum(this%fates(nc)%bc_out(s)%litt_flux_lig_n_si(1:nlevdecomp)*this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp))

         end select

      end do

   end subroutine UpdateLitterFluxes

   !--------------------------------------------------------------------------------------

   subroutine wrap_update_hlmfates_dyn(this, nc, bounds_clump,      &
         canopystate_inst, frictionvel_inst, is_initing_from_restart)

      ! ---------------------------------------------------------------------------------
      ! This routine handles the updating of vegetation canopy diagnostics, (such as lai)
      ! that either requires HLM boundary conditions (like snow accumulation) or
      ! provides boundary conditions (such as vegetation fractional coverage)
      ! ---------------------------------------------------------------------------------

     implicit none
     class(hlm_fates_interface_type), intent(inout) :: this
     type(bounds_type),intent(in)                   :: bounds_clump
     integer                 , intent(in)           :: nc
     type(canopystate_type)  , intent(inout)        :: canopystate_inst
     type(frictionvel_type)  , intent(inout)        :: frictionvel_inst

     ! is this being called during a read from restart sequence (if so then use the restarted fates
     ! snow depth variable rather than the CLM variable).
     logical                 , intent(in)           :: is_initing_from_restart

     logical :: dispersal_flag ! local flag to pass to the inside of the site loop

     integer :: npatch  ! number of patches in each site
     integer :: ifp     ! index FATES patch
     integer :: p       ! HLM patch index
     integer :: s       ! site index
     integer :: c       ! column index
     integer :: g       ! gridcell index

     real(r8) :: areacheck


     associate(                                &
         tlai => canopystate_inst%tlai_patch , &
         elai => canopystate_inst%elai_patch , &
         tsai => canopystate_inst%tsai_patch , &
         esai => canopystate_inst%esai_patch , &
         htop => canopystate_inst%htop_patch , &
         hbot => canopystate_inst%hbot_patch , &
         z0m  => frictionvel_inst%z0m_patch  , & ! Output: [real(r8) (:)   ] momentum roughness length (m)
         displa => canopystate_inst%displa_patch, &
         dleaf_patch => canopystate_inst%dleaf_patch, &
         snow_depth => col_ws%snow_depth, &
         frac_sno_eff => col_ws%frac_sno_eff, &
         frac_veg_nosno_alb => canopystate_inst%frac_veg_nosno_alb_patch)


       ! Process input boundary conditions to FATES
       ! --------------------------------------------------------------------------------
       do s=1,this%fates(nc)%nsites
          c = this%f2hmap(nc)%fcolumn(s)
          this%fates(nc)%bc_in(s)%snow_depth_si   = snow_depth(c)
          this%fates(nc)%bc_in(s)%frac_sno_eff_si = frac_sno_eff(c)
       end do

       ! Only update the fates internal snow burial if this is not a restart
       if (.not. is_initing_from_restart) then
          call UpdateFatesAvgSnowDepth(this%fates(nc)%sites,this%fates(nc)%bc_in)
       end if

       ! Canopy diagnostics for FATES
       call canopy_summarization(this%fates(nc)%nsites, &
            this%fates(nc)%sites,  &
            this%fates(nc)%bc_in)

       ! Canopy diagnostic outputs for HLM, including LUC
       call update_hlm_dynamics(this%fates(nc)%nsites, &
            this%fates(nc)%sites,  &
            this%f2hmap(nc)%fcolumn, &
            this%fates(nc)%bc_out )

       !---------------------------------------------------------------------------------
       ! Diagnose water storage in canopy if hydraulics is on
       if ( use_fates_planthydro ) then
          !pass the water storage in plants back to the HLM
          do s = 1, this%fates(nc)%nsites
             c = this%f2hmap(nc)%fcolumn(s)
             col_ws%total_plant_stored_h2o(c) = &
                  this%fates(nc)%bc_out(s)%plant_stored_h2o_si
          end do
       end if
       !---------------------------------------------------------------------------------

       ! Convert FATES dynamics into HLM usable information
       ! Initialize weighting variables (note FATES is the only HLM module
       ! that uses "is_veg" and "is_bareground".  The entire purpose of these
       ! variables is to inform patch%wtcol(p).  wt_ed is imposed on wtcol,
       ! but only for FATES columns.

       ! Check if seed dispersal mode is 'turned on', if not return to calling procedure
       if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then
          ! zero the outgoing seed array
          this%fates_seed%outgoing_local(:,:) = 0._r8
          dispersal_flag = .false.
          if (IsItDispersalTime()) dispersal_flag = .true.
       end if
       
       do s = 1,this%fates(nc)%nsites

          c = this%f2hmap(nc)%fcolumn(s)
          g = col_pp%gridcell(c)
          
          ! Accumulate seeds from sites to the gridcell local outgoing buffer
          if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then
             if (dispersal_flag) this%fates_seed%outgoing_local(g,:) = this%fates(nc)%sites(s)%seed_out(:)
          end if

          veg_pp%is_veg(col_pp%pfti(c):col_pp%pftf(c))        = .false.
          veg_pp%is_bareground(col_pp%pfti(c):col_pp%pftf(c)) = .false.
          veg_pp%wt_ed(col_pp%pfti(c):col_pp%pftf(c))         = 0.0_r8


          ! Other modules may have AI's we only flush values
          ! that are on the naturally vegetated columns
          elai(col_pp%pfti(c):col_pp%pftf(c)) = 0.0_r8
          esai(col_pp%pfti(c):col_pp%pftf(c)) = 0.0_r8
          hbot(col_pp%pfti(c):col_pp%pftf(c)) = 0.0_r8

          if(use_fates_sp)then
            canopystate_inst%tlai_hist_patch(col_pp%pfti(c):col_pp%pftf(c)) = 0.0_r8
            canopystate_inst%tsai_hist_patch(col_pp%pfti(c):col_pp%pftf(c)) = 0.0_r8
            canopystate_inst%htop_hist_patch(col_pp%pfti(c):col_pp%pftf(c)) = 0.0_r8
          else
            tlai(col_pp%pfti(c):col_pp%pftf(c)) = 0.0_r8
            tsai(col_pp%pfti(c):col_pp%pftf(c)) = 0.0_r8
            htop(col_pp%pfti(c):col_pp%pftf(c)) = 0.0_r8
          endif

          ! FATES does not dictate bare-ground so turbulent
          ! variables are not over-written.
          z0m(col_pp%pfti(c)+1:col_pp%pftf(c)) = 0.0_r8
          displa(col_pp%pfti(c)+1:col_pp%pftf(c)) = 0.0_r8
          dleaf_patch(col_pp%pfti(c)+1:col_pp%pftf(c)) = 0.0_r8

          frac_veg_nosno_alb(col_pp%pfti(c):col_pp%pftf(c)) = 0

          ! Set the bareground patch indicator
          veg_pp%is_bareground(col_pp%pfti(c)) = .true.
          npatch = this%fates(nc)%sites(s)%youngest_patch%patchno

          ! Precision errors on the canopy_fraction_pa sum, even small (e-12)
          ! do exist, and can create potentially negetive bare-soil fractions
          ! (ie -1e-12 or smaller). Even though this is effectively zero,
          ! it can generate weird logic scenarios in the ctsm/elm code, so we
          ! protext it here with a lower bound of 0.0_r8.

          veg_pp%wt_ed(col_pp%pfti(c)) = max(0.0_r8, &
               1.0_r8 - sum(this%fates(nc)%bc_out(s)%canopy_fraction_pa(1:npatch)) )

          ! initialize SP mode pft order index to 0.  Below ground is the 0th patch
          veg_pp%sp_pftorder_index(col_pp%pfti(c)) = 0
          areacheck = veg_pp%wt_ed(col_pp%pfti(c))

          do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno

             p = ifp+col_pp%pfti(c)

             ! bc_out(s)%canopy_fraction_pa(ifp) is the area fraction
             ! the site's total ground area that is occupied by the
             ! area footprint of the current patch's vegetation canopy

             veg_pp%is_veg(p) = .true.
             veg_pp%wt_ed(p)  = this%fates(nc)%bc_out(s)%canopy_fraction_pa(ifp)

             areacheck = areacheck + veg_pp%wt_ed(p)

             elai(p) = this%fates(nc)%bc_out(s)%elai_pa(ifp)
             esai(p) = this%fates(nc)%bc_out(s)%esai_pa(ifp)
             hbot(p) = this%fates(nc)%bc_out(s)%hbot_pa(ifp)

             if(use_fates_sp)then
               canopystate_inst%tlai_hist_patch(p) = this%fates(nc)%bc_out(s)%tlai_pa(ifp)
               canopystate_inst%tsai_hist_patch(p) = this%fates(nc)%bc_out(s)%tsai_pa(ifp)
               canopystate_inst%htop_hist_patch(p) = this%fates(nc)%bc_out(s)%htop_pa(ifp)
             else
               tlai(p) = this%fates(nc)%bc_out(s)%tlai_pa(ifp)
               tsai(p) = this%fates(nc)%bc_out(s)%tsai_pa(ifp)
               htop(p) = this%fates(nc)%bc_out(s)%htop_pa(ifp)
             endif

             if(use_fates_sp.and.abs(canopystate_inst%tlai_hist_patch(p) - &
                this%fates(nc)%bc_out(s)%tlai_pa(ifp)).gt.1e-09)then
                write(iulog,*) 'fates lai not like hlm lai',tlai(p),this%fates(nc)%bc_out(s)%tlai_pa(ifp),ifp
                call endrun(msg=errMsg(sourcefile, __LINE__))
             endif

             frac_veg_nosno_alb(p) = this%fates(nc)%bc_out(s)%frac_veg_nosno_alb_pa(ifp)

             ! Note that while we pass the following values at this point
             ! we have to send the same values after each time-step because
             ! the HLM keeps changing the value and re-setting, so we
             ! re-send instead of re-set. See alm_fates%TransferZ0mDisp()
             z0m(p)    = this%fates(nc)%bc_out(s)%z0m_pa(ifp)
             displa(p) = this%fates(nc)%bc_out(s)%displa_pa(ifp)
             dleaf_patch(p) = this%fates(nc)%bc_out(s)%dleaf_pa(ifp)
             
             


          end do

          if(abs(areacheck - 1.0_r8).gt.1.e-9_r8)then
            write(iulog,*) 'area wrong in interface',areacheck - 1.0_r8
            call endrun(msg=errMsg(sourcefile, __LINE__))
          endif

       end do
     end associate
   end subroutine wrap_update_hlmfates_dyn

   ! ====================================================================================

   subroutine restart( this, bounds_proc, ncid, flag, &
                             canopystate_inst, frictionvel_inst, soilstate_inst )

      ! ---------------------------------------------------------------------------------
      ! The ability to restart the model is handled through three different types of calls
      ! "Define" the variables in the restart file, we "read" those variables into memory
      ! or "write" data into the file from memory.  This subroutine accomodates all three
      ! of those modes through the "flag" argument.  FATES as an external model also
      ! requires an initialization step, where we set-up the dimensions, allocate and
      ! flush the memory space that is used to transfer data in and out of the file.  This
      ! Only occurs once, where as the define step occurs every time a file is opened.
      !
      ! Note: canopystate_inst are arguments only because following
      ! the reading of variables, it is necessary to update diagnostics of the canopy
      ! throug the interface call alm_fates%wrap_update_hlmfates_dyn() which requires
      ! this information from the HLM.
      ! ---------------------------------------------------------------------------------


     use FatesConstantsMod, only : fates_long_string_length
     use FatesIODimensionsMod, only: fates_bounds_type
     use FatesIOVariableKindMod, only : site_r8, site_int, cohort_r8, cohort_int
     use EDMainMod, only :        ed_update_site
     use FatesInterfaceTypesMod, only:  fates_maxElementsPerSite

      implicit none

      ! Arguments

      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type)              , intent(in)    :: bounds_proc
      type(file_desc_t)              , intent(inout) :: ncid    ! netcdf id
      character(len=*)               , intent(in)    :: flag
      type(canopystate_type)         , intent(inout) :: canopystate_inst
      type(frictionvel_type)         , intent(inout) :: frictionvel_inst
      type(soilstate_type)           , intent(inout) :: soilstate_inst

      ! Locals
      type(bounds_type) :: bounds_clump
      integer           :: nc
      integer           :: nclumps
      type(fates_bounds_type) :: fates_bounds
      type(fates_bounds_type) :: fates_clump
      integer                 :: c   ! HLM column index
      integer                 :: s   ! Fates site index
      integer                 :: g   ! HLM grid index
      integer                 :: p   ! HLM patch index
      integer                 :: ft  ! plant functional type
      integer                 :: dk_index
      character(len=fates_long_string_length) :: ioname
      integer                 :: nvar
      integer                 :: nlevsoil
      integer                 :: ivar
      logical                 :: readvar

      logical, save           :: initialized = .false.


      nclumps = get_proc_clumps()

      ! ---------------------------------------------------------------------------------
      ! note (rgk: 11-2016) The history and restart intialization process assumes
      ! that the number of site/columns active is a static entity.  Thus
      ! we only allocate the mapping tables for the column/sites we start with.
      ! If/when we start having dynamic column/sites (for reasons uknown as of yet)
      ! we will need to re-evaluate the allocation of the mapping tables so they
      ! can be unallocated,reallocated and set every time a new column/site is spawned
      ! ---------------------------------------------------------------------------------

      ! ---------------------------------------------------------------------------------
      ! Only initialize the FATES restart structures the first time it is called
      ! Note that the allocations involved with initialization are static.
      ! This is because the array spaces for IO span the entire column, patch and cohort
      ! range on the proc.
      ! With DYNAMIC LANDUNITS or SPAWNING NEW OR CULLING OLD SITES:
      ! we will in that case have to de-allocate, reallocate and then re-set the mapping
      ! tables:  this%fates_restart%restart_map(nc)
      ! I think that is it...
      ! ---------------------------------------------------------------------------------

      ! Set the FATES global time and date variables
      call GetAndSetTime

      if(.not.initialized) then

         initialized=.true.

         ! ------------------------------------------------------------------------------
         ! PART I: Set FATES DIMENSIONING INFORMATION
         ! ------------------------------------------------------------------------------

         call hlm_bounds_to_fates_bounds(bounds_proc, fates_bounds)

         call this%fates_restart%Init(nclumps, fates_bounds)

         ! Define the bounds on the first dimension for each thread
         !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,fates_clump)
         do nc = 1,nclumps
            call get_clump_bounds(nc, bounds_clump)

            ! thread bounds for patch
            call hlm_bounds_to_fates_bounds(bounds_clump, fates_clump)
            call this%fates_restart%SetThreadBoundsEach(nc, fates_clump)
         end do
         !$OMP END PARALLEL DO

         !$OMP PARALLEL DO PRIVATE (nc,s,c,g)
         do nc = 1,nclumps

            allocate(this%fates_restart%restart_map(nc)%site_index(this%fates(nc)%nsites))
            allocate(this%fates_restart%restart_map(nc)%cohort1_index(this%fates(nc)%nsites))
            do s=1,this%fates(nc)%nsites
               c = this%f2hmap(nc)%fcolumn(s)
               this%fates_restart%restart_map(nc)%site_index(s)   = c
               g = col_pp%gridcell(c)
               this%fates_restart%restart_map(nc)%cohort1_index(s) = (g-1)*fates_maxElementsPerSite + 1
            end do

         end do
         !$OMP END PARALLEL DO

         ! ------------------------------------------------------------------------------------
         ! PART II: USE THE JUST DEFINED DIMENSIONS TO ASSEMBLE THE VALID IO TYPES
         ! INTERF-TODO: THESE CAN ALL BE EMBEDDED INTO A SUBROUTINE IN HISTORYIOMOD
         ! ------------------------------------------------------------------------------------
         call this%fates_restart%assemble_restart_output_types()


         ! ------------------------------------------------------------------------------------
         ! PART III: DEFINE THE LIST OF OUTPUT VARIABLE OBJECTS, AND REGISTER THEM WITH THE
         ! HLM ACCORDING TO THEIR TYPES
         ! ------------------------------------------------------------------------------------
         call this%fates_restart%initialize_restart_vars()

      end if

      ! ---------------------------------------------------------------------------------
      ! If we are writing, we must loop through our linked list structures and transfer the
      ! information in the linked lists (FATES state memory) to the output vectors.
      ! ---------------------------------------------------------------------------------

      if(flag=='write')then
         !$OMP PARALLEL DO PRIVATE (nc)
         do nc = 1, nclumps
            if (this%fates(nc)%nsites>0) then
               call this%fates_restart%set_restart_vectors(nc,this%fates(nc)%nsites, &
                                                           this%fates(nc)%sites)
            end if
         end do
         !$OMP END PARALLEL DO
      end if

      ! ---------------------------------------------------------------------------------
      ! In all cases, iterate through the list of variable objects
      ! and either define, write or read to the NC buffer
      ! This seems strange, but keep in mind that the call to restartvar()
      ! has a different function in all three cases.
      ! ---------------------------------------------------------------------------------

      nvar = this%fates_restart%num_restart_vars()
      do ivar = 1, nvar

         associate( vname => this%fates_restart%rvars(ivar)%vname, &
              vunits      => this%fates_restart%rvars(ivar)%units,   &
              vlong       => this%fates_restart%rvars(ivar)%long )

           dk_index = this%fates_restart%rvars(ivar)%dim_kinds_index
           ioname = trim(this%fates_restart%dim_kinds(dk_index)%name)

           select case(trim(ioname))
           case(cohort_r8)

              call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
                    xtype=ncd_double,dim1name=trim('cohort'),long_name=trim(vlong), &
                    units=trim(vunits),interpinic_flag='interp', &
                    data=this%fates_restart%rvars(ivar)%r81d,readvar=readvar)

           case(site_r8)

              call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
                    xtype=ncd_double,dim1name=trim('column'),long_name=trim(vlong), &
                    units=trim(vunits),interpinic_flag='interp', &
                    data=this%fates_restart%rvars(ivar)%r81d,readvar=readvar)

           case(cohort_int)

              call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
                    xtype=ncd_int,dim1name=trim('cohort'),long_name=trim(vlong), &
                    units=trim(vunits),interpinic_flag='interp', &
                    data=this%fates_restart%rvars(ivar)%int1d,readvar=readvar)

           case(site_int)

              call restartvar(ncid=ncid, flag=flag, varname=trim(vname), &
                    xtype=ncd_int,dim1name=trim('column'),long_name=trim(vlong), &
                    units=trim(vunits),interpinic_flag='interp', &
                    data=this%fates_restart%rvars(ivar)%int1d,readvar=readvar)

           case default
              write(iulog,*) 'A FATES iotype was created that was not registerred'
              write(iulog,*) 'in CLM.:',trim(ioname)
              call endrun(msg=errMsg(sourcefile, __LINE__))
           end select

         end associate
      end do

      ! ---------------------------------------------------------------------------------
      ! If we are in a read mode, then we have just populated the sparse vectors
      ! in the IO object list. The data in these vectors needs to be transferred
      ! to the linked lists to populate the state memory.
      ! ---------------------------------------------------------------------------------

      if(flag=='read')then

         !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,s)
         do nc = 1, nclumps
            if (this%fates(nc)%nsites>0) then

               call get_clump_bounds(nc, bounds_clump)

               ! ------------------------------------------------------------------------
               ! Convert newly read-in vectors into the FATES namelist state variables
               ! ------------------------------------------------------------------------
               call this%fates_restart%create_patchcohort_structure(nc, &
                    this%fates(nc)%nsites, this%fates(nc)%sites, this%fates(nc)%bc_in,&
                    this%fates(nc)%bc_out)

               call this%fates_restart%get_restart_vectors(nc, this%fates(nc)%nsites, &
                    this%fates(nc)%sites )



               ! I think ed_update_site and update_hlmfates_dyn are doing some similar
               ! update type stuff, should consolidate (rgk 11-2016)
               do s = 1,this%fates(nc)%nsites

                  c = this%f2hmap(nc)%fcolumn(s)
                  this%fates(nc)%bc_in(s)%max_rooting_depth_index_col = &
                       min(this%fates(nc)%bc_in(s)%nlevsoil, canopystate_inst%altmax_lastyear_indx_col(c))

                  ! When restarting the model, this subroutine has several
                  ! procedures that are incremental or don't need to be performed for
                  ! during the restart sequence. For the prior, we don't want the restarted
                  ! run to call these routines more than would had been called during
                  ! a continuous simulation period, as it would change results. So
                  ! we pass in the "is_restarting=.true." flag so we can bypass those procedures

                  call ed_update_site( this%fates(nc)%sites(s), &
                        this%fates(nc)%bc_in(s), &
                        this%fates(nc)%bc_out(s), &
                        is_restarting = .true.)

                  ! This call sends internal fates variables into the
                  ! output boundary condition structures. Note: this is called
                  ! internally in fates dynamics as well.

                  call FluxIntoLitterPools(this%fates(nc)%sites(s), &
                       this%fates(nc)%bc_in(s), &
                       this%fates(nc)%bc_out(s))
               end do

               if(use_fates_sp)then
                  do s = 1,this%fates(nc)%nsites
                     c = this%f2hmap(nc)%fcolumn(s)
                     do ft = surfpft_lb,surfpft_ub !set of pfts in HLM
                        ! here we are mapping from P space in the HLM to FT space in the sp_input arrays.
                        p = ft + col_pp%pfti(c) ! for an FT of 1 we want to use
                        this%fates(nc)%bc_in(s)%hlm_sp_tlai(ft) = canopystate_inst%tlai_patch(p)
                        this%fates(nc)%bc_in(s)%hlm_sp_tsai(ft) = canopystate_inst%tsai_patch(p)
                        this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = canopystate_inst%htop_patch(p)
                        if(canopystate_inst%htop_patch(p).lt.1.0e-20)then ! zero htop causes inifinite/nans. This is
                           this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = 0.01_r8
                        endif
                     end do ! p
                  end do ! c
                end if ! SP

               ! ------------------------------------------------------------------------
               ! Re-populate all the hydraulics variables that are dependent
               ! on the key hydro state variables and plant carbon/geometry
               ! ------------------------------------------------------------------------
               if (use_fates_planthydro) then

                  do s = 1,this%fates(nc)%nsites
                     c = this%f2hmap(nc)%fcolumn(s)
                     nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil
                     this%fates(nc)%bc_in(s)%hksat_sisl(1:nlevsoil) = &
                          soilstate_inst%hksat_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil) = &
                          soilstate_inst%watsat_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%watres_sisl(1:nlevsoil) = &
                          soilstate_inst%watmin_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil) = &
                          soilstate_inst%sucsat_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil) = &
                          soilstate_inst%bsw_col(c,1:nlevsoil)

                     this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil) = &
                          col_ws%h2osoi_liq(c,1:nlevsoil)

                  end do

                  call RestartHydrStates(this%fates(nc)%sites,  &
                       this%fates(nc)%nsites, &
                       this%fates(nc)%bc_in,  &
                       this%fates(nc)%bc_out)
               end if

               ! ------------------------------------------------------------------------
               ! Update diagnostics of FATES ecosystem structure used in HLM.
               ! ------------------------------------------------------------------------
               call this%wrap_update_hlmfates_dyn(nc,bounds_clump, &
                     canopystate_inst,frictionvel_inst, .true.)

               ! ------------------------------------------------------------------------
               ! Update the 3D patch level radiation absorption fractions
               ! ------------------------------------------------------------------------
               call this%fates_restart%update_3dpatch_radiation(this%fates(nc)%nsites, &
                                                                this%fates(nc)%sites, &
                                                                this%fates(nc)%bc_out)

               ! ------------------------------------------------------------------------
               ! Update history IO fields that depend on ecosystem dynamics
               ! ------------------------------------------------------------------------

               if(fates_history_dimlevel(2)>0) then
                  call fates_hist%flush_hvars(nc,upfreq_in=group_dyna_simple)
                  do s = 1,this%fates(nc)%nsites
                     call fates_hist%zero_site_hvars(this%fates(nc)%sites(s), &
                          upfreq_in=group_dyna_simple)
                  end do
                  if(fates_history_dimlevel(2)>1) then
                     call fates_hist%flush_hvars(nc,upfreq_in=group_dyna_complx)
                     do s = 1,this%fates(nc)%nsites
                        call fates_hist%zero_site_hvars(this%fates(nc)%sites(s), &
                             upfreq_in=group_dyna_complx)
                     end do
                  end if
               end if
               
               call fates_hist%update_history_dyn( nc, &
                    this%fates(nc)%nsites, &
                    this%fates(nc)%sites,  &
                    this%fates(nc)%bc_in)

            end if
         end do
         !$OMP END PARALLEL DO

         ! Disperse seeds
         if (fates_seeddisp_cadence /= fates_dispersal_cadence_none) then
            call this%WrapGlobalSeedDispersal(is_restart_flag=.true.)
         end if

      end if

      return
   end subroutine restart

   !=====================================================================================

   subroutine init_coldstart(this, canopystate_inst, soilstate_inst, frictionvel_inst)

     ! Arguments
     class(hlm_fates_interface_type), intent(inout) :: this
     type(canopystate_type)         , intent(inout) :: canopystate_inst
     type(soilstate_type)           , intent(inout) :: soilstate_inst
     type(frictionvel_type)  , intent(inout)        :: frictionvel_inst

     ! locals
     integer                                        :: nclumps
     integer                                        :: nc
     type(bounds_type)                              :: bounds_clump
     ! locals
     real(r8) :: vol_ice
     real(r8) :: eff_porosity
     integer :: nlevsoil
     integer :: j, i
     integer :: s
     integer :: c, g
     integer :: p   ! HLM patch index
     integer :: ft  ! plant functional type
     logical  :: do_landuse_update        ! local flag to pass transitions update to fates

     ! Set the FATES global time and date variables
     call GetAndSetTime

     nclumps = get_proc_clumps()

     !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,s,c,j,vol_ice,eff_porosity)
     do nc = 1, nclumps

        if ( this%fates(nc)%nsites>0 ) then

           call get_clump_bounds(nc, bounds_clump)

           do s = 1,this%fates(nc)%nsites
              call init_site_vars(this%fates(nc)%sites(s),this%fates(nc)%bc_in(s),this%fates(nc)%bc_out(s) )
              call zero_site(this%fates(nc)%sites(s))
           end do

           call set_site_properties(this%fates(nc)%nsites, &
                                    this%fates(nc)%sites, &
                                    this%fates(nc)%bc_in)

           ! ----------------------------------------------------------------------------
           ! Initialize satellite phenology values if turned on
           ! ----------------------------------------------------------------------------
            if(use_fates_sp)then
               do s = 1,this%fates(nc)%nsites
                  c = this%f2hmap(nc)%fcolumn(s)
                  do ft = surfpft_lb,surfpft_ub !set of pfts in HLM
                     ! here we are mapping from P space in the HLM to FT space in the sp_input arrays.
                     p = ft + col_pp%pfti(c) ! for an FT of 1 we want to use
                     this%fates(nc)%bc_in(s)%hlm_sp_tlai(ft) = canopystate_inst%tlai_patch(p)
                     this%fates(nc)%bc_in(s)%hlm_sp_tsai(ft) = canopystate_inst%tsai_patch(p)
                     this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = canopystate_inst%htop_patch(p)
                     if(canopystate_inst%htop_patch(p).lt.1.0e-20)then ! zero htop causes inifinite/nans. This is
                        this%fates(nc)%bc_in(s)%hlm_sp_htop(ft) = 0.01_r8
                     endif
                  end do ! p
               end do ! c
            end if ! SP

           ! ----------------------------------------------------------------------------
           ! Initialize Hydraulics Code if turned on
           ! Called prior to init_patches(). Site level rhizosphere shells must
           ! be set prior to cohort initialization.
           ! ----------------------------------------------------------------------------
           if (use_fates_planthydro) then

              do s = 1,this%fates(nc)%nsites
                 c = this%f2hmap(nc)%fcolumn(s)
                 nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

                 this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil) = &
                      soilstate_inst%watsat_col(c,1:nlevsoil)

                 this%fates(nc)%bc_in(s)%watres_sisl(1:nlevsoil) = &
                      soilstate_inst%watmin_col(c,1:nlevsoil)

                 this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil) = &
                      soilstate_inst%sucsat_col(c,1:nlevsoil)

                 this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil) = &
                      soilstate_inst%bsw_col(c,1:nlevsoil)

                 this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil) = &
                      col_ws%h2osoi_liq(c,1:nlevsoil)

                 this%fates(nc)%bc_in(s)%hksat_sisl(1:nlevsoil) = &
                       soilstate_inst%hksat_col(c,1:nlevsoil)

                 do j = 1, nlevsoil
                    vol_ice = min(soilstate_inst%watsat_col(c,j), &
                          col_ws%h2osoi_ice(c,j)/(col_pp%dz(c,j)*denice))
                    eff_porosity = max(0.01_r8,soilstate_inst%watsat_col(c,j)-vol_ice)
                    this%fates(nc)%bc_in(s)%eff_porosity_sl(j) = eff_porosity
                 end do

              end do

              call HydrSiteColdStart(this%fates(nc)%sites,this%fates(nc)%bc_in)
           end if

           ! Transfer initial values to fates
           if (use_fates_luh) then
              do s = 1,this%fates(nc)%nsites
                 c = this%f2hmap(nc)%fcolumn(s)
                 g = col_pp%gridcell(c)

                 this%fates(nc)%bc_in(s)%hlm_luh_states = landuse_states(:,g)
                 this%fates(nc)%bc_in(s)%hlm_luh_state_names = landuse_state_varnames
                 this%fates(nc)%bc_in(s)%hlm_luh_transitions = landuse_transitions(:,g)
                 this%fates(nc)%bc_in(s)%hlm_luh_transition_names = landuse_transition_varnames

                 if (trim(fates_harvest_mode) == fates_harvest_luh_area .or. &
                     trim(fates_harvest_mode) == fates_harvest_luh_mass) then
                    this%fates(nc)%bc_in(s)%hlm_harvest_rates = landuse_harvest(:,g)
                    this%fates(nc)%bc_in(s)%hlm_harvest_catnames = landuse_harvest_varnames
                    this%fates(nc)%bc_in(s)%hlm_harvest_units = landuse_harvest_units
                 end if

              end do
           end if

           ! Initialize patches
           call init_patches(this%fates(nc)%nsites, this%fates(nc)%sites, &
                             this%fates(nc)%bc_in)


           do s = 1,this%fates(nc)%nsites

              c = this%f2hmap(nc)%fcolumn(s)
              this%fates(nc)%bc_in(s)%max_rooting_depth_index_col = this%fates(nc)%bc_in(s)%nlevdecomp


              call ed_update_site(this%fates(nc)%sites(s), &
                   this%fates(nc)%bc_in(s), &
                   this%fates(nc)%bc_out(s), &
                   is_restarting = .false.)

              ! This call sends internal fates variables into the
              ! output boundary condition structures. Note: this is called
              ! internally in fates dynamics as well.

              call FluxIntoLitterPools(this%fates(nc)%sites(s), &
                   this%fates(nc)%bc_in(s), &
                   this%fates(nc)%bc_out(s))
           end do

           ! ------------------------------------------------------------------------
           ! Update diagnostics of FATES ecosystem structure used in HLM.
           ! ------------------------------------------------------------------------
           call this%wrap_update_hlmfates_dyn(nc,bounds_clump, &
                canopystate_inst,frictionvel_inst, .false.)

           ! ------------------------------------------------------------------------
           ! Update history IO fields that depend on ecosystem dynamics
           ! ------------------------------------------------------------------------
           if(fates_history_dimlevel(2)>0) then
              call fates_hist%flush_hvars(nc,upfreq_in=group_dyna_simple)
              do s = 1,this%fates(nc)%nsites
                 call fates_hist%zero_site_hvars(this%fates(nc)%sites(s), &
                      upfreq_in=group_dyna_simple)
              end do
              if(fates_history_dimlevel(2)>1) then
                 call fates_hist%flush_hvars(nc,upfreq_in=group_dyna_complx)
                 do s = 1,this%fates(nc)%nsites
                    call fates_hist%zero_site_hvars(this%fates(nc)%sites(s), &
                         upfreq_in=group_dyna_complx)
                 end do
              end if
           end if
           
           call fates_hist%update_history_dyn( nc, &
                this%fates(nc)%nsites,                 &
                this%fates(nc)%sites, &
                this%fates(nc)%bc_in)

        end if
     end do
     !$OMP END PARALLEL DO

   end subroutine init_coldstart

   ! ======================================================================================

   subroutine wrap_sunfrac(this,bounds_clump,top_af_inst,canopystate_inst)

      ! ---------------------------------------------------------------------------------
      ! This interface function is a wrapper call on ED_SunShadeFracs. The only
      ! returned variable is a patch vector, fsun_patch, which describes the fraction
      ! of the canopy that is exposed to sun.
      ! ---------------------------------------------------------------------------------

      implicit none

      ! Input Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type)              , intent(in)    :: bounds_clump

      ! direct and diffuse downwelling radiation (W/m2)
      type(topounit_atmospheric_flux),intent(in)     :: top_af_inst

      ! Input/Output Arguments to CLM
      type(canopystate_type),intent(inout) :: canopystate_inst

      ! Local Variables
      integer  :: p                           ! global index of the host patch
      integer  :: g                           ! global index of the host gridcell
      integer  :: t                           ! global index of the host topounit
      integer  :: c                           ! global index of the host column

      integer  :: s                           ! FATES site index
      integer  :: ifp                         ! FATEs patch index
                                              ! this is the order increment of patch
                                              ! on the site
      integer  :: nc                          ! clump index

      type(fates_patch_type), pointer :: cpatch  ! c"urrent" patch  INTERF-TODO: SHOULD
                                              ! BE HIDDEN AS A FATES PRIVATE

      associate( forc_solad => top_af_inst%solad, &
                 forc_solai => top_af_inst%solai, &
                 fsun       => canopystate_inst%fsun_patch, &
                 laisun     => canopystate_inst%laisun_patch, &
                 laisha     => canopystate_inst%laisha_patch )

        nc = bounds_clump%clump_index
        ! -------------------------------------------------------------------------------
        ! Convert input BC's
        ! The sun-shade calculations are performed only on FATES patches
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           t = col_pp%topounit(c)
           g = col_pp%gridcell(c)

           do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
              this%fates(nc)%bc_in(s)%solad_parb(ifp,:) = forc_solad(t,:)
              this%fates(nc)%bc_in(s)%solai_parb(ifp,:) = forc_solai(t,:)
           end do
        end do

        ! -------------------------------------------------------------------------------
        ! Call FATES public function to calculate internal sun/shade structures
        ! as well as total patch sun/shade fraction output boundary condition
        ! -------------------------------------------------------------------------------

        call FatesSunShadeFracs(this%fates(nc)%nsites, &
             this%fates(nc)%sites,  &
             this%fates(nc)%bc_in,  &
             this%fates(nc)%bc_out)

        ! -------------------------------------------------------------------------------
        ! Transfer the FATES output boundary condition for canopy sun/shade fraction
        ! to the HLM
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
              p = ifp+col_pp%pfti(c)
              fsun(p)   = this%fates(nc)%bc_out(s)%fsun_pa(ifp)
              laisun(p) = this%fates(nc)%bc_out(s)%laisun_pa(ifp)
              laisha(p) = this%fates(nc)%bc_out(s)%laisha_pa(ifp)
           end do
        end do

      end associate

   end subroutine wrap_sunfrac

   ! ===================================================================================

   subroutine prep_canopyfluxes(this, bounds_clump )

     ! ----------------------------------------------------------------------
     ! the main function for calculating photosynthesis is called within a
     ! loop based on convergence.  Some intitializations, including
     ! canopy resistance must be intitialized before the loop
     ! The photosyns_ structure is currently unused, leaving it for now
     ! in case we want to do any value initializing in future.
     ! ----------------------------------------------------------------------

     ! Arguments
     class(hlm_fates_interface_type), intent(inout) :: this
     type(bounds_type)              , intent(in)    :: bounds_clump

     ! locals
     integer                                        :: c,s
     integer                                        :: nc

     nc = bounds_clump%clump_index
     do s = 1, this%fates(nc)%nsites
        ! filter flag == 1 means that this patch has not been called for photosynthesis
        this%fates(nc)%bc_in(s)%filter_photo_pa(:) = 1

        ! set transpiration input boundary condition to zero. The exposed
        ! vegetation filter may not even call every patch.
        if (use_fates_planthydro) then
           this%fates(nc)%bc_in(s)%qflx_transp_pa(:) = 0._r8
        end if

     end do
  end subroutine prep_canopyfluxes

   ! ====================================================================================

   subroutine wrap_btran(this,bounds_clump,fn,filterc,soilstate_inst, &
                         energyflux_inst,  &
                         soil_water_retention_curve)

      ! ---------------------------------------------------------------------------------
      ! This subroutine calculates btran for FATES, this will be an input boundary
      ! condition for FATES photosynthesis/transpiration.
      !
      ! This subroutine also calculates rootr
      !
      ! ---------------------------------------------------------------------------------

      

      implicit none

      ! Arguments
      class(hlm_fates_interface_type), intent(inout) :: this
      type(bounds_type)              , intent(in)    :: bounds_clump
      integer                , intent(in)            :: fn
      integer                , intent(in)            :: filterc(fn) ! This is a list of
                                                                        ! columns with exposed veg
      type(soilstate_type)   , intent(inout)         :: soilstate_inst
      type(energyflux_type)  , intent(inout)         :: energyflux_inst
      class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

      ! local variables
      real(r8) :: smp_node ! Soil suction potential, negative, [mm]
      real(r8) :: s_node
      integer  :: s
      integer  :: c
      integer  :: j
      integer  :: ifp
      integer  :: p
      integer  :: nlevsoil
      integer  :: nc

      associate(&
         sucsat      => soilstate_inst%sucsat_col           , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
         watsat      => soilstate_inst%watsat_col           , & ! Input:  [real(r8) (:,:) ]  volumetric soil water at saturation (porosity)
         bsw         => soilstate_inst%bsw_col              , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
         eff_porosity => soilstate_inst%eff_porosity_col    , & ! Input:  [real(r8) (:,:) ]  effective porosity = porosity - vol_ice
         t_soisno     => col_es%t_soisno                    , & ! Input:  [real(r8) (:,:) ]  soil temperature (Kelvin)
         h2osoi_liqvol => col_ws%h2osoi_liqvol , & ! Input: [real(r8) (:,:) ]  liquid volumetric moisture, will be used for BeTR
         btran       => energyflux_inst%btran_patch         , & ! Output: [real(r8) (:)   ]  transpiration wetness factor (0 to 1)
         btran2       => energyflux_inst%btran2_patch       , & ! Output: [real(r8) (:)   ]
         rresis      => energyflux_inst%rresis_patch        , & ! Output: [real(r8) (:,:) ]  root resistance by layer (0-1)  (nlevgrnd)
         rootr       => soilstate_inst%rootr_patch          & ! Output: [real(r8) (:,:) ]  Fraction of water uptake in each layer
         )


        nc = bounds_clump%clump_index

        ! -------------------------------------------------------------------------------
        ! Convert input BC's
        ! Critical step: a filter is being passed in that dictates which columns have
        ! exposed vegetation (above snow).  This is necessary, because various hydrologic
        ! variables like h2osoi_liqvol are not calculated and will have uninitialized
        ! values outside this list.
        !
        ! bc_in(s)%filter_btran      (this is in, but is also used in this subroutine)
        !
        ! We also filter a second time within this list by determining which soil layers
        ! have conditions for active uptake based on soil moisture and temperature. This
        ! must be determined by FATES (science stuff).  But the list of layers and patches
        ! needs to be passed back to the interface, because it then needs to request
        ! suction on these layers via CLM/ALM functions.  We cannot wide-swath calculate
        ! this on all layers, because values with no moisture or low temps will generate
        ! unstable values and cause sigtraps.
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

           ! Check to see if this column is in the exposed veg filter
           if( any(filterc==c) )then

              this%fates(nc)%bc_in(s)%filter_btran = .true.
              do j = 1,nlevsoil
                 this%fates(nc)%bc_in(s)%tempk_sl(j)         = t_soisno(c,j)
                 this%fates(nc)%bc_in(s)%h2o_liqvol_sl(j)    = h2osoi_liqvol(c,j)
                 this%fates(nc)%bc_in(s)%eff_porosity_sl(j)  = eff_porosity(c,j)
                 this%fates(nc)%bc_in(s)%watsat_sl(j)        = watsat(c,j)
              end do

           else
              this%fates(nc)%bc_in(s)%filter_btran = .false.
              this%fates(nc)%bc_in(s)%tempk_sl(:)         = -999._r8
              this%fates(nc)%bc_in(s)%h2o_liqvol_sl(:)    = -999._r8
              this%fates(nc)%bc_in(s)%eff_porosity_sl(:)  = -999._r8
              this%fates(nc)%bc_in(s)%watsat_sl(:)        = -999._r8
           end if

        end do

        ! -------------------------------------------------------------------------------
        ! This function evaluates the ground layer to determine if
        ! root water uptake can happen, and soil suction should even
        ! be calculated.  We ask FATES for a boundary condition output
        ! logical because we don't want science calculations in the interface
        ! yet... hydrology (suction calculation) is provided by the host
        ! so we need fates to tell us where to calculate suction
        ! but not calculate it itself. Yeah, complicated, but thats life.
        ! -------------------------------------------------------------------------------
        call get_active_suction_layers(this%fates(nc)%nsites, &
             this%fates(nc)%sites,  &
             this%fates(nc)%bc_in,  &
             this%fates(nc)%bc_out)

        ! Now that the active layers of water uptake have been decided by fates
        ! Calculate the suction that is passed back to fates
        ! Note that the filter_btran is unioned with active_suction_sl

        do s = 1, this%fates(nc)%nsites
           c = this%f2hmap(nc)%fcolumn(s)
           nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

           do j = 1,nlevsoil
              if(this%fates(nc)%bc_out(s)%active_suction_sl(j)) then
                 s_node = max(h2osoi_liqvol(c,j)/eff_porosity(c,j),0.01_r8)
                 call soil_water_retention_curve%soil_suction( soilstate_inst%sucsat_col(c,j), &
                       s_node, &
                       soilstate_inst%bsw_col(c,j), &
                       smp_node)

                 ! Non-fates places a maximum (which is a negative upper bound) on smp

                 this%fates(nc)%bc_in(s)%smp_sl(j)           = smp_node
              end if
           end do
        end do

        ! -------------------------------------------------------------------------------
        ! Suction and active uptake layers calculated, lets calculate uptake (btran)
        ! This will calculate internals, as well as output boundary conditions:
        ! btran, rootr
        ! -------------------------------------------------------------------------------

        call btran_ed(this%fates(nc)%nsites, &
             this%fates(nc)%sites,  &
             this%fates(nc)%bc_in,  &
             this%fates(nc)%bc_out)

        ! -------------------------------------------------------------------------------
        ! Convert output BC's
        ! For CLM/ALM this wrapper provides return variables that should
        ! be similar to that of calc_root_moist_stress().  However,
        ! CLM/ALM-FATES simulations will no make use of rresis, btran or btran2
        ! outside of FATES. We do not have code in place to calculate btran2 or
        ! rresis right now, so we force to bad.  We have btran calculated so we
        ! pass it in case people want diagnostics.  rootr is actually the only
        ! variable that will be used, as it is needed to help distribute the
        ! the transpiration sink to the appropriate layers. (RGK)
        ! -------------------------------------------------------------------------------

        do s = 1, this%fates(nc)%nsites
           nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil
           c = this%f2hmap(nc)%fcolumn(s)
           do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno

              p = ifp+col_pp%pfti(c)

              do j = 1,nlevsoil

                 rresis(p,j) = -999.9  ! We do not calculate this correctly
                 ! it should not thought of as valid output until we decide to.
                 rootr(p,j)  = this%fates(nc)%bc_out(s)%rootr_pasl(ifp,j)
                 btran(p)    = this%fates(nc)%bc_out(s)%btran_pa(ifp)
                 btran2(p)   = -999.9  ! Not available, force to nonsense

              end do
           end do
        end do
      end associate

   end subroutine wrap_btran

   ! ====================================================================================

   subroutine wrap_photosynthesis(this, bounds_clump, fn, filterp, &
         esat_tv, eair, oair, cair, rb, dayl_factor,             &
         atm2lnd_inst, canopystate_inst, photosyns_inst)

    use shr_log_mod       , only : errMsg => shr_log_errMsg
    use abortutils        , only : endrun
    use decompMod         , only : bounds_type
    use elm_varcon        , only : rgas, tfrz, namep
    use elm_varctl        , only : iulog
    use quadraticMod      , only : quadratic

    !
    ! !ARGUMENTS:
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type)      , intent(in)            :: bounds_clump
    integer                , intent(in)            :: fn                          ! size of pft filter
    integer                , intent(in)            :: filterp(fn)                 ! pft filter
    real(r8)               , intent(in)            :: esat_tv(bounds_clump%begp: )      ! saturation vapor pressure at t_veg (Pa)
    real(r8)               , intent(in)            :: eair( bounds_clump%begp: )        ! vapor pressure of canopy air (Pa)
    real(r8)               , intent(in)            :: oair( bounds_clump%begp: )        ! Atmospheric O2 partial pressure (Pa)
    real(r8)               , intent(in)            :: cair( bounds_clump%begp: )        ! Atmospheric CO2 partial pressure (Pa)
    real(r8)               , intent(in)            :: rb( bounds_clump%begp: )          ! boundary layer resistance (s/m)
    real(r8)               , intent(in)            :: dayl_factor( bounds_clump%begp: ) ! scalar (0-1) for daylength
    type(atm2lnd_type)     , intent(in)            :: atm2lnd_inst
    type(canopystate_type) , intent(inout)         :: canopystate_inst
    type(photosyns_type)   , intent(inout)         :: photosyns_inst

    integer                                        :: nlevsoil
    integer                                        :: s,t,c,p,ifp,j,icp,nc
    real(r8)                                       :: dtime

    call t_startf('edpsn')
    associate(&
          t_soisno  => col_es%t_soisno               , &
          t_veg     => veg_es%t_veg                  , &
          tgcm      => veg_es%thm                    , &
          forc_pbot => top_as%pbot                   , &
          rssun     => photosyns_inst%rssun_patch    , &
          rssha     => photosyns_inst%rssha_patch    , &
          psnsun    => photosyns_inst%psnsun_patch   , &
          psnsha    => photosyns_inst%psnsha_patch)


      nc = bounds_clump%clump_index

      do s = 1, this%fates(nc)%nsites

         c = this%f2hmap(nc)%fcolumn(s)
         t = col_pp%topounit(c)

         nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

         do j = 1,nlevsoil
            this%fates(nc)%bc_in(s)%t_soisno_sl(j)   = t_soisno(c,j)  ! soil temperature (Kelvin)
         end do
         this%fates(nc)%bc_in(s)%forc_pbot           = forc_pbot(t)   ! atmospheric pressure (Pa)

         do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno

            p = ifp+col_pp%pfti(c)

            ! Check to see if this patch is in the filter
            ! Note that this filter is most likely changing size, and getting smaller
            ! and smaller as more patch have converged on solution
            if( any(filterp==p) )then

               ! This filter is flushed to 1 before the canopyflux stability iterator
               ! It is set to status 2 if it is an active patch within the iterative loop
               ! After photosynthesis is called, it is upgraded to 3 if it was called.
               ! After all iterations we can evaluate which patches have a final flag
               ! of 3 to check if we missed any.

               this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) = 2
               this%fates(nc)%bc_in(s)%dayl_factor_pa(ifp) = dayl_factor(p) ! scalar (0-1) for daylength
               this%fates(nc)%bc_in(s)%esat_tv_pa(ifp)     = esat_tv(p)     ! saturation vapor pressure at t_veg (Pa)
               this%fates(nc)%bc_in(s)%eair_pa(ifp)        = eair(p)        ! vapor pressure of canopy air (Pa)
               this%fates(nc)%bc_in(s)%oair_pa(ifp)        = oair(p)        ! Atmospheric O2 partial pressure (Pa)
               this%fates(nc)%bc_in(s)%cair_pa(ifp)        = cair(p)        ! Atmospheric CO2 partial pressure (Pa)
               this%fates(nc)%bc_in(s)%rb_pa(ifp)          = rb(p)          ! boundary layer resistance (s/m)
               this%fates(nc)%bc_in(s)%t_veg_pa(ifp)       = t_veg(p)       ! vegetation temperature (Kelvin)
               this%fates(nc)%bc_in(s)%tgcm_pa(ifp)        = tgcm(p)        ! air temperature at agcm
                                                                            ! reference height (kelvin)
            end if
         end do
      end do

      dtime = real(get_step_size(),r8)

      ! Call photosynthesis

      call FatesPlantRespPhotosynthDrive (this%fates(nc)%nsites, &
                                this%fates(nc)%sites,  &
                                this%fates(nc)%bc_in,  &
                                this%fates(nc)%bc_out, &
                                dtime)

      ! Perform a double check to see if all patches on naturally vegetated columns
      ! were activated for photosynthesis
      ! ---------------------------------------------------------------------------------
      do icp = 1,fn
         p = filterp(icp)
         c = veg_pp%column(p)
         s = this%f2hmap(nc)%hsites(c)
         ! do if structure here and only pass natveg columns
         ifp = p-col_pp%pfti(c)
         if(this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) /= 2)then
            write(iulog,*) 'Not all patches on the natveg column in the photosynthesis'
            write(iulog,*) 'filter ran photosynthesis: s, p, icp, ifp: ',s,p,icp,ifp
            call endrun(msg=errMsg(sourcefile, __LINE__))
         else
            this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) = 3
            rssun(p) = this%fates(nc)%bc_out(s)%rssun_pa(ifp)
            rssha(p) = this%fates(nc)%bc_out(s)%rssha_pa(ifp)

            ! These fields are marked with a bad-value flag
            photosyns_inst%psnsun_patch(p)   = spval
            photosyns_inst%psnsha_patch(p)   = spval
         end if
      end do

    end associate
    call t_stopf('edpsn')

 end subroutine wrap_photosynthesis

 ! ======================================================================================

 subroutine wrap_accumulatefluxes(this, bounds_clump, fn, filterp)

   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type)              , intent(in)    :: bounds_clump
   integer                        , intent(in)    :: fn                   ! size of pft filter
   integer                        , intent(in)    :: filterp(fn)          ! pft filter

   ! Locals
   integer                                        :: s,c,p,ifp,icp
   real(r8)                                       :: dtime
   integer                                        :: nc

   nc = bounds_clump%clump_index
    ! Run a check on the filter
    do icp = 1,fn
       p = filterp(icp)
       c = veg_pp%column(p)
       s = this%f2hmap(nc)%hsites(c)
       ifp = p-col_pp%pfti(c)
       if(this%fates(nc)%bc_in(s)%filter_photo_pa(ifp) /= 3)then
          write(iulog,*) 'Not all patches on the natveg column in the canopys'
          write(iulog,*) 'filter ran canopy fluxes: s, p, icp, ifp: ',s,p,icp,ifp
          call endrun(msg=errMsg(sourcefile, __LINE__))
       end if
    end do


    dtime = real(get_step_size(),r8)

    call  AccumulateFluxes_ED(this%fates(nc)%nsites,  &
                               this%fates(nc)%sites, &
                               this%fates(nc)%bc_in,  &
                               this%fates(nc)%bc_out, &
                               dtime)
    return
 end subroutine wrap_accumulatefluxes

 ! ======================================================================================

 subroutine wrap_WoodProducts(this, bounds_clump, fc, filterc)

   use FatesConstantsMod     , only : g_per_kg

   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type)              , intent(in)    :: bounds_clump
   integer                        , intent(in)    :: fc                   ! size of column filter
   integer                        , intent(in)    :: filterc(fc)          ! column filter
   
   ! Locacs
   integer                                        :: s,c,icc
   integer                                        :: nc

   associate(&
         gpp     => col_cf%gpp    , &
         ar     => col_cf%ar    , &
         hrv_deadstemc_to_prod10c     => col_cf%hrv_deadstemc_to_prod10c    , &
         hrv_deadstemc_to_prod100c    => col_cf%hrv_deadstemc_to_prod100c)
 
    nc = bounds_clump%clump_index
    ! Loop over columns
    do icc = 1,fc
       c = filterc(icc)
       s = this%f2hmap(nc)%hsites(c)

       ! Pass harvested wood products to ELM variable
       hrv_deadstemc_to_prod10c(c)  = this%fates(nc)%bc_out(s)%hrv_deadstemc_to_prod10c
       hrv_deadstemc_to_prod100c(c) = this%fates(nc)%bc_out(s)%hrv_deadstemc_to_prod100c

       ! Pass LUC related C fluxes which are calculated in FATES [gC m-2 s-1]
       gpp(c) = this%fates(nc)%bc_out(s)%gpp_site*g_per_kg
       ar(c) = this%fates(nc)%bc_out(s)%ar_site*g_per_kg

    end do

    end associate
    return
 end subroutine wrap_WoodProducts

 ! ======================================================================================
 
 subroutine wrap_canopy_radiation(this, bounds_clump, &
         num_vegsol, filter_vegsol, coszen, surfalb_inst)


    ! Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),  intent(in)             :: bounds_clump
    ! filter for vegetated pfts with coszen>0
    integer            , intent(in)            :: num_vegsol
    integer            , intent(in)            :: filter_vegsol(num_vegsol)
    ! cosine solar zenith angle for next time step
    real(r8)           , intent(in)            :: coszen( bounds_clump%begp: )
    type(surfalb_type) , intent(inout)         :: surfalb_inst

    ! locals
    integer                                    :: s,c,p,ifp,icp,nc

    associate(&
         albgrd_col   =>    surfalb_inst%albgrd_col         , & !in
         albgri_col   =>    surfalb_inst%albgri_col         , & !in
         albd         =>    surfalb_inst%albd_patch         , & !out
         albi         =>    surfalb_inst%albi_patch         , & !out
         fabd         =>    surfalb_inst%fabd_patch         , & !out
         fabi         =>    surfalb_inst%fabi_patch         , & !out
         ftdd         =>    surfalb_inst%ftdd_patch         , & !out
         ftid         =>    surfalb_inst%ftid_patch         , & !out
         ftii         =>    surfalb_inst%ftii_patch)            !out

    nc = bounds_clump%clump_index

    do s = 1, this%fates(nc)%nsites

       c = this%f2hmap(nc)%fcolumn(s)
       do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno

          p = ifp+col_pp%pfti(c)

          if( any(filter_vegsol==p) )then

             this%fates(nc)%bc_in(s)%filter_vegzen_pa(ifp) = .true.
             this%fates(nc)%bc_in(s)%coszen_pa(ifp)  = coszen(p)
             this%fates(nc)%bc_in(s)%albgr_dir_rb(:) = albgrd_col(c,:)
             this%fates(nc)%bc_in(s)%albgr_dif_rb(:) = albgri_col(c,:)

             if (veg_es%t_veg(p) <= tfrz) then
                this%fates(nc)%bc_in(s)%fcansno_pa(ifp) = veg_ws%fwet(p)
             else
                this%fates(nc)%bc_in(s)%fcansno_pa(ifp) = 0._r8
             end if
             
          else

             this%fates(nc)%bc_in(s)%filter_vegzen_pa(ifp) = .false.

          end if

       end do
    end do

    call FatesNormalizedCanopyRadiation(this%fates(nc)%nsites,  &
         this%fates(nc)%sites, &
         this%fates(nc)%bc_in,  &
         this%fates(nc)%bc_out)

    ! Pass FATES BC's back to HLM
    ! -----------------------------------------------------------------------------------
    do icp = 1,num_vegsol
       p = filter_vegsol(icp)
       c = veg_pp%column(p)
       s = this%f2hmap(nc)%hsites(c)
       ! do if structure here and only pass natveg columns
       ifp = p-col_pp%pfti(c)

       if(.not.this%fates(nc)%bc_in(s)%filter_vegzen_pa(ifp) )then
          write(iulog,*) 's,p,ifp',s,p,ifp
          write(iulog,*) 'Not all patches on the natveg column were passed to canrad',veg_pp%sp_pftorder_index(p)
          call endrun(msg=errMsg(sourcefile, __LINE__))
       else
          albd(p,:) = this%fates(nc)%bc_out(s)%albd_parb(ifp,:)
          albi(p,:) = this%fates(nc)%bc_out(s)%albi_parb(ifp,:)
          fabd(p,:) = this%fates(nc)%bc_out(s)%fabd_parb(ifp,:)
          fabi(p,:) = this%fates(nc)%bc_out(s)%fabi_parb(ifp,:)
          ftdd(p,:) = this%fates(nc)%bc_out(s)%ftdd_parb(ifp,:)
          ftid(p,:) = this%fates(nc)%bc_out(s)%ftid_parb(ifp,:)
          ftii(p,:) = this%fates(nc)%bc_out(s)%ftii_parb(ifp,:)
       end if
    end do

  end associate

 end subroutine wrap_canopy_radiation

 ! ======================================================================================
 
 subroutine WrapGlobalSeedDispersal(this,is_restart_flag)
   
   ! Call mpi procedure to provide the global seed output distribution array to every gridcell.
   ! This could be conducted with a more sophisticated halo-type structure or distributed graph.
   
   use decompMod,              only : procinfo
   use spmdMod,                only : MPI_REAL8, mpicom
   use FatesDispersalMod,      only : lneighbors, neighbor_type
   use FatesInterfaceTypesMod, only : numpft_fates => numpft
   
   ! Arguments
   class(hlm_fates_interface_type), intent(inout) :: this
   logical, optional                              :: is_restart_flag

   ! Local
   integer :: numg ! total number of gridcells across all processors
   integer :: ier  ! error code
   integer :: g    ! gridcell index

   logical :: set_restart_flag ! local logical variable to pass to IsItDispersalTime

#ifdef _OPENMP
   logical, external :: omp_in_parallel
#endif
       
   type (neighbor_type), pointer :: neighbor

   ! Check to see if we are not in a threaded region.  Fail the run if this returns true.
#ifdef _OPENMP
   if (omp_in_parallel()) then
      call endrun(msg='elmfates interface error: MPI routine called within threaded region'//&
           errMsg(sourcefile, __LINE__))
   end if
#endif

   ! This should only be run once per day
   if(is_beg_curr_day()) then

      ! If WrapGlobalSeedDispersal is being called at the end a fates restart call,
      ! pass .false. to the set_dispersed_flag to avoid updating the
      ! global dispersal date
      set_restart_flag = .true.
      if (present(is_restart_flag)) then
         if (is_restart_flag) set_restart_flag = .false.
      end if

      call t_startf('fates-seed-mpi_allgatherv')

      if (IsItDispersalTime(setdispersedflag=set_restart_flag)) then

         ! Re-initialize incoming seed buffer for this time step
         this%fates_seed%incoming_global(:,:) = 0._r8
         this%fates_seed%outgoing_global(:,:) = 0._r8

         ! Distribute and sum outgoing seed data from all nodes to all nodes
         call MPI_Allgatherv(this%fates_seed%outgoing_local, procinfo%ncells*numpft_fates, MPI_REAL8, &
              this%fates_seed%outgoing_global, this%fates_seed%ncells_array*numpft_fates, this%fates_seed%begg_array, &
              MPI_REAL8, mpicom, ier)

        if (ier /= 0) then
            call endrun(msg='elmfates interface error: MPI_Allgatherv failed'//&
                 errMsg(sourcefile, __LINE__))
         end if

         ! zero outgoing local for all gridcells outside threaded region now that we've passed them out
         this%fates_seed%outgoing_local(:,:) = 0._r8

         ! Calculate the current gridcell incoming seed for each gridcell index
         ! This should be conducted outside of a threaded region to provide access to
         ! the neighbor%gindex which might not be available in via the clumped index
         call get_proc_global(ng=numg)
         do g = 1, numg
            neighbor => lneighbors(g)%first_neighbor
            do while (associated(neighbor))
               ! This also applies the same neighborhood distribution scheme to all pfts
               ! This needs to have a per pft density probability value
               this%fates_seed%incoming_global(g,:) = this%fates_seed%incoming_global(g,:) + &
                    this%fates_seed%outgoing_global(neighbor%gindex,:) * &
                    neighbor%density_prob(:) / lneighbors(g)%neighbor_count
               neighbor => neighbor%next_neighbor
            end do
         end do
      endif
   endif
   call t_stopf('fates-seed-mpi_allgatherv')
   
 end subroutine WrapGlobalSeedDispersal
 
 ! ======================================================================================

 subroutine WrapUpdateFatesSeedInOut(this,bounds_clump)

    ! This subroutine pass seed_id_global to bc_in and reset seed_out
   
    ! Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),  intent(in)                 :: bounds_clump
    
    integer  :: g                           ! global index of the host gridcell
    integer  :: c                           ! global index of the host column
    integer  :: s                           ! FATES site index
    integer  :: nc                          ! clump index

    call t_startf('fates-seed-disperse')

    nc = bounds_clump%clump_index

    ! Check that it is the beginning of the current dispersal time step
    if (IsItDispersalTime()) then
       do s = 1, this%fates(nc)%nsites
          c = this%f2hmap(nc)%fcolumn(s)
          g = col_pp%gridcell(c)

          ! assuming equal area for all sites, seed_id_global in [kg/grid/day], seed_in in [kg/site/day]
          this%fates(nc)%sites(s)%seed_in(:) = this%fates_seed%incoming_global(g,:)
          this%fates(nc)%sites(s)%seed_out(:) = 0._r8  ! reset seed_out
       end do
    else
       ! if it is not the dispersing time, pass in zero
       this%fates(nc)%sites(s)%seed_in(:) = 0._r8
    end if

    call t_stopf('fates-seed-disperse')

 end subroutine WrapUpdateFatesSeedInOut

 ! ======================================================================================

 subroutine wrap_update_hifrq_hist(this, bounds_clump )



    ! Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),  intent(in)                 :: bounds_clump

    ! locals
    real(r8) :: dtime
    integer  :: nstep
    logical  :: is_beg_day
    integer  :: s,c,nc

    associate(&
        hr            => col_cf%hr,      & ! (gC/m2/s) total heterotrophic respiration
        totsomc       => col_cs%totsomc, & ! (gC/m2) total soil organic matter carbon
        totlitc       => col_cs%totlitc)   ! (gC/m2) total litter carbon in BGC pools

      nc = bounds_clump%clump_index
      dtime = real(get_step_size(),r8)

      ! Summarize Net Fluxes
      do s = 1, this%fates(nc)%nsites
         c = this%f2hmap(nc)%fcolumn(s)
         this%fates(nc)%bc_in(s)%tot_het_resp = hr(c)
         this%fates(nc)%bc_in(s)%tot_somc     = totsomc(c)
         this%fates(nc)%bc_in(s)%tot_litc     = totlitc(c)
      end do

      ! Update history variables that track these variables
      call fates_hist%update_history_hifrq(nc, &
           this%fates(nc)%nsites,  &
           this%fates(nc)%sites,   &
           this%fates(nc)%bc_in,   &
           this%fates(nc)%bc_out,  &
           dtime)


    end associate
end subroutine wrap_update_hifrq_hist

 ! ======================================================================================


 subroutine TransferZ0mDisp(this,bounds_clump,frictionvel_inst,canopystate_inst)

    ! Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                   :: bounds_clump
    type(canopystate_type)  , intent(inout)        :: canopystate_inst
    type(frictionvel_type)  , intent(inout)        :: frictionvel_inst

    ! Locals
    integer :: ci   ! Current clump index
    integer :: s    ! Site index
    integer :: c    ! Column index
    integer :: ifp  ! Fates patch index
    integer :: p    ! CLM patch index

    ci = bounds_clump%clump_index

    do s = 1, this%fates(ci)%nsites
       c = this%f2hmap(ci)%fcolumn(s)

       frictionvel_inst%z0m_patch(col_pp%pfti(c)+1:col_pp%pftf(c)) = 0.0_r8
       canopystate_inst%displa_patch(col_pp%pfti(c)+1:col_pp%pftf(c)) = 0.0_r8

       do ifp = 1, this%fates(ci)%sites(s)%youngest_patch%patchno
          p = ifp+col_pp%pfti(c)
          frictionvel_inst%z0m_patch(p) = this%fates(ci)%bc_out(s)%z0m_pa(ifp)
          canopystate_inst%displa_patch(p) = this%fates(ci)%bc_out(s)%displa_pa(ifp)
       end do
    end do

    return
 end subroutine TransferZ0mDisp

  !-----------------------------------------------------------------------

 subroutine InterpFileInputs(this, bounds)
   !
   ! !DESCRIPTION:
   ! Interpolate inputs from files
   !
   ! NOTE(wjs, 2016-02-23) Stuff done here could probably be done at the end of
   ! InitEachTimeStep, rather than in this separate routine, except for the
   ! fact that
   ! (currently) this Interp stuff is done with proc bounds rather thna clump
   ! bounds. I
   ! think that is needed so that you don't update a given stream multiple
   ! times. If we
   ! rework the handling of threading / clumps so that there is a separate
   ! object for
   ! each clump, then I think this problem would disappear - at which point we
   ! could
   ! remove this Interp routine, moving its body to the end of
   ! InitEachTimeStep.
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type), intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:

   character(len=*), parameter :: subname = 'InterpFileInputs'
   !-----------------------------------------------------------------------

   call t_startf('fates_interpfileinputs')

   call this%fates_fire_data_method%FireInterp(bounds)

   call t_stopf('fates_interpfileinputs')

 end subroutine InterpFileInputs

 ! ======================================================================================
 
 subroutine Init2(this, bounds, NLFilename)
   !
   ! !DESCRIPTION:
   ! Initialization after subgrid weights are determined
   !
   ! This copy should only be called if use_fates is .true.
   !
   ! !USES:
   ! use CNStateType, only : cnstate_type
   !
   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type), intent(in)                   :: bounds
   character(len=*), intent(in) :: NLFilename ! Namelist filename
   ! type(cnstate_type), intent(in)                 :: cnstate_vars

   !
   ! !LOCAL VARIABLES:

   character(len=*), parameter :: subname = 'Init2'
   !-----------------------------------------------------------------------

   call t_startf('fates_init2')

   write(iulog,*) 'Init2: calling FireInit'
   call this%fates_fire_data_method%FireInit(bounds, NLFilename)

   call t_stopf('fates_init2')

 end subroutine Init2

 ! ======================================================================================
 
 subroutine InitAccBuffer(this, bounds)
   !
   ! !DESCRIPTION:
   ! Initialized any accumulation buffers needed for FATES
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type), intent(in)                  :: bounds
   !
   ! !LOCAL VARIABLES:

   character(len=*), parameter :: subname = 'InitAccBuffer'
   !-----------------------------------------------------------------------

   call t_startf('fates_initaccbuff')

   call this%fates_fire_data_method%InitAccBuffer( bounds )

   call t_stopf('fates_initaccbuff')

 end subroutine InitAccBuffer

! ======================================================================================
 
 subroutine InitAccVars(this, bounds)
   !
   ! !DESCRIPTION:
   ! Initialized any accumulation variables needed for FATES
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type), intent(in) :: bounds
   !
   ! !LOCAL VARIABLES:

   character(len=*), parameter :: subname = 'InitAccVars'
   !-----------------------------------------------------------------------

   call t_startf('fates_initaccvars')

   call this%fates_fire_data_method%InitAccVars( bounds )

   call t_stopf('fates_initaccvars')

 end subroutine InitAccVars

 ! ======================================================================================
 
 subroutine UpdateAccVars(this, bounds_proc)
  
   ! !DESCRIPTION:
   ! Update any accumulation variables needed for FATES
   !
   ! !USES:
   !
   ! !ARGUMENTS:
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type), intent(in)                  :: bounds_proc
   !
 
   character(len=*), parameter :: subname = 'UpdateAccVars'
   !-----------------------------------------------------------------------

   call t_startf('fates_updateaccvars')

   call this%fates_fire_data_method%UpdateAccVars( bounds_proc )
   
   call t_stopf('fates_updateaccvars')

 end subroutine UpdateAccVars
 
 ! ======================================================================================

 subroutine WrapUpdateFatesRmean(this, nc)
   
   class(hlm_fates_interface_type), intent(inout) :: this
   integer,intent(in) :: nc
   
   ! !LOCAL VARIABLES:
   integer                     :: s,c,p,ifp  ! indices and loop counters
   
   do s = 1, this%fates(nc)%nsites
      c = this%f2hmap(nc)%fcolumn(s)
      do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
         p = ifp+col_pp%pfti(c)
         this%fates(nc)%bc_in(s)%t_veg_pa(ifp) = veg_es%t_veg(p)
      end do
   end do

   call UpdateFatesRMeansTStep(this%fates(nc)%sites,this%fates(nc)%bc_in,this%fates(nc)%bc_out)
   
  end subroutine WrapUpdateFatesRmean

 ! ======================================================================================
 
 subroutine init_history_io(this,bounds_proc)

   use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp

   use FatesConstantsMod, only : fates_short_string_length, fates_long_string_length
   use FatesIOVariableKindMod, only : site_r8, site_soil_r8, site_size_pft_r8
   use FatesIOVariableKindMod, only : site_size_r8, site_pft_r8, site_age_r8
   use FatesIOVariableKindMod, only : site_fuel_r8, site_cwdsc_r8, site_scag_r8
   use FatesIOVariableKindMod, only : site_scagpft_r8, site_agepft_r8, site_agefuel_r8
   use FatesIOVariableKindMod, only : site_height_r8, site_elem_r8, site_elpft_r8
   use FatesIOVariableKindMod, only : site_elcwd_r8, site_elage_r8
   use FatesIOVariableKindMod, only : site_coage_r8, site_coage_pft_r8
   use FatesIOVariableKindMod, only : site_can_r8, site_cnlf_r8, site_cnlfpft_r8
   use FatesIOVariableKindMod, only : site_cdpf_r8, site_cdsc_r8, site_cdam_r8
   use FatesIOVariableKindMod, only : site_landuse_r8, site_lulu_r8
   use FatesIODimensionsMod, only : fates_bounds_type


   ! Arguments
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type),intent(in)                   :: bounds_proc  ! Currently "proc"


   ! Locals
   type(bounds_type)                              :: bounds_clump
   integer :: nvar  ! number of IO variables found
   integer :: ivar  ! variable index 1:nvar
   integer :: nc    ! thread counter 1:nclumps
   integer :: nclumps ! number of threads on this proc
   integer :: s     ! FATES site index
   integer :: c     ! ALM/CLM column index
   character(len=fates_short_string_length) :: dim2name
   character(len=fates_long_string_length) :: ioname
   integer :: d_index, dk_index

   type(fates_bounds_type) :: fates_bounds
   type(fates_bounds_type) :: fates_clump

   ! This routine initializes the types of output variables
   ! not the variables themselves, just the types
   ! ---------------------------------------------------------------------------------

   nclumps = get_proc_clumps()

   ! ------------------------------------------------------------------------------------
   ! PART I: Set FATES DIMENSIONING INFORMATION
   !
   ! -------------------------------------------------------------------------------
   ! Those who wish add variables that require new dimensions, please
   ! see FATES: FatesHistoryInterfaceMod.F90.  Dimension types are defined at the top of the
   ! module, and a new explicitly named instance of that type should be created.
   ! With this new dimension, a new output type/kind can contain that dimension.
   ! A new type/kind can be added to the dim_kinds structure, which defines its members
   ! in created in init_dim_kinds_maps().  Make sure to increase the size of fates_num_dim_kinds.
   ! A type/kind of output is defined by the data type (ie r8,int,..)
   ! and the dimensions.  Keep in mind that 3D variables (or 4D if you include time)
   ! are not really supported in CLM/ALM right now.  There are ways around this
   ! limitations by creating combined dimensions, for instance the size+pft dimension
   ! "scpf"
   ! ------------------------------------------------------------------------------------

   call hlm_bounds_to_fates_bounds(bounds_proc, fates_bounds)

   call fates_hist%Init(nclumps, fates_bounds)

   ! Define the bounds on the first dimension for each thread
   !$OMP PARALLEL DO PRIVATE (nc,bounds_clump,fates_clump)
   do nc = 1,nclumps

      call get_clump_bounds(nc, bounds_clump)

      ! thread bounds for patch
      call hlm_bounds_to_fates_bounds(bounds_clump, fates_clump)
      call fates_hist%SetThreadBoundsEach(nc, fates_clump)
   end do
   !$OMP END PARALLEL DO

   ! ------------------------------------------------------------------------------------
   ! PART II: USE THE JUST DEFINED DIMENSIONS TO ASSEMBLE THE VALID IO TYPES
   ! INTERF-TODO: THESE CAN ALL BE EMBEDDED INTO A SUBROUTINE IN HISTORYIOMOD
   ! ------------------------------------------------------------------------------------
   call fates_hist%assemble_history_output_types()

   ! ------------------------------------------------------------------------------------
   ! PART III: DEFINE THE LIST OF OUTPUT VARIABLE OBJECTS, AND REGISTER THEM WITH THE
   ! HLM ACCORDING TO THEIR TYPES
   ! ------------------------------------------------------------------------------------
   call fates_hist%initialize_history_vars()
   nvar = fates_hist%num_history_vars()

   call CrossRefHistoryFields()
   
   do ivar = 1, nvar

      associate( vname    => fates_hist%hvars(ivar)%vname, &
                 vunits   => fates_hist%hvars(ivar)%units,   &
                 vlong    => fates_hist%hvars(ivar)%long, &
                 vdefault => fates_hist%hvars(ivar)%use_default, &
                 vavgflag => fates_hist%hvars(ivar)%avgflag)

        dk_index = fates_hist%hvars(ivar)%dim_kinds_index
        ioname = trim(fates_hist%dim_kinds(dk_index)%name)

        select case(trim(ioname))
        case(site_r8)
           call hist_addfld1d(fname=trim(vname),units=trim(vunits),         &
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_col=fates_hist%hvars(ivar)%r81d,      &
                              default=trim(vdefault))

       case(site_soil_r8, site_size_pft_r8, site_size_r8, site_pft_r8, &
             site_age_r8, site_height_r8, site_fuel_r8, site_cwdsc_r8, &
             site_can_r8,site_cnlf_r8, site_cnlfpft_r8, site_scag_r8, &
             site_scagpft_r8, site_agepft_r8, site_elem_r8, site_elpft_r8, &
             site_elcwd_r8, site_elage_r8, site_coage_r8, site_coage_pft_r8, &
             site_agefuel_r8,site_cdsc_r8, site_cdpf_r8, site_cdam_r8, &
             site_landuse_r8, site_lulu_r8)

           d_index = fates_hist%dim_kinds(dk_index)%dim2_index
           dim2name = fates_hist%dim_bounds(d_index)%name
           call hist_addfld2d(fname=trim(vname),units=trim(vunits),         &
                              type2d=trim(dim2name),                        &
                              avgflag=trim(vavgflag),long_name=trim(vlong), &
                              ptr_col=fates_hist%hvars(ivar)%r82d,     &
                              default=trim(vdefault))

        case default
           write(iulog,*) 'A FATES iotype was created that was not registerred'
           write(iulog,*) 'in CLM.:',trim(ioname)
           call endrun(msg=errMsg(sourcefile, __LINE__))
        end select

      end associate
   end do
 end subroutine init_history_io

 ! ======================================================================================

 subroutine init_soil_depths(this, nc)

    ! Input Arguments
    class(hlm_fates_interface_type), intent(inout) :: this
    integer,intent(in)                             :: nc   ! Clump

    ! Locals
    integer :: s  ! site index
    integer :: c  ! column index
    integer :: j  ! Depth index
    integer :: nlevsoil

    do s = 1, this%fates(nc)%nsites

       c = this%f2hmap(nc)%fcolumn(s)
       nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

       this%fates(nc)%bc_in(s)%zi_sisl(0:nlevsoil)    = col_pp%zi(c,0:nlevsoil)
       this%fates(nc)%bc_in(s)%dz_sisl(1:nlevsoil)    = col_pp%dz(c,1:nlevsoil)
       this%fates(nc)%bc_in(s)%z_sisl(1:nlevsoil)     = col_pp%z(c,1:nlevsoil)
       this%fates(nc)%bc_in(s)%dz_decomp_sisl(1:nlevdecomp) = &
            dzsoi_decomp(1:nlevdecomp)

       if (use_vertsoilc) then
          do j=1,nlevsoil
             this%fates(nc)%bc_in(s)%decomp_id(j) = j
             ! Check to make sure that dz = dz_decomp_sisl when vertical soil dynamics
             ! are active
             if(abs(this%fates(nc)%bc_in(s)%dz_decomp_sisl(j)-this%fates(nc)%bc_in(s)%dz_sisl(j))>1.e-10_r8)then
                write(iulog,*) 'when vertical soil decomp dynamics are on'
                write(iulog,*) 'fates assumes that the decomposition depths equal the soil depths'
                write(iulog,*) 'layer: ',j
                write(iulog,*) 'dz_decomp_sisl(j): ',this%fates(nc)%bc_in(s)%dz_decomp_sisl(j)
                write(iulog,*) 'dz_sisl(j): ',this%fates(nc)%bc_in(s)%dz_sisl(j)
                call endrun(msg=errMsg(sourcefile, __LINE__))
             end if
          end do
       else
          do j=1,nlevsoil
             this%fates(nc)%bc_in(s)%decomp_id(j) = 1
          end do
       end if


    end do

    return
 end subroutine init_soil_depths

 ! ======================================================================================

 subroutine ComputeRootSoilFlux(this, bounds_clump, num_filterc, filterc, &
       soilstate_inst)

    class(hlm_fates_interface_type), intent(inout) :: this
    type(bounds_type),intent(in)                   :: bounds_clump
    integer,intent(in)                             :: num_filterc
    integer,intent(in)                             :: filterc(num_filterc)
    type(soilstate_type), intent(inout)            :: soilstate_inst

    ! locals
    integer :: s
    integer :: c
    integer :: l
    integer :: nc
    integer :: num_filter_fates
    integer :: num_filter_hydroc
    integer :: nlevsoil


    if( .not. use_fates_planthydro ) return

    nc = bounds_clump%clump_index

    ! Perform a check that the number of columns submitted to fates for
    ! root water sink is the same that was expected in the hydrology filter
    num_filter_hydroc = 0
    do s = 1,num_filterc
       l = col_pp%landunit(filterc(s))
       if (lun_pp%itype(l) == istsoil ) then
          num_filter_hydroc = num_filter_hydroc + 1
       end if
    end do

    num_filter_fates = 0
    do s = 1, this%fates(nc)%nsites
       c = this%f2hmap(nc)%fcolumn(s)
       if(col_pp%active(c)) then
          num_filter_fates = num_filter_fates+1
       end if
    end do

    if(num_filter_fates .ne. num_filter_hydroc )then
       write(iulog,*) 'The HLM list of natural veg columns during root water transfer'
       write(iulog,*) 'is not the same size as the fates site list?'
       write(iulog,*) 'num_filter_fates: ',num_filter_hydroc
       write(iulog,*) 'nsites (active): ',num_filter_fates
       call endrun(msg=errMsg(sourcefile, __LINE__))
    end if

    do s = 1, this%fates(nc)%nsites
       c = this%f2hmap(nc)%fcolumn(s)
       nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

       ! This is the water removed from the soil layers by roots (or added)
       col_wf%qflx_rootsoi(c,1:nlevsoil) = &
            this%fates(nc)%bc_out(s)%qflx_soil2root_sisl(1:nlevsoil)

       ! This is the total amount of water transferred to surface runoff
       ! (this is generated potentially from supersaturating soils
       ! (currently this is unnecessary)
       ! col_wf%qflx_drain_vr_col(c,1:nlevsoil) = &
       !           this%fates(nc)%bc_out(s)%qflx_ro_sisl(1:nlevsoil)

    end do

 end subroutine ComputeRootSoilFlux

 ! ======================================================================================

 subroutine wrap_hydraulics_drive(this, bounds_clump, &
                                 fn, filterp, &
                                 soilstate_inst, &
                                 solarabs_inst, &
                                 energyflux_inst)


   implicit none
   class(hlm_fates_interface_type), intent(inout) :: this
   type(bounds_type),intent(in)                   :: bounds_clump
   integer, intent(in)                            :: fn
   integer, intent(in)                            :: filterp(fn)
   type(soilstate_type)    , intent(inout)        :: soilstate_inst
   type(solarabs_type)     , intent(in)           :: solarabs_inst
   type(energyflux_type)   , intent(inout)        :: energyflux_inst


   ! locals
   integer :: s
   integer :: c
   integer :: j
   integer :: f    ! filter loop index
   integer :: ifp
   integer :: p
   integer :: nc
   real(r8) :: dtime
   integer  :: nlevsoil


   if ( .not.use_fates_planthydro ) return

   nc = bounds_clump%clump_index
   dtime = real(get_step_size(),r8)

   ! Prepare Input Boundary Conditions
   ! ------------------------------------------------------------------------------------

   do s = 1, this%fates(nc)%nsites
      c = this%f2hmap(nc)%fcolumn(s)
      nlevsoil = this%fates(nc)%bc_in(s)%nlevsoil

      this%fates(nc)%bc_in(s)%smpmin_si                 = &
            soilstate_inst%smpmin_col(c)
      this%fates(nc)%bc_in(s)%watsat_sisl(1:nlevsoil)    = &
            soilstate_inst%watsat_col(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%watres_sisl(1:nlevsoil)    = &
           soilstate_inst%watmin_col(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%sucsat_sisl(1:nlevsoil)     = &
            soilstate_inst%sucsat_col(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%bsw_sisl(1:nlevsoil)        = &
            soilstate_inst%bsw_col(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%h2o_liq_sisl(1:nlevsoil)    = &
            col_ws%h2osoi_liq(c,1:nlevsoil)
      this%fates(nc)%bc_in(s)%eff_porosity_sl(1:nlevsoil) = &
            soilstate_inst%eff_porosity_col(c,1:nlevsoil)

      do ifp = 1, this%fates(nc)%sites(s)%youngest_patch%patchno
         p = ifp+col_pp%pfti(c)
         this%fates(nc)%bc_in(s)%swrad_net_pa(ifp) = solarabs_inst%fsa_patch(p)
         this%fates(nc)%bc_in(s)%lwrad_net_pa(ifp) = energyflux_inst%eflx_lwrad_net_patch(p)
      end do
   end do

   ! The exposed vegetation filter "filterp" dictates which patches
   ! had their transpiration updated during canopy_fluxes(). Patches
   ! not in the filter had been zero'd during prep_canopyfluxes().

   do f = 1,fn
      p = filterp(f)
      c = veg_pp%column(p)
      s = this%f2hmap(nc)%hsites(c)
      ifp = p - col_pp%pfti(c)
      this%fates(nc)%bc_in(s)%qflx_transp_pa(ifp) = veg_wf%qflx_tran_veg(p)
   end do

   ! Call Fates Hydraulics
   ! ------------------------------------------------------------------------------------


   call hydraulics_drive(this%fates(nc)%nsites, &
            this%fates(nc)%sites,  &
            this%fates(nc)%bc_in,  &
            this%fates(nc)%bc_out, &
            dtime)

   ! Prepare Output Boundary Conditions
   ! ------------------------------------------------------------------------------------

   do s = 1, this%fates(nc)%nsites
      c = this%f2hmap(nc)%fcolumn(s)
      col_ws%total_plant_stored_h2o(c) = &
            this%fates(nc)%bc_out(s)%plant_stored_h2o_si

   end do




   ! Update History Buffers that need to be updated after hydraulics calls

   call fates_hist%update_history_hydraulics(nc, &
         this%fates(nc)%nsites, &
         this%fates(nc)%sites, &
         this%fates(nc)%bc_in, &
         dtime)


   return
 end subroutine wrap_hydraulics_drive

 ! ======================================================================================

 subroutine hlm_bounds_to_fates_bounds(hlm, fates)

   use FatesIODimensionsMod, only : fates_bounds_type
   use FatesInterfaceTypesMod, only : nlevsclass_fates => nlevsclass
   use FatesInterfaceTypesMod, only : nlevage_fates    => nlevage
   use FatesInterfaceTypesMod, only : nlevheight_fates => nlevheight
   use FatesInterfaceTypesMod, only : nlevdamage_fates => nlevdamage
   use FatesFuelClassesMod,        only : nfc_fates   => num_fuel_classes
   use FatesLitterMod,    only : ncwd_fates       => ncwd
   use EDParamsMod,       only : nlevleaf_fates   => nlevleaf
   use EDParamsMod,       only : nclmax_fates     => nclmax
   use FatesInterfaceTypesMod, only : numpft_fates     => numpft
   use FatesInterfaceTypesMod, only : nlevcoage
   use FatesConstantsMod, only : n_landuse_cats

   implicit none

   type(bounds_type), intent(in)        :: hlm
   type(fates_bounds_type), intent(out) :: fates

   fates%cohort_begin = hlm%begcohort
   fates%cohort_end = hlm%endcohort

   fates%column_begin = hlm%begc
   fates%column_end = hlm%endc

   fates%soil_begin = 1
   fates%soil_end = nlevsoi

   fates%sizepft_class_begin = 1
   fates%sizepft_class_end = nlevsclass_fates * numpft_fates

   fates%size_class_begin = 1
   fates%size_class_end = nlevsclass_fates

   fates%pft_class_begin = 1
   fates%pft_class_end = numpft_fates

   fates%age_class_begin = 1
   fates%age_class_end = nlevage_fates

   fates%sizeage_class_begin = 1
   fates%sizeage_class_end   = nlevsclass_fates * nlevage_fates

   fates%fuel_begin = 1
   fates%fuel_end = nfc_fates

   fates%cdpf_begin = 1
   fates%cdpf_end = nlevdamage_fates * numpft_fates * nlevsclass_fates

   fates%cdsc_begin = 1
   fates%cdsc_end = nlevdamage_fates * nlevsclass_fates
   
   fates%cdam_begin = 1
   fates%cdam_end = nlevdamage_fates
   
   fates%cwdsc_begin = 1
   fates%cwdsc_end = ncwd_fates

   fates%can_begin = 1
   fates%can_end = nclmax_fates

   fates%cnlf_begin = 1
   fates%cnlf_end = nlevleaf_fates * nclmax_fates

   fates%cnlfpft_begin = 1
   fates%cnlfpft_end = nlevleaf_fates * nclmax_fates * numpft_fates

   fates%height_begin = 1
   fates%height_end = nlevheight_fates

   fates%agepft_class_begin = 1
   fates%agepft_class_end   = nlevage_fates * numpft_fates

   fates%sizeagepft_class_begin = 1
   fates%sizeagepft_class_end   = nlevsclass_fates * nlevage_fates * numpft_fates

   fates%elem_begin = 1
   fates%elem_end   = num_elements

   fates%elpft_begin = 1
   fates%elpft_end   = num_elements * numpft_fates

   fates%elcwd_begin = 1
   fates%elcwd_end   = num_elements * ncwd_fates

   fates%elage_begin = 1
   fates%elage_end   = num_elements * nlevage_fates

   fates%coagepf_class_begin = 1
   fates%coagepf_class_end = nlevcoage * numpft_fates

   fates%coage_class_begin = 1
   fates%coage_class_end = nlevcoage

   fates%agefuel_begin = 1
   fates%agefuel_end   = nlevage_fates * nfc_fates

   fates%landuse_begin = 1
   fates%landuse_end   = n_landuse_cats

   fates%lulu_begin = 1
   fates%lulu_end   = n_landuse_cats * n_landuse_cats

 end subroutine hlm_bounds_to_fates_bounds

! ======================================================================================

 subroutine GetAndSetTime()

   ! CLM MODULES
   use elm_time_manager  , only : get_days_per_year, &
                                  get_curr_date,     &
                                  get_ref_date,      &
                                  timemgr_datediff

   ! FATES MODULES
   use FatesInterfaceMod     , only : SetFatesTime

   ! LOCAL VARIABLES
   integer  :: yr                       ! year (0, ...)
   integer  :: mon                      ! month (1, ..., 12)
   integer  :: day                      ! day of month (1, ..., 31)
   integer  :: sec                      ! seconds of the day
   integer  :: current_year
   integer  :: current_month
   integer  :: current_day
   integer  :: current_tod
   integer  :: current_date
   integer  :: jan01_curr_year
   integer  :: reference_date
   integer  :: days_per_year
   real(r8) :: model_day
   real(r8) :: day_of_year


   ! Get the current date and determine the set the start of the current year
   call get_curr_date(current_year,current_month,current_day,current_tod)
   current_date = current_year*10000 + current_month*100 + current_day
   jan01_curr_year = current_year*10000 + 100 + 1

   ! Get the reference date components and compute the date
   call get_ref_date(yr, mon, day, sec)
   reference_date = yr*10000 + mon*100 + day

   ! Get the defined number of days per year
   days_per_year = get_days_per_year()

   ! Determine the model day
   call timemgr_datediff(reference_date, sec, current_date, current_tod, model_day)

   ! Determine the current DOY
   call timemgr_datediff(jan01_curr_year,0,current_date,sec,day_of_year)

   ! Set the FATES global time variables
   call SetFatesTime(current_year, current_month, &
                     current_day, current_tod, &
                     current_date, reference_date, &
                     model_day, floor(day_of_year), &
                     days_per_year, 1.0_r8/dble(days_per_year))

 end subroutine GetAndSetTime

! ======================================================================================

 subroutine GetLandusePFTData(bounds, landuse_pft_file, landuse_pft_map, landuse_bareground)

   ! !DESCRIPTION:
   ! Read in static landuse x pft file

   ! !USES:
   use fileutils , only : getfil
   use ncdio_pio , only : file_desc_t, ncd_io, ncd_inqdlen
   use ncdio_pio , only : ncd_pio_openfile, ncd_pio_closefile
   use decompMod , only : BOUNDS_LEVEL_PROC
   use elm_varcon, only : grlnd
   use FatesConstantsMod, only : fates_unset_r8


   ! !ARGUMENTS:
   type(bounds_type), intent(in)        :: bounds            ! proc-level bounds
   character(len=*) , intent(in)        :: landuse_pft_file  ! name of file containing static landuse x pft information
   real(r8), allocatable, intent(inout) :: landuse_pft_map(:,:,:)
   real(r8), allocatable, intent(inout) :: landuse_bareground(:)

   ! !LOCAL VARIABLES
   integer            :: varnum                    ! variable number
   integer            :: dimid, dimlen             ! dimension id number and length
   integer            :: ier                       ! error id
   character(len=256) :: locfn                     ! local file name
   type(file_desc_t)  :: ncid                      ! netcdf id
   real(r8), pointer  :: arraylocal(:,:)           ! local array for reading fraction data
   real(r8), pointer  :: arraylocal_bareground(:)  ! local array for reading bareground data
   logical            :: readvar                   ! true => variable is on dataset
   !character(len=16), parameter :: grlnd  = 'lndgrid'      ! name of lndgrid

   integer, parameter :: dim_landuse_pft = 14

   ! Land use name arrays
   character(len=10), parameter  :: landuse_pft_map_varnames(num_landuse_pft_vars) = &
                    [character(len=10)  :: 'frac_primr','frac_secnd','frac_pastr','frac_range'] !need to move 'frac_surf' to a different variable

   character(len=*), parameter :: subname = 'GetLandusePFTData'

   !-----------------------------------------------------------------------

   ! Check to see if the landuse file name has been provided
   ! Note: getfile checks this as well
   if (masterproc) then
      write(iulog,*) 'Attempting to read landuse x pft data .....'
      if (landuse_pft_file == ' ') then
         write(iulog,*)'landuse_pft_file must be specified'
         call endrun(msg=errMsg(__FILE__, __LINE__))
      endif
   endif

   ! Initialize the landuse x pft arrays and initialize to unset
    allocate(landuse_pft_map(bounds%begg:bounds%endg,dim_landuse_pft,num_landuse_pft_vars),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for landuse_pft_map'//errMsg(__FILE__, __LINE__))
    end if
    landuse_pft_map = fates_unset_r8

    allocate(landuse_bareground(bounds%begg:bounds%endg),stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for landuse_bareground'//errMsg(__FILE__, __LINE__))
    end if
    landuse_bareground = fates_unset_r8


   ! Get the local filename and open the file
   call getfil(landuse_pft_file, locfn, 0)
   call ncd_pio_openfile (ncid, trim(locfn), 0)

   ! Check that natpft dimension on the file matches the target array dimensions
   call ncd_inqdlen(ncid, dimid, dimlen, 'natpft')
   if (dimlen /= dim_landuse_pft) then
      write(iulog,*) 'natpft dimensions on the landuse x pft file do not match target array size'
      call endrun(msg=errMsg(__FILE__, __LINE__))
   end if

   ! Allocate a temporary array since ncdio expects a pointer
   allocate(arraylocal(bounds%begg:bounds%endg,dim_landuse_pft))
   allocate(arraylocal_bareground(bounds%begg:bounds%endg))

   ! Read the landuse x pft data from file
   do varnum = 1, num_landuse_pft_vars
      call ncd_io(ncid=ncid, varname=landuse_pft_map_varnames(varnum), flag='read', &
                  data=arraylocal, dim1name=grlnd, readvar=readvar)
      if (.not. readvar) &
         call endrun(msg='ERROR: '//trim(landuse_pft_map_varnames(varnum))// &
                         ' NOT on landuse x pft file'//errMsg(__FILE__, __LINE__))
      landuse_pft_map(bounds%begg:bounds%endg,:,varnum) = arraylocal(bounds%begg:bounds%endg,:)
   end do

   ! Read the bareground data from file.  This is per gridcell only.
   call ncd_io(ncid=ncid, varname='frac_brgnd', flag='read', &
               data=arraylocal_bareground, dim1name=grlnd, readvar=readvar)
   if (.not. readvar) call endrun(msg='ERROR: frac_brgnd NOT on landuse x pft file'//errMsg(__FILE__, __LINE__))
   landuse_bareground(bounds%begg:bounds%endg) = arraylocal_bareground(bounds%begg:bounds%endg)

   ! Deallocate the temporary local array point and close the file
   deallocate(arraylocal)
   deallocate(arraylocal_bareground)
   call ncd_pio_closefile(ncid)

 end subroutine GetLandusePFTData

end module ELMFatesInterfaceMod
