module BeTRSimulation
  !
  ! !DESCRIPTION:
  !  BeTR simulation base class.
  !
  !  BeTR simulation class are API definitions, mapping data
  !  structures from a specific LSM, e.g. CLM, ALM, into BeTR data
  !  structures.
  !

  use abortutils     , only : endrun
  use elm_varctl     , only : iulog, use_cn
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use tracer_varcon  , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use BeTR_decompMod , only : betr_bounds_type
  use decompMod      , only : bounds_type
#if (defined SBETR)
  use PatchType      , only : patch_type
  use ColumnType     , only : column_type
  use LandunitType   , only : landunit_type
  use CNCarbonFluxType  , only : carbonflux_type
  use PfCarbonFluxType  , only : pf_carbonflux_type
  use WaterStateType , only : waterstate_type
  use WaterfluxType     , only : waterflux_type
  use TemperatureType   , only : temperature_type
  use PfTemperatureType   , only : pf_temperature_type
  use PfWaterfluxType     , only : pf_waterflux_type
#else
  use ColumnType     , only : column_type => column_physical_properties
  use VegetationType , only : patch_type  => vegetation_physical_properties
  use LandunitType   , only : landunit_type => landunit_physical_properties
  use ColumnDataType      , only : carbonflux_type => column_carbon_flux
  use ColumnDataType , only : waterstate_type => column_water_state
  use ColumnDataType , only : waterflux_type => column_water_flux
  use VegetationDataType, only : veg_es, veg_wf
  use ColumnDataType    , only : col_es, col_ws, col_wf
  use ColumnDataType, only : temperature_type=> column_energy_state
  use VegetationDataType  , only : pf_carbonflux_type => vegetation_carbon_flux
  use VegetationDataType, only : pf_temperature_type => vegetation_energy_state
  use VegetationDataType, only : pf_waterflux_type => vegetation_water_flux
#endif

  ! !USES:
  use BetrType                 , only : betr_type, create_betr_type
  use betr_ctrl                , only : max_betr_hist_type, betr_offline, max_betr_rest_type, biulog
  use betr_constants           , only : betr_string_length
  use betr_constants           , only : betr_filename_length
  use betr_regression_module   , only : betr_regression_type
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type, create_betr_biogeophys_input
  use BeTR_biogeoStateType     , only : betr_biogeo_state_type, create_betr_biogeo_state
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type, create_betr_biogeoFlux
  use BeTR_TimeMod             , only : betr_time_type
  use betr_varcon              , only : betr_maxpatch_pft
  use BetrStatusType           , only : betr_status_type, create_betr_status_type
  use BetrStatusSimType        , only : betr_status_sim_type, create_betr_status_sim_type
  use betr_columnType          , only : betr_column_type, create_betr_column_type
  use betr_patchType           , only : betr_patch_type, create_betr_patch_type
  use betr_varcon              , only : spval => bspval, ispval=>bispval
  use BeTRHistVarType          , only : betr_hist_var_type
  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: betr_simulation_type
     type(betr_type)                     , public, pointer  :: betr(:)  => null()
     type(betr_biogeophys_input_type)    , public, pointer  :: biophys_forc(:) => null()
     type(betr_biogeo_state_type)        , public, pointer  :: biogeo_state(:) => null()
     type(betr_biogeo_flux_type)         , public, pointer  :: biogeo_flux(:) => null()
     type(betr_status_type)              , public, pointer  :: bstatus(:) => null()
     type(betr_status_sim_type)          , public, pointer  :: bsimstatus
     logical                             , public, pointer :: active_col(:) => null()
     type(betr_column_type)              , public, pointer :: betr_col(:) => null()
     type(betr_patch_type)               , public, pointer :: betr_pft(:) => null()
     type(betr_bounds_type)              , public, pointer :: betr_bounds(:) => null()
     character(len=betr_filename_length) , private :: base_filename
     character(len=betr_filename_length) , private :: hist_filename
     character(len=64)                   , private :: case_id

     type(betr_regression_type), private          :: regression
     type(betr_time_type), public                 :: betr_time
     integer, public                              :: num_soilp
     integer, public, allocatable                 :: filter_soilp(:)
     integer, public                              :: num_jtops
     integer, public, allocatable                 :: jtops(:)
     integer, public                              :: num_surfc
     integer, public                              :: num_parcols
     integer, public, allocatable                 :: filter_soilc(:)
     integer, public :: spinup_count
     type(betr_hist_var_type), allocatable :: state_hist1d_var(:)
     type(betr_hist_var_type), allocatable :: state_hist2d_var(:)
     type(betr_hist_var_type), allocatable :: flux_hist1d_var(:)
     type(betr_hist_var_type), allocatable :: flux_hist2d_var(:)
     integer :: num_hist_state1d
     integer :: num_hist_state2d
     integer :: num_hist_flux1d
     integer :: num_hist_flux2d

     integer :: num_rest_state1d
     integer :: num_rest_state2d
     real(r8), pointer :: hist_states_2d(:,:,:)
     real(r8), pointer :: hist_states_1d(:,:)
     real(r8), pointer :: rest_states_2d(:,:,:)
     real(r8), pointer :: rest_states_1d(:,:)
     real(r8), pointer :: hist_fluxes_1d(:,:)
     real(r8), pointer :: hist_fluxes_2d(:,:,:)
     real(r8), pointer :: hist_fluxes_1d_accum(:,:)
     real(r8), pointer :: hist_fluxes_2d_accum(:,:,:)
     real(r8), pointer :: scalaravg_col(:)
     real(r8), pointer :: dom_scalar_col(:)
     logical,  private :: active_soibgc
     real(r8), private :: hist_naccum
     integer , private :: hist_record

     ! FIXME(bja, 201603) most of these types should be private!

     ! NOTE(bja, 201603) BeTR types only, no LSM specific types here!

   contains
     procedure, public :: BeTRInit
     procedure, public :: BeTRSetFilter
     procedure, public :: SetClock
     procedure, public :: InitOnline              => BeTRSimulationInit
     procedure, public :: Init                    => BeTRSimulationInitOffline
     procedure, public :: ConsistencyCheck        => BeTRSimulationConsistencyCheck
     procedure, public :: PreDiagSoilColWaterFlux => BeTRSimulationPreDiagSoilColWaterFlux
     procedure, public :: DiagnoseDtracerFreezeThaw=> BeTRSimulationDiagnoseDtracerFreezeThaw
     procedure, public :: DiagAdvWaterFlux        => BeTRSimulationDiagAdvWaterFlux
     procedure, public :: DiagDrainWaterFlux      => BeTRSimulationDiagDrainWaterFlux
     procedure, public :: BeginSnowLayerAdjst     => BeTRSimulationBeginTracerSnowLayerAdjst
     procedure, public :: EndSnowLayerAdjst       => BeTRSimulationEndTracerSnowLayerAdjst
     procedure, public :: CombineSnowLayers       => BeTRSimulationCombineSnowLayers
     procedure, public :: DvideSnowLayers         => BeTRSimulationDvideSnowLayers
     procedure, public :: StepWithoutDrainage     => BeTRSimulationStepWithoutDrainage
     procedure, public :: StepWithDrainage        => BeTRSimulationStepWithDrainage
     procedure, public :: BeginMassBalanceCheck   => BeTRSimulationBeginMassBalanceCheck
     procedure, public :: MassBalanceCheck        => BeTRSimulationMassBalanceCheck
     procedure, public :: BeTRSetBiophysForcing   => BeTRSimulationSetBiophysForcing
     procedure, public :: RetrieveBiogeoFlux      => BeTRSimulationRetrieveBiogeoFlux
     procedure, public :: DiagnoseLnd2atm         => BeTRSimulationDiagnoseLnd2atm
     procedure, public :: CreateOfflineHistory    => hist_htapes_create
     procedure, public :: WriteOfflineHistory     => hist_write
     procedure, public :: SetSpinup               => BeTRSimulationSetSpinup
     procedure, public :: ReadParams              => BeTRSimulationReadParams
     procedure, public :: WriteRegressionOutput
     !the following are used to interact with lsm
     procedure, public :: do_soibgc
     procedure, public :: BeTRRestart             => BeTRSimulationRestart
     procedure, public :: BeTRRestartOffline      => BeTRSimulationRestartOffline
     procedure, public :: BeTRRestartOpen         => BeTRSimulationRestartOpen
     procedure, public :: BeTRRestartClose        => BeTRSimulationRestartClose
     procedure, private :: BeTRCreateHistory      => BeTRSimulationCreateHistory
     procedure, public  :: HistRetrieval          => BeTRSimulationHistRetrieval
     procedure, private :: BeTRRetrieveHistoryState    => BeTRSimulationRetrieveHistoryState
     procedure, private :: BeTRRetrieveHistoryFlux    => BeTRSimulationRetrieveHistoryFlux
     procedure, public :: BeTRSetcps              => BeTRSimulationSetcps
     procedure, public :: BeTRSetBounds           => BeTRSimulationSetBounds
     procedure, private:: hist_create_states
     procedure, private:: hist_create_fluxes
     procedure, private:: hist_output_states
     procedure, private:: hist_output_fluxes
     procedure, private:: RestAlloc               => BeTRSimulationRestartAlloc
     procedure, private:: HistAlloc               => BeTRSimulationHistoryAlloc
     procedure, private:: set_activecol
     procedure, public :: do_regress_test
     procedure, private :: hist_flux_accum
  end type betr_simulation_type

contains

  subroutine SetClock(this, dtime, nelapstep)
  implicit none
  class(betr_simulation_type)              , intent(inout) :: this
  real(r8), intent(in) :: dtime
  integer, intent(in) :: nelapstep


  call this%betr_time%setClock(dtime, nelapstep)

  end subroutine SetClock
  !-------------------------------------------------------------------------------
  subroutine BeTRSimulationInit(this, bounds, lun, col, pft, waterstate, namelist_buffer, masterproc)
    !
    ! DESCRIPTIONS
    ! Dummy routine for inheritance purposes. don't use.
    !
    !USES
    use betr_constants , only : betr_namelist_buffer_size, betr_filename_length
    implicit none

    class(betr_simulation_type)              , intent(inout) :: this
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(inout) :: col
    type(patch_type)                         , intent(in) :: pft
    type(bounds_type)                        , intent(in)    :: bounds
    character(len=*)                         , intent(in)    :: namelist_buffer
    type(waterstate_type)                    , intent(inout) :: waterstate
    logical,                        optional , intent(in)    :: masterproc
    character(len=*), parameter :: subname = 'BeTRSimulationInit'

    call endrun(msg="ERROR "//subname//" unimplemented. "//errmsg(mod_filename, __LINE__))

    if (this%num_surfc > 0)                  continue
    if (bounds%begc > 0)                     continue
    if (size(waterstate%h2osoi_liq) > 0) continue
    if (len(namelist_buffer) > 0)            continue
  end subroutine BeTRSimulationInit

  !-------------------------------------------------------------------------------
  subroutine BeTRSimulationInitOffline(this, bounds, lun, col, pft, waterstate, &
       namelist_buffer, base_filename, case_id)
    !
    ! DESCRIPTIONS
    ! Dummy routine for inheritance purposes. don't use.
    !
    !USES
    use betr_constants , only : betr_namelist_buffer_size, betr_filename_length
    implicit none

    class(betr_simulation_type)              , intent(inout) :: this
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(inout) :: col
    type(patch_type)                         , intent(in) :: pft
    type(bounds_type)                        , intent(in)    :: bounds
    type(waterstate_type)                    , intent(inout) :: waterstate
    character(len=*)                         , intent(in)    :: namelist_buffer
    character(len=*)                         , intent(in)    :: base_filename
    character(len=*)                         , intent(in)    :: case_id

    character(len=*), parameter :: subname = 'BeTRSimulationInit'

    call endrun(msg="ERROR "//subname//" unimplemented. "//errmsg(mod_filename, __LINE__))

    if (this%num_surfc > 0)                  continue
    if (bounds%begc > 0)                     continue
    if (size(waterstate%h2osoi_liq) > 0) continue
    if (len(base_filename) > 0)              continue
    if (len(namelist_buffer) > 0)            continue
    call this%betr_time%Init(namelist_buffer)
  end subroutine BeTRSimulationInitOffline
!-------------------------------------------------------------------------------
  subroutine BeTRSetFilter(this, maxpft_per_col, nsoilorder, boffline)
  !
  !DESCRIPTION
  ! set betr filter, only used for standalone applicaitons
    use betr_ctrl                , only : betr_offline
  use betr_varcon                , only : betr_max_soilorder
  implicit none
  !ARGUMENTS
  class(betr_simulation_type), intent(inout) :: this
  integer, intent(in) :: maxpft_per_col
  integer, optional, intent(in) :: nsoilorder
  logical, optional, intent(in) :: boffline
  integer :: p

  !by default, surface litter layer is off
  this%num_jtops = 1
  allocate(this%jtops(this%num_jtops))
  this%jtops(:) = 1

  this%num_parcols = 1
  this%num_surfc = 1
  allocate(this%filter_soilc(this%num_surfc))
  this%filter_soilc(:) = 1

  this%num_soilp = maxpft_per_col
  if(this%num_soilp>0)allocate(this%filter_soilp(this%num_soilp))
  do p = 1, maxpft_per_col
    this%filter_soilp(p) = p
  enddo
  if(present(boffline))then
    betr_offline=boffline
  else
    betr_offline=.true.
  endif
  betr_maxpatch_pft = maxpft_per_col
  if(present(nsoilorder))then
    betr_max_soilorder= nsoilorder
  else
    betr_max_soilorder=1
  endif
  end subroutine BeTRSetFilter

!-------------------------------------------------------------------------------
  subroutine BeTRSimulationReadParams(this, bounds)

  use ncdio_pio               , only :  file_desc_t
  use ncdio_pio               , only : ncd_pio_closefile, ncd_pio_openfile, &
                                         file_desc_t, ncd_inqdid, ncd_inqdlen
  use ApplicationsFactory      , only : AppLoadParameters, AppCopyParas
  use BetrStatusType           , only : betr_status_type
  use decompMod                , only : bounds_type
  use tracer_varcon            , only : bgc_param_file
  use fileutils                , only : getfil
  use spmdMod                  , only : masterproc
  implicit none
  class(betr_simulation_type), intent(inout) :: this

  type(bounds_type), intent(in) :: bounds
  !temporary variables
  type(betr_status_type)   :: bstatus
  type(betr_bounds_type)   :: betr_bounds
  integer :: c
  character(len=256) :: locfn ! local file name
  type(file_desc_t)  :: ncid  ! pio netCDF file id

  !open file for parameter reading
  call getfil (bgc_param_file, locfn, 0)
  if (masterproc) then
    write(iulog,*) 'read betr bgc parameter file '//trim(locfn)
  endif
  call ncd_pio_openfile (ncid, trim(locfn), 0)

  !read in parameters
  call AppLoadParameters(ncid, bstatus)

  call ncd_pio_closefile(ncid)

  if(bstatus%check_status())then
    call endrun(msg=bstatus%print_msg())
  endif

  !grid horizontal bounds
  call this%BeTRSetBounds(betr_bounds)

  if(this%num_parcols==1)then
    call AppCopyParas(1, 1, bstatus)
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call this%betr(c)%SetParCols(1)
    enddo
  else
    call AppCopyParas(bounds%begc, bounds%endc, bstatus)
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call this%betr(c)%SetParCols(c)
    enddo
  endif
  if(bstatus%check_status())then
    call endrun(msg=bstatus%print_msg())
  endif

  do c = bounds%begc, bounds%endc
    if(.not. this%active_col(c))cycle
    call this%betr(c)%UpdateParas(betr_bounds, this%bstatus(c))
    if(this%bstatus(c)%check_status())then
      call this%bsimstatus%setcol(c)
      call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
      exit
    endif
  enddo
  if(this%bsimstatus%check_status())call endrun(msg=this%bsimstatus%print_msg())

  end subroutine BeTRSimulationReadParams
!-------------------------------------------------------------------------------
  subroutine set_activecol(this, col)

  implicit none
    class(betr_simulation_type)              , intent(inout) :: this
    type(column_type)                        , intent(inout) :: col

    logical :: do_debug
    integer :: cc
    do_debug=.false.
    if(do_debug)then
      cc=1
      col%active(:) =.false.
      col%active(cc)=.true.
      this%betr(cc)%tracers%debug=.true.
    endif
  !    col%active(c_act)=.true.

  end subroutine set_activecol
!-------------------------------------------------------------------------------
  subroutine BeTRInit(this, bounds, lun, col, pft, waterstate, namelist_buffer, &
     base_filename, case_id, masterproc)
    !
    ! DESCRIPTION
    ! initialize BeTR
    !
    !!USES
    use betr_constants , only : betr_namelist_buffer_size
    use betr_constants , only : betr_filename_length
    use betr_varcon    , only : betr_maxpatch_pft
    use landunit_varcon, only : istsoil, istcrop
    implicit none
    !ARGUMENTS
    class(betr_simulation_type)              , intent(inout) :: this
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(inout) :: col
    type(patch_type)                         , intent(in) :: pft
    type(waterstate_type)                    , intent(in) :: waterstate
    character(len=*)                         , intent(in) :: namelist_buffer
    character(len=*)               , optional, intent(in) :: base_filename
    character(len=*)               , optional, intent(in) :: case_id
    logical,                      optional   , intent(in) :: masterproc
    !TEMPORARY VARIABLES
    character(len=*), parameter :: subname = 'BeTRInit'
    type(betr_bounds_type) :: betr_bounds
    integer :: c, l, c_l
    logical :: asoibgc
    !print*,'base_filename',trim(base_filename)

    this%hist_record=0
    this%active_soibgc=.false.
    biulog = iulog
    if(present(base_filename))then
      this%base_filename = base_filename
    else
      this%base_filename = 'sbetr'
    endif

    if(present(case_id))then
      this%case_id=case_id
    else
      this%case_id='exp0'
    endif
    if(present(masterproc))then
      call this%betr_time%Init(namelist_buffer, masterproc)
    else
      call this%betr_time%Init(namelist_buffer)
    endif
    !allocate memory
    allocate(this%betr(bounds%begc:bounds%endc))
    allocate(this%biophys_forc(bounds%begc:bounds%endc))
    allocate(this%biogeo_flux(bounds%begc:bounds%endc))
    allocate(this%biogeo_state(bounds%begc:bounds%endc))
    allocate(this%bstatus(bounds%begc:bounds%endc))
    allocate(this%betr_col(bounds%begc:bounds%endc))
    allocate(this%betr_pft(bounds%begc:bounds%endc))
    allocate(this%active_col(bounds%begc:bounds%endc))
    allocate(this%bsimstatus)

    call this%bsimstatus%reset()

    !grid horizontal bounds
    call this%BeTRSetBounds(betr_bounds)

    do c = bounds%begc, bounds%endc
      l = col%landunit(c)
      call this%biophys_forc(c)%Init(betr_bounds)

      call this%betr_col(c)%Init(betr_bounds)

      call this%betr_pft(c)%Init(betr_bounds)

      if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
        this%active_col(c) = .true.
      else
        this%active_col(c) = .false.
      endif
    enddo
    call this%BeTRSetcps(bounds, col, pft)
    !the following pass in of waterstate_vars is needed for doing h2o isotopes
    call this%BeTRSetBiophysForcing(bounds, col, pft, betr_bounds%lbj, betr_bounds%ubj, &
      waterstate_vars = waterstate)

    do c = bounds%begc, bounds%endc
      call this%betr(c)%Init(namelist_buffer, betr_bounds, this%betr_col(c), &
          this%biophys_forc(c), asoibgc, this%bstatus(c))
      if(c==bounds%begc)this%active_soibgc=asoibgc
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo

    call this%set_activecol(col)

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call this%betr(c)%UpdateParas(betr_bounds, this%bstatus(c))
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo

    !allocate spinup factor
    if(this%do_soibgc())then
      allocate(this%scalaravg_col(bounds%begc:bounds%endc));this%scalaravg_col(:) = 0._r8
      allocate(this%dom_scalar_col(bounds%begc:bounds%endc));this%dom_scalar_col(:)=1._r8
      c_l=1
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle
        call this%betr(c)%set_bgc_spinup(betr_bounds, 1,  betr_bounds%ubj, this%biophys_forc(c))
        this%dom_scalar_col(c)=this%biophys_forc(c)%dom_scalar_col(c_l)
      enddo
    endif

    if(this%bsimstatus%check_status())call endrun(msg=this%bsimstatus%print_msg())

    do c = bounds%begc, bounds%endc

      call this%biogeo_state(c)%Init(betr_bounds, this%active_soibgc)

      call this%biogeo_flux(c)%Init(betr_bounds, this%active_soibgc)
    enddo
    !identify variables that are used for history output
    c = bounds%begc
    call this%betr(c)%get_hist_size(this%num_hist_state1d, this%num_hist_state2d, &
      this%num_hist_flux1d, this%num_hist_flux2d)

    !allocate memory for history variables
    call this%HistAlloc(betr_bounds%lbj, betr_bounds%ubj, bounds)

    call this%betr(c)%get_hist_info(this%num_hist_state1d, this%num_hist_state2d, &
      this%num_hist_flux1d, this%num_hist_flux2d, &
      this%state_hist1d_var(1:this%num_hist_state1d), this%state_hist2d_var(1:this%num_hist_state2d), &
      this%flux_hist1d_var(1:this%num_hist_flux1d), this%flux_hist2d_var(1:this%num_hist_flux2d))

    if(this%bstatus(c)%check_status())then
      call this%bsimstatus%setcol(c)
      call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
      if(this%bsimstatus%check_status())call endrun(msg=this%bsimstatus%print_msg())
    endif

    if(betr_offline)then
      call this%CreateOfflineHistory(bounds, betr_bounds%ubj, &
         this%num_hist_state1d, this%num_hist_state2d, &
            this%num_hist_flux1d, this%num_hist_flux2d)
    else
      call this%BeTRCreateHistory(bounds, betr_bounds%ubj, &
         this%num_hist_state1d, this%num_hist_state2d, &
            this%num_hist_flux1d, this%num_hist_flux2d)
    endif
    !identify restart variables
    call this%betr(c)%get_restartvar_size(this%num_rest_state1d, this%num_rest_state2d)

    !allocate memory of passing restart variables
    call this%RestAlloc(betr_bounds%lbj, betr_bounds%ubj, bounds)

    if(present(base_filename)) then
      call this%regression%Init(base_filename, namelist_buffer, this%bsimstatus)
      if(this%bsimstatus%check_status())call endrun(msg=this%bsimstatus%print_msg())
    endif

  end subroutine BeTRInit
  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationRestartAlloc(this, lbj, ubj, bounds)
  implicit none
  !ARGUMENTS
  class(betr_simulation_type)         , intent(inout) :: this
  integer                             , intent(in)    :: lbj, ubj
  type(bounds_type)                   , intent(in)    :: bounds

  integer :: begc, endc

  begc = bounds%begc; endc=bounds%endc

  allocate(this%rest_states_1d(begc:endc, 1:this%num_rest_state1d))
  this%rest_states_1d(:,:)=spval
  allocate(this%rest_states_2d(begc:endc, 1:ubj, 1:this%num_rest_state2d))
  this%rest_states_2d(:,:,:)=spval

  end subroutine BeTRSimulationRestartAlloc
  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationHistoryAlloc(this,lbj, ubj, bounds)
  implicit none
  !ARGUMENTS
  class(betr_simulation_type)         , intent(inout) :: this
  integer                             , intent(in)    :: lbj, ubj
  type(bounds_type)                   , intent(in)    :: bounds

  integer :: begc, endc

  begc = bounds%begc; endc=bounds%endc
  !state variables
  if(this%num_hist_state1d>0)allocate(this%state_hist1d_var(this%num_hist_state1d))
  if(this%num_hist_state2d>0)allocate(this%state_hist2d_var(this%num_hist_state2d))

  if(this%num_hist_state2d>0)allocate(this%hist_states_2d(begc:endc, 1:ubj, 1:this%num_hist_state2d))
  if(this%num_hist_state1d>0)allocate(this%hist_states_1d(begc:endc, 1:this%num_hist_state1d))

  !flux variables
  if(this%num_hist_flux1d>0)allocate(this%flux_hist1d_var(this%num_hist_flux1d))
  if(this%num_hist_flux2d>0)allocate(this%flux_hist2d_var(this%num_hist_flux2d))

  if(this%num_hist_flux2d>0)allocate(this%hist_fluxes_2d(begc:endc, 1:ubj, 1:this%num_hist_flux2d))
  if(this%num_hist_flux1d>0)allocate(this%hist_fluxes_1d(begc:endc, 1:this%num_hist_flux1d))

  if(betr_offline)then
    if(this%num_hist_flux2d>0)allocate(this%hist_fluxes_2d_accum(begc:endc, 1:ubj, 1:this%num_hist_flux2d))
    if(this%num_hist_flux1d>0)allocate(this%hist_fluxes_1d_accum(begc:endc, 1:this%num_hist_flux1d))
    this%hist_fluxes_1d_accum(:,:) = 0._r8
    this%hist_fluxes_2d_accum(:,:,:) = 0._r8
    this%hist_naccum = 0._r8
  endif
  end subroutine BeTRSimulationHistoryAlloc


  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationRestartOpen(this, fname, flag, ncid)
    !
    !! DESCRIPTION
    ! open restart file, note it is only used in the stand alone mode
    ! note
    ! !USES:
    use bncdio_pio, only : file_desc_t, ncd_nowrite, ncd_pio_openfile, ncd_pio_createfile
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    character(len=*), intent(in) :: fname
    type(file_desc_t)          , intent(out) :: ncid ! netcdf id
    character(len=*), intent(in) :: flag

    integer :: c
    print*,'open restart file ',trim(fname), ' for ',trim(flag)
    if(trim(flag)=='read')then
      call ncd_pio_openfile(ncid, trim(fname), ncd_nowrite)
    elseif(trim(flag)=='write')then
      call ncd_pio_createfile(ncid, trim(fname))
    endif
    print*,'creating file succeeded'

  end subroutine BeTRSimulationRestartOpen


  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationRestartClose(this, ncid)
    !
    !! DESCRIPTION
    ! initialize for restart run
    ! !USES:
    use bncdio_pio, only : file_desc_t, ncd_pio_closefile
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(file_desc_t)          , intent(inout) :: ncid ! netcdf id

    call ncd_pio_closefile(ncid)
  end subroutine BeTRSimulationRestartClose
  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithoutDrainage(this, bounds, col, pft)
  !DESCRPTION
  !interface for StepWithoutDrainage
  !
  ! USES
    use SoilStateType     , only : soilstate_type
    use ChemStateType     , only : chemstate_type
    use atm2lndType       , only : atm2lnd_type
    use SoilHydrologyType , only : soilhydrology_type
    use CanopyStateType   , only : canopystate_type
    use BeTR_TimeMod      , only : betr_time_type
    use pftvarcon         , only : crop
    implicit none
  !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds ! bounds
    type(column_type)           , intent(in)    :: col ! column type
    type(patch_type)            , intent(in)    :: pft

    ! remove compiler warnings about unused dummy args
    if (this%num_surfc > 0)                           continue
    if (this%betr_time%tstep > 0)                          continue
    if (bounds%begc > 0)                              continue
  !  if (size(col%z) > 0)                              continue

  end subroutine BeTRSimulationStepWithoutDrainage
  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationDiagnoseLnd2atm(this, bounds,  col, lnd2atm_vars)
   !
   !DESCRIPTION
   !interface for using diagnose land fluxes to atm and river copmonents
   !
   !USES
    use MathfuncMod   , only : safe_div
    use lnd2atmType    , only : lnd2atm_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    type(column_type)           , intent(in)    :: col ! column type
    type(lnd2atm_type)          , intent(inout) :: lnd2atm_vars

    ! remove compiler warnings about unused dummy args
    if (this%num_surfc > 0) continue
    if (bounds%begc > 0)    continue
!    if (size(col%z) > 0)    continue

  end subroutine BeTRSimulationDiagnoseLnd2atm
  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationStepWithDrainage(this, bounds,  col)
   !
   !DESCRIPTION
   !interface for using StepWithDrainage
   !
   !USES
    use MathfuncMod   , only : safe_div
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    type(column_type)           , intent(in)    :: col ! column type

    ! remove compiler warnings about unused dummy args
    if (this%num_surfc > 0) continue
    if (bounds%begc > 0)    continue
!    if (size(col%z) > 0)    continue

  end subroutine BeTRSimulationStepWithDrainage

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationBeginMassBalanceCheck(this, bounds)
    !DESCRIPTION
    !stage tracer mass balance check
    !
    !USES
    use TracerBalanceMod, only : begin_betr_tracer_massbalance
    implicit none
    !ARGUMENTS
    class(betr_simulation_type), intent(inout)   :: this
    type(bounds_type), intent(in) :: bounds

    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj
    integer :: c

    !set lbj and ubj
    call this%BeTRSetBounds(betr_bounds)

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      call begin_betr_tracer_massbalance(betr_bounds,             &
         this%betr_col(c), this%num_surfc, this%filter_soilc,     &
         this%num_soilp, this%filter_soilp, &
         this%betr(c)%tracers, this%betr(c)%tracerstates,         &
         this%betr(c)%tracerfluxes, this%bstatus(c))
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=trim(this%bsimstatus%print_msg()))
  end  subroutine BeTRSimulationBeginMassBalanceCheck
  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationMassBalanceCheck(this, bounds)
   !DESCRIPTION
   ! do tracer mass balance check
   !
   !USES
    use TracerBalanceMod , only : betr_tracer_massbalance_check
    use BeTR_TimeMod     , only : betr_time_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    !TEMPORARY VARIABLES
    type(betr_bounds_type) :: betr_bounds
    integer :: c
    logical :: ldebug
    !set lbj and ubj
    call this%BeTRSetBounds(betr_bounds)

    do c = bounds%begc, bounds%endc
      if(this%betr(c)%skip_mass_bal_check() .and. this%betr_time%is_first_step())cycle
      if(.not. this%active_col(c))cycle
      call betr_tracer_massbalance_check(this%betr_time, betr_bounds,           &
         this%betr_col(c), this%num_surfc, this%filter_soilc,           &
         this%betr(c)%tracers, this%betr(c)%tracerstates,                &
         this%betr(c)%tracerfluxes, this%bstatus(c))!, ldebug)
      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
    enddo
    if(this%bsimstatus%check_status()) then
      print*,this%bsimstatus%cindex
      print*,trim(this%bsimstatus%print_msg())
      call endrun(msg=trim(this%bsimstatus%print_msg()))
    endif
  end subroutine BeTRSimulationMassBalanceCheck

!-------------------------------------------------------------------------------

  subroutine hist_htapes_create(this, bounds, betr_nlevtrc_soil, &
     num_hist_state1d, num_hist_state2d, num_hist_flux1d, num_hist_flux2d)
  !
  ! DESCRIPTIONS
  ! create history file and define output variables, only for standalone applicaitons
  !
  ! USES
    use netcdf          , only : nf90_float
    use bncdio_pio       , only : file_desc_t,ncd_nowrite
    use bncdio_pio       , only : ncd_pio_createfile
    use bncdio_pio       , only : ncd_pio_closefile
    use bncdio_pio       , only : ncd_enddef
    use bncdio_pio       , only : ncd_defvar
    use bncdio_pio       , only : ncd_putvar
    use bncdio_pio       , only : get_dim_len
    use bhistFileMod    , only : hist_file_create, hist_def_fld1d, hist_def_fld2d
    use betr_varcon     , only : bspval
    use betr_columnType , only : betr_column_type
    use betr_ctrl       , only : continue_run
    !
    !ARGUMENTS
    implicit none
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds               ! bounds
    integer           , intent(in) :: betr_nlevtrc_soil
    integer           ,     intent(in)   :: num_hist_state1d
    integer           ,     intent(in)   :: num_hist_state2d
    integer           ,     intent(in)   :: num_hist_flux1d
    integer           ,     intent(in)   :: num_hist_flux2d

  !TEMPORARY VARIABLES
    integer                     :: jj, kk, c
    integer :: ncol
    type(file_desc_t)           :: ncid
    character(len=*), parameter :: subname = 'hist_htapes_create'

    c = 1
    associate(                                                     &
         ntracers          => this%betr(c)%tracers%ntracers,          &
         ngwmobile_tracers => this%betr(c)%tracers%ngwmobile_tracers, &
         is_volatile       => this%betr(c)%tracers%is_volatile,       &
         is_h2o            => this%betr(c)%tracers%is_h2o,            &
         is_isotope        => this%betr(c)%tracers%is_isotope,        &
         volatileid        => this%betr(c)%tracers%volatileid,        &
         tracernames       => this%betr(c)%tracers%tracernames        &
         )

    ncol = bounds%endc-bounds%begc + 1
    if(len(trim(this%case_id))>0)then
      this%hist_filename = trim(this%base_filename) //'.'//trim(this%case_id)// '.output.nc'
    else
      this%hist_filename = trim(this%base_filename) //'.output.nc'
    endif
    if(continue_run)then
      this%hist_record = get_dim_len(this%hist_filename, 'time')
      return
    endif
    call ncd_pio_createfile(ncid, this%hist_filename)

    call hist_file_create(ncid, betr_nlevtrc_soil, ncol)

    call ncd_defvar(ncid, "ZSOI", nf90_float, &
        dim1name="ncol",dim2name="levtrc",    &
        long_name="grid center"           ,    &
        units="m", missing_value=bspval, fill_value=bspval)

    call  hist_def_fld2d(ncid, varname="QFLX_ADV", nf90_type=nf90_float,dim1name="ncol",     &
           dim2name="levtrc",long_name="advective flux / velocity", units="m/s")

    write(iulog,*)'hist_create_states'
    call this%hist_create_states(bounds, betr_nlevtrc_soil, num_hist_state1d, num_hist_state2d, ncid=ncid)

    write(iulog,*)'hist_create_fluxes'
    call this%hist_create_fluxes(bounds, betr_nlevtrc_soil, num_hist_flux1d, num_hist_flux2d, ncid=ncid)

    call ncd_enddef(ncid)
    call ncd_putvar(ncid,"ZSOI",1,this%betr_col(c)%z(1:1,1:betr_nlevtrc_soil))
    call ncd_pio_closefile(ncid)

    end associate
  end subroutine hist_htapes_create

  !-------------------------------------------------------------------------------
  subroutine hist_write(this, bounds, ubj, numf, filter, velocity)
    !
    ! DESCRIPTION
    ! output hist file, only for standalone applications
    !
    ! USES
    use bncdio_pio    , only : file_desc_t
    use bncdio_pio    , only : ncd_pio_openfile_for_write
    use bncdio_pio    , only : ncd_putvar
    use bncdio_pio    , only : ncd_pio_closefile
    use BeTR_TimeMod , only : betr_time_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: ubj
    integer                     , intent(in)    :: numf
    integer                     , intent(in)    :: filter(:)
    real(r8)                    , intent(in)    :: velocity(:, :)
    !TEMPORARY VARIABLES
    type(file_desc_t)           :: ncid
    integer                     :: jj
    integer                     :: c
    type(betr_bounds_type)     :: betr_bounds
    real(r8) :: timef
    character(len=*), parameter :: subname='hist_write'
      c = 1
    associate(                                                        &
         ntracers          => this%betr(c)%tracers%ntracers,          &
         ngwmobile_tracers => this%betr(c)%tracers%ngwmobile_tracers, &
         is_volatile       => this%betr(c)%tracers%is_volatile,       &
         is_h2o            => this%betr(c)%tracers%is_h2o,            &
         is_isotope        => this%betr(c)%tracers%is_isotope,        &
         volatileid        => this%betr(c)%tracers%volatileid,        &
         tracernames       => this%betr(c)%tracers%tracernames        &
         )
      call this%betr_time%print_model_time_stamp(iulog)

      call this%HistRetrieval(numf, filter)

      call this%BeTRSetBounds(betr_bounds)

      if(this%betr_time%its_time_to_histflush())then

        this%hist_record=this%hist_record+1

        call ncd_pio_openfile_for_write(ncid, this%hist_filename)

        timef=this%betr_time%get_cur_timef()/86400._r8;call ncd_putvar(ncid, "time", this%hist_record, timef)

        do c = bounds%begc, bounds%endc
          call ncd_putvar(ncid, 'QFLX_ADV', this%hist_record, velocity(c:c, 1:ubj))
        enddo

        call this%hist_output_states(ncid, this%hist_record, betr_bounds, numf, filter, ubj, &
            this%num_hist_state1d, this%num_hist_state2d)

        call this%hist_output_fluxes(ncid, this%hist_record, betr_bounds, numf, filter, ubj, &
           this%num_hist_flux1d, this%num_hist_flux2d)

        call ncd_pio_closefile(ncid)
      else
        call this%hist_flux_accum(betr_bounds, numf, filter, ubj, &
           this%num_hist_flux1d, this%num_hist_flux2d)
      endif

    end associate
  end subroutine hist_write

  !---------------------------------------------------------------------------------

  subroutine BeTRSimulationHistRetrieval(this,numf, filter)
  !
  !DESCRIPTION
  !Retrieve records for history files
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer, intent(in) :: numf
   integer, intent(in) :: filter(:)

  call this%BeTRRetrieveHistoryState(numf, filter)

  call this%BeTRRetrieveHistoryFlux(numf, filter)

  end subroutine BeTRSimulationHistRetrieval

  !---------------------------------------------------------------------------------
  subroutine BeTRSimulationConsistencyCheck(this, &
     bounds, ubj, num_surfc, filter_soilc, waterstate_vars)
  ! DESCRIPTION
  ! Do consistency check, can be overwritten for varies purpose
  !
  ! USES
  implicit none
   !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_surfc ! number of columns in column filter_soilc
    integer                     , intent(in)    :: filter_soilc(:) ! column filter_soilc
    integer                     , intent(in)    :: ubj
    type(Waterstate_Type)       , intent(in)    :: waterstate_vars ! water state variables

    ! remove compiler warnings
    if (this%num_surfc > 0)                         continue
    if (bounds%begc > 0)                            continue
    if (ubj > 0)                                    continue
    if (num_surfc > 0)                              continue
    if (size(filter_soilc) > 0)                     continue
    if (associated(waterstate_vars%h2osoi_liq)) continue
  end subroutine BeTRSimulationConsistencyCheck

  !---------------------------------------------------------------------------------

  subroutine WriteRegressionOutput(this, velocity)
  !DESCRIPTION
  ! output for regression test
  !
  ! USES
    use bshr_kind_mod  , only : r8 => shr_kind_r8
    use betr_constants , only : betr_string_length
    implicit none
   !ARGUMENTS
    class(betr_simulation_type) , intent(inout) :: this
    real(r8)                    , intent(in)    :: velocity(:, :)
  !TEMPORARY VARIABLES
    integer                           :: jj, tt, begc, endc, c
    character(len=betr_string_length) :: category
    character(len=betr_string_length) :: name
    integer :: loc_id

    ! FIXME(bja, 201603) should we output units as well...?
    begc = 1
    endc = 1
    c= 1
    if (this%regression%write_regression_output) then
       call this%regression%OpenOutput()
       ! NOTE(bja, 201603) currently we are allocating all tracer
       ! state vars all the time.
       do tt = 1, this%betr(c)%tracers%ntracers
          if (tt <= this%betr(c)%tracers%ngwmobile_tracers) then
             category = 'concentration'
             name = trim(this%betr(c)%tracers%tracernames(tt)) // '_total_aqueous_conc'
             call this%regression%WriteData(category, name, &
                  this%betr(c)%tracerstates%tracer_conc_mobile_col(begc, :, tt))
          end if
          if (tt <= this%betr(c)%tracers%nvolatile_tracers) then
             if(.not. this%betr(c)%tracers%is_volatile(tt))cycle
             category = 'pfraction'
             loc_id=this%betr(c)%tracers%volatileid(tt)
             name = trim(this%betr(c)%tracers%tracernames(tt)) // '_gas_partial_fraction'
             call this%regression%WriteData(category, name, &
                  this%betr(c)%tracerstates%tracer_P_gas_frac_col(begc, :, loc_id))
          end if
       end do

       name = 'total_gas_pressure'
       category = 'pressure'
       call this%regression%WriteData(category, name, &
            this%betr(c)%tracerstates%tracer_P_gas_col(begc, :))

       name = 'advective flux'
       category = 'velocity'
       call this%regression%WriteData(category, name, &
            velocity(begc, :))

       call this%regression%CloseOutput()
    end if
  end subroutine WriteRegressionOutput

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBiophysForcing(this, bounds,  col, pft, lbj, ubj, &
    carbonflux_vars, pf_carbonflux_vars, waterstate_vars,  waterflux_vars, pf_waterflux_vars, &
    temperature_vars, pf_temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
    chemstate_vars, soilstate_vars)
  !DESCRIPTION
  !pass in biogeophysical variables for running betr
  !USES
  use SoilStateType     , only : soilstate_type
  use ChemStateType     , only : chemstate_type
  use atm2lndType       , only : atm2lnd_type
  use SoilHydrologyType , only : soilhydrology_type
  use CanopyStateType   , only : canopystate_type
  use MathfuncMod       , only : isnan => bisnan
  use pftvarcon         , only : noveg
  use betr_varcon       , only : denh2o => bdenh2o, denice => bdenice
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout)        :: this
  type(bounds_type)           , intent(in)           :: bounds
  type(patch_type)            , intent(in)           :: pft
  type(column_type)           , intent(in)           :: col ! column type
  integer                     , intent(in)           :: lbj, ubj
  type(carbonflux_type)       , optional, intent(in) :: carbonflux_vars
  type(pf_carbonflux_type), optional, intent(in) :: pf_carbonflux_vars
  type(Waterstate_Type)       , optional, intent(in) :: Waterstate_vars
  type(waterflux_type)        , optional, intent(in) :: waterflux_vars
  type(pf_waterflux_type)        , optional, intent(in) :: pf_waterflux_vars
  type(temperature_type)      , optional, intent(in) :: temperature_vars
  type(pf_temperature_type)      , optional, intent(in) :: pf_temperature_vars
  type(soilhydrology_type)    , optional, intent(in) :: soilhydrology_vars
  type(atm2lnd_type)          , optional, intent(in) :: atm2lnd_vars
  type(canopystate_type)      , optional, intent(in) :: canopystate_vars
  type(chemstate_type)        , optional, intent(in) :: chemstate_vars
  type(soilstate_type)        , optional, intent(in) :: soilstate_vars

  !TEMPORARY VARIABLES
  integer :: p, pi, cc, c, l, pp
  integer :: npft_loc
  cc = 1
  do c = bounds%begc, bounds%endc
    if(.not. this%active_col(c))cycle
    this%biophys_forc(c)%stwl(cc)=0  !by default this is set zero layers of standing water
    if(present(pf_carbonflux_vars))then
      npft_loc = ubound(pf_carbonflux_vars%annsum_npp,1)-lbound(pf_carbonflux_vars%annsum_npp,1)+1
      if(col%pfti(c) /= lbound(pf_carbonflux_vars%annsum_npp,1) .and. npft_loc/=col%npfts(c))then
        do pi = 1, betr_maxpatch_pft
          this%biophys_forc(c)%annsum_npp_patch(pi) = 0._r8
          this%biophys_forc(c)%agnpp_patch(pi) = 0._r8
          this%biophys_forc(c)%bgnpp_patch(pi)  = 0._r8
        enddo
      else
        if(use_cn)then
          pp = 0
          do pi = 1, betr_maxpatch_pft
            if (pi <= col%npfts(c)) then
              p = col%pfti(c) + pi - 1
              if (pft%active(p) .and. pft%itype(p)/=noveg) then
                pp = pp + 1
                this%biophys_forc(c)%annsum_npp_patch(pp) = pf_carbonflux_vars%annsum_npp(p)
                this%biophys_forc(c)%agnpp_patch(pp)      = pf_carbonflux_vars%agnpp(p)
                this%biophys_forc(c)%bgnpp_patch(pp)      = pf_carbonflux_vars%bgnpp(p)
                this%biophys_forc(c)%tempavg_agnpp_patch(pp)= pf_carbonflux_vars%tempavg_agnpp(p)
                this%biophys_forc(c)%tempavg_bgnpp_patch(pp)= pf_carbonflux_vars%tempavg_bgnpp(p)
                this%biophys_forc(c)%annavg_agnpp_patch(pp) = pf_carbonflux_vars%annavg_agnpp(p)
                this%biophys_forc(c)%annavg_bgnpp_patch(pp) = pf_carbonflux_vars%annavg_bgnpp(p)
              endif
            endif
          enddo
        else
          npft_loc = ubound(pf_carbonflux_vars%annsum_npp,1)-lbound(pf_carbonflux_vars%annsum_npp,1)+1
          if(col%pfti(c) /= lbound(pf_carbonflux_vars%annsum_npp,1) .and. npft_loc/=col%npfts(c))then
            do pi = 1, betr_maxpatch_pft
              this%biophys_forc(c)%annsum_npp_patch(pi) = 0._r8
              this%biophys_forc(c)%agnpp_patch(pi) = 0._r8
              this%biophys_forc(c)%bgnpp_patch(pi)  = 0._r8
            enddo
          endif
        endif
      endif
    endif
    !assign waterstate
    if(present(waterstate_vars))then
      if(col%snl(c)<0)then
        this%biophys_forc(c)%h2osno_liq_col(cc,col%snl(c)+1:0) = waterstate_vars%h2osoi_liq(c,col%snl(c)+1:0)
        this%biophys_forc(c)%h2osno_ice_col(cc,col%snl(c)+1:0) = waterstate_vars%h2osoi_ice(c,col%snl(c)+1:0)
      endif
      this%biophys_forc(c)%finundated_col(cc)            = waterstate_vars%finundated(c)
      this%biophys_forc(c)%frac_h2osfc_col(cc)           = waterstate_vars%frac_h2osfc(c)
      this%biophys_forc(c)%h2osoi_liq_col(cc,lbj:ubj)    = waterstate_vars%h2osoi_liq(c,lbj:ubj)
      this%biophys_forc(c)%h2osoi_ice_col(cc,lbj:ubj)    = waterstate_vars%h2osoi_ice(c,lbj:ubj)
      this%biophys_forc(c)%h2osoi_icevol_col(cc,lbj:ubj) = waterstate_vars%h2osoi_icevol(c,lbj:ubj)
      do l = lbj, ubj
        this%biophys_forc(c)%h2osoi_icevol_col(cc,l)     = waterstate_vars%h2osoi_ice(c,l)/(col%dz(c,l)*denice)
        this%biophys_forc(c)%h2osoi_liqvol_col(cc,l)     = waterstate_vars%h2osoi_liq(c,l)/(col%dz(c,l)*denh2o)
        this%biophys_forc(c)%h2osoi_liqvol_col(cc,l)     = max(0.01_r8,this%biophys_forc(c)%h2osoi_liqvol_col(cc,l))
        this%biophys_forc(c)%h2osoi_vol_col(cc,l)        = this%biophys_forc(c)%h2osoi_liqvol_col(cc,l) + &
                                                           this%biophys_forc(c)%h2osoi_icevol_col(cc,l)
      enddo
      this%biophys_forc(c)%air_vol_col(cc,lbj:ubj)       = waterstate_vars%air_vol(c,lbj:ubj)
      this%biophys_forc(c)%smp_l_col(cc,lbj:ubj)         = waterstate_vars%smp_l(c,lbj:ubj)

!      this%biophys_forc(c)%rho_vap(cc,lbj:ubj)           = waterstate_vars%rho_vap_col(c,lbj:ubj)
!      this%biophys_forc(c)%rhvap_soi(cc,lbj:ubj)         = waterstate_vars%rhvap_soi_col(c,lbj:ubj)

    endif
    if(present(waterflux_vars))then
      this%biogeo_flux(c)%qflx_infl_col(cc)             = waterflux_vars%qflx_infl(c)
      this%biogeo_flux(c)%qflx_totdrain_col(cc)         = waterflux_vars%qflx_totdrain(c)
      this%biogeo_flux(c)%qflx_gross_evap_soil_col(cc)  = waterflux_vars%qflx_gross_evap_soil(c)
      this%biogeo_flux(c)%qflx_gross_infl_soil_col(cc)  = waterflux_vars%qflx_gross_infl_soil(c)
      this%biophys_forc(c)%qflx_surf_col(cc)            = waterflux_vars%qflx_surf(c)
      this%biophys_forc(c)%qflx_dew_grnd_col(cc)        = waterflux_vars%qflx_dew_grnd(c)
      this%biophys_forc(c)%qflx_dew_snow_col(cc)        = waterflux_vars%qflx_dew_snow(c)
      this%biophys_forc(c)%qflx_sub_snow_vol_col(cc)    = waterflux_vars%qflx_sub_snow_vol(c)
      this%biophys_forc(c)%qflx_sub_snow_col(cc)        = waterflux_vars%qflx_sub_snow(c)
      this%biophys_forc(c)%qflx_h2osfc2topsoi_col(cc)   = waterflux_vars%qflx_h2osfc2topsoi(c)
      this%biophys_forc(c)%qflx_snow2topsoi_col(cc)     = waterflux_vars%qflx_snow2topsoi(c)
      this%biophys_forc(c)%qflx_rootsoi_col(cc,lbj:ubj) = waterflux_vars%qflx_rootsoi(c,lbj:ubj)*1.e-3_r8
      this%biophys_forc(c)%qflx_runoff_col(cc)          = waterflux_vars%qflx_runoff_betr(c)  !mm/s
    endif
    if(present(pf_waterflux_vars))then
      pp = 0
      do pi = 1, betr_maxpatch_pft
       if (pi <= col%npfts(c)) then
         p = col%pfti(c) + pi - 1
         if (pft%active(p) .and. pft%itype(p)/=noveg) then
           pp = pp + 1
           this%biophys_forc(c)%qflx_tran_veg_patch(pp)     = pf_waterflux_vars%qflx_tran_veg(p)
           this%biophys_forc(c)%qflx_rootsoi_frac_patch(pp,lbj:ubj) = pf_waterflux_vars%qflx_rootsoi_frac(p,lbj:ubj)
         endif
       endif
      enddo
    endif
    if(present(temperature_vars))then
      this%biophys_forc(c)%t_soi_10cm(cc)           = temperature_vars%t_soi10cm(c)
      this%biophys_forc(c)%t_soisno_col(cc,lbj:ubj) = temperature_vars%t_soisno(c,lbj:ubj)
      if(col%snl(c)<0)then
        this%biophys_forc(c)%t_snow_col(cc,col%snl(c)+1:0) = temperature_vars%t_soisno(c,col%snl(c)+1:0)
      endif
    endif
    if(present(pf_temperature_vars))then
      pp = 0
      do pi = 1, betr_maxpatch_pft
        if (pi <= col%npfts(c)) then
          p = col%pfti(c) + pi - 1
          if (pft%active(p) .and. pft%itype(p)/=noveg) then
            pp = pp + 1
            this%biophys_forc(c)%t_veg_patch(pp)         = pf_temperature_vars%t_veg(p)
          endif
        endif
      enddo
    endif
    if(present(soilhydrology_vars))then
      this%biophys_forc(c)%qflx_bot_col(cc)        = soilhydrology_vars%qcharge_col(c)
      this%biophys_forc(c)%fracice_col(cc,lbj:ubj) = soilhydrology_vars%fracice_col(c,lbj:ubj)
    endif

    if(present(atm2lnd_vars))then
      this%biophys_forc(c)%forc_pbot_downscaled_col(cc) = atm2lnd_vars%forc_pbot_downscaled_col(c)
      this%biophys_forc(c)%forc_t_downscaled_col(cc)    = atm2lnd_vars%forc_t_downscaled_col(c)
    endif

    if(present(canopystate_vars))then
      this%biophys_forc(c)%altmax_col(cc)          = canopystate_vars%altmax_col(c)
      this%biophys_forc(c)%altmax_lastyear_col(cc) = canopystate_vars%altmax_lastyear_col(c)
      pp = 0
      do pi = 1, betr_maxpatch_pft
        if (pi <= col%npfts(c)) then
          p = col%pfti(c) + pi - 1
          if (pft%active(p) .and. pft%itype(p)/=noveg) then
            pp = pp + 1
            this%biophys_forc(c)%lbl_rsc_h2o_patch(pp) = canopystate_vars%lbl_rsc_h2o_patch(p)
            this%biophys_forc(c)%elai_patch(pp)        = canopystate_vars%elai_patch(p)
          endif
        endif
      enddo
    endif
    if(present(chemstate_vars))then
      this%biophys_forc(c)%soil_pH(cc,lbj:ubj) = chemstate_vars%soil_pH(c,lbj:ubj)
    endif
    if(present(soilstate_vars))then
      this%biophys_forc(c)%bsw_col(cc,lbj:ubj)          = soilstate_vars%bsw_col(c,lbj:ubj)
      this%biophys_forc(c)%watsat_col(cc,lbj:ubj)       = soilstate_vars%watsat_col(c,lbj:ubj)
      this%biophys_forc(c)%eff_porosity_col(cc,lbj:ubj) = soilstate_vars%eff_porosity_col(c,lbj:ubj)
      this%biophys_forc(c)%cellorg_col(cc,lbj:ubj)      = soilstate_vars%cellorg_col(c,lbj:ubj)
      this%biophys_forc(c)%cellclay_col(cc,lbj:ubj)     = soilstate_vars%cellclay_col(c,lbj:ubj)
      this%biophys_forc(c)%cellsand_col(cc,lbj:ubj)     = soilstate_vars%cellsand_col(c,lbj:ubj)
      this%biophys_forc(c)%bd_col(cc,lbj:ubj)           = soilstate_vars%bd_col(c,lbj:ubj)
      this%biophys_forc(c)%watfc_col(cc,lbj:ubj)        = soilstate_vars%watfc_col(c,lbj:ubj)
      this%biophys_forc(c)%sucsat_col(cc,lbj:ubj)       = soilstate_vars%sucsat_col(c,lbj:ubj)
      pp = 0
      do pi = 1, betr_maxpatch_pft
        if (pi <= col%npfts(c)) then
          p = col%pfti(c) + pi - 1
          if (pft%active(p) .and. pft%itype(p)/=noveg) then
            pp = pp + 1
            this%biophys_forc(c)%rootfr_patch(pp,lbj:ubj) = soilstate_vars%rootfr_patch(p,lbj:ubj)
          endif
        endif
      enddo
    endif
  enddo
  end subroutine BeTRSimulationSetBiophysForcing

  !------------------------------------------------------------------------
  subroutine BeTRSimulationRetrieveBiogeoFlux(this, bounds, lbj,ubj, carbonflux_vars,  &
    waterflux_vars)
  ! DESCRIPTIONS
  ! update and return fluxes, this eventually will be expanded to
  ! include other fluxes
  ! USES
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout)           :: this
  type(bounds_type)           , intent(in)              :: bounds
  integer                     , intent(in)              :: lbj, ubj
  type(carbonflux_type)       , optional, intent(inout) :: carbonflux_vars
  type(waterflux_type)        , optional, intent(inout) :: waterflux_vars

  integer :: begp, begc, endp, endc
  integer :: p, c, cc
  cc = 1
  if(present(carbonflux_vars))then
    !do nothing
  endif
  if(present(waterflux_vars))then
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      waterflux_vars%qflx_infl(c)            = this%biogeo_flux(c)%qflx_infl_col(cc)
      waterflux_vars%qflx_adv(c,lbj-1:ubj)   = this%biogeo_flux(c)%qflx_adv_col(cc,lbj-1:ubj)
      waterflux_vars%qflx_totdrain(c)        = this%biogeo_flux(c)%qflx_totdrain_col(cc)
      waterflux_vars%qflx_gross_evap_soil(c) = this%biogeo_flux(c)%qflx_gross_evap_soil_col(cc)
      waterflux_vars%qflx_gross_infl_soil(c) = this%biogeo_flux(c)%qflx_gross_infl_soil_col(cc)
      waterflux_vars%qflx_drain_vr(c,1:ubj)  = this%biogeo_flux(c)%qflx_drain_vr_col(cc,1:ubj)
    enddo
  endif

  end subroutine BeTRSimulationRetrieveBiogeoFlux

  !------------------------------------------------------------------------
  subroutine BeTRSimulationPreDiagSoilColWaterFlux(this, num_nolakec, filter_nolakec)
  !DESCRIPTION
  !prepare for water flux diagnosis. it is called before diagnosing the advective fluxes and applying
  !freeze-thaw tracer partition.
  !
  !USES

  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
   integer                    , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
   integer                    , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points

  !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c


   call this%BeTRSetBounds(betr_bounds)

   do fc= 1, num_nolakec
     c = filter_nolakec(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%pre_diagnose_soilcol_water_flux(betr_bounds, this%num_surfc, &
       this%filter_soilc, this%biophys_forc(c))
   enddo
  end subroutine BeTRSimulationPreDiagSoilColWaterFlux

  !------------------------------------------------------------------------
  subroutine BeTRSimulationDiagnoseDtracerFreezeThaw(this, bounds, num_nolakec, filter_nolakec, col, lun)
  !
  ! DESCRIPTION
  ! aqueous tracer partition based on freeze-thaw
  !
  ! USES
  implicit none
  !
  ! Arguments
  class(betr_simulation_type), intent(inout)   :: this
  type(bounds_type)     , intent(in) :: bounds
  integer               , intent(in) :: num_nolakec                        ! number of column non-lake points in column filter
  integer               , intent(in) :: filter_nolakec(:)                  ! column filter for non-lake points
!  type(waterstate_type), intent(in) :: waterstate_vars
  type(column_type)     , intent(in) :: col                                ! column type
  type(landunit_type)   , intent(in)  :: lun

  !temporary variables
  type(betr_bounds_type)     :: betr_bounds
  integer :: fc, c

  call this%BeTRSetBounds(betr_bounds)

  call this%BeTRSetcps(bounds, col)

  do fc = 1, num_nolakec
    c = filter_nolakec(fc)
    if(.not. this%active_col(c))cycle
    call this%betr(c)%diagnose_dtracer_freeze_thaw(betr_bounds, this%num_surfc, this%filter_soilc,  &
      this%biophys_forc(c))
  enddo
  end subroutine BeTRSimulationDiagnoseDtracerFreezeThaw

  !------------------------------------------------------------------------
  subroutine BeTRSimulationDiagAdvWaterFlux(this, num_hydrologyc, &
    filter_hydrologyc)

  !DESCRIPTION
  ! diagnose water fluxes for tracer advection
  !
  ! USES
  !
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_hydrologyc                        ! number of column non-lake points in column filter
   integer                     , intent(in)    :: filter_hydrologyc(:)                  ! column filter for non-lake points

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c, j

   call this%BeTRSetBounds(betr_bounds)
   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%diagnose_advect_water_flux(this%betr_time,              &
       betr_bounds, this%num_surfc, this%filter_soilc,                         &
       this%biophys_forc(c), this%biogeo_flux(c))
   enddo

  end subroutine BeTRSimulationDiagAdvWaterFlux
  !------------------------------------------------------------------------
  subroutine BeTRSimulationDiagDrainWaterFlux(this, num_hydrologyc, filter_hydrologyc)
  !DESCRIPTION
  ! diagnose water fluxes due to subsurface drainage
  !
  ! USES
  !
    use SoilHydrologyType , only : soilhydrology_type
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_hydrologyc                        ! number of column non-lake points in column filter
   integer                     , intent(in)    :: filter_hydrologyc(:)                  ! column filter for non-lake points

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%diagnose_drainage_water_flux(this%betr_time, &
       betr_bounds, this%num_surfc, this%filter_soilc,      &
       this%biophys_forc(c), this%biogeo_flux(c))
  enddo
  end subroutine BeTRSimulationDiagDrainWaterFlux
  !------------------------------------------------------------------------
  subroutine BeTRSimulationBeginTracerSnowLayerAdjst(this,  num_snowc, filter_snowc)
  !DESCRIPTION
  !prepare for tracer adjustment in snow layers
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_snowc
   integer                     , intent(in)    :: filter_snowc(:)

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%Enter_tracer_LayerAdjustment(betr_bounds, this%betr_col(c), &
       this%num_surfc, this%filter_soilc)
   enddo
  end subroutine BeTRSimulationBeginTracerSnowLayerAdjst
  !------------------------------------------------------------------------
  subroutine BeTRSimulationEndTracerSnowLayerAdjst(this, num_snowc, filter_snowc)
  !DESCRIPTION
  !wrap up tracer adjustment in snow layers
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   integer                     , intent(in)    :: num_snowc
   integer                     , intent(in)    :: filter_snowc(:)

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%Exit_tracer_LayerAdjustment(betr_bounds, this%betr_col(c), &
       this%num_surfc, this%filter_soilc)
   enddo

  end subroutine BeTRSimulationEndTracerSnowLayerAdjst
  !------------------------------------------------------------------------
  subroutine BeTRSimulationDvideSnowLayers(this, bounds, num_snowc, filter_snowc, divide_matrix)
  !DESCRIPTIONS
  !redistribute tracer in snow layers due to division
  !
  !USES
  use elm_varpar, only : nlevsno
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds               ! bounds
   integer                     , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                     , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)                    , intent(in)    :: divide_matrix(bounds%begc:bounds%endc , 1:nlevsno , 1:nlevsno )

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%tracer_DivideSnowLayers(betr_bounds, this%betr_col(c),this%num_surfc, &
       this%filter_soilc, divide_matrix(c:c,:,:), this%bstatus(c))
     if(this%bstatus(c)%check_status())then
       call this%bsimstatus%setcol(c)
       call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
       exit
     endif
   enddo
  if(this%bsimstatus%check_status()) &
    call endrun(msg=trim(this%bsimstatus%print_msg()))
  end subroutine BeTRSimulationDvideSnowLayers

  !------------------------------------------------------------------------
  subroutine BeTRSimulationCombineSnowLayers(this, bounds, num_snowc, filter_snowc, combine_matrix)
  !DESCRIPTIONS
  !redistribute tracer in snow layers due to division
  !
  !USES
  use elm_varpar, only : nlevsno
  use betr_varcon, only : betr_maxpatch_pft
  implicit none
  !ARGUMENTS
   class(betr_simulation_type) , intent(inout) :: this
   type(bounds_type)           , intent(in)    :: bounds               ! bounds
   integer                     , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                     , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)                    , intent(in)    :: combine_matrix(bounds%begc:bounds%endc,-nlevsno+1:1 ,-nlevsno+1:1 )

   !TEMPORARY VARIABLES
   type(betr_bounds_type)     :: betr_bounds
   integer :: fc, c

   call this%BeTRSetBounds(betr_bounds)

   do fc = 1, num_snowc
     c = filter_snowc(fc)
     if(.not. this%active_col(c))cycle
     call this%betr(c)%tracer_CombineSnowLayers(betr_bounds, this%betr_col(c),this%num_surfc,&
       this%filter_soilc, combine_matrix(c:c,:,:),this%bstatus(c))
     if(this%bstatus(c)%check_status())then
       call this%bsimstatus%setcol(c)
       call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
       exit
     endif
   enddo
   if(this%bsimstatus%check_status()) &
    call endrun(msg=trim(this%bsimstatus%print_msg()))
  end subroutine BeTRSimulationCombineSnowLayers

  !------------------------------------------------------------------------
  subroutine hist_create_fluxes(this, bounds, betr_nlevtrc_soil, num_flux1d, num_flux2d, ncid)
  !
  !DESCRIPTION
  !create history file for betr fluxes
  !
  use histFileMod         , only: hist_addfld1d, hist_addfld2d
  use bhistFileMod        , only : hist_def_fld2d , hist_def_fld1d
  use bncdio_pio           , only : file_desc_t, ncd_float
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: betr_nlevtrc_soil
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d
  type(file_desc_t) ,   optional,  intent(inout)   :: ncid
  !local variables
  integer :: jj, begc, endc

  character(len=*), parameter :: subname = 'hist_create_fluxes'

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  if(betr_offline .and. (.not. present(ncid)))then
    call endrun(msg="ncid not defined in "//subname//errmsg(mod_filename, __LINE__))
  endif

  begc=bounds%begc; endc=bounds%endc

  do jj = 1, num_flux2d

    if(betr_offline)then
      if(trim(this%flux_hist2d_var(jj)%use_default)/='inactive')&
        call hist_def_fld2d (ncid, varname=this%flux_hist2d_var(jj)%varname, &
            nf90_type=ncd_float, dim1name = "ncol",&
            dim2name="levtrc", long_name=this%flux_hist2d_var(jj)%long_name, &
            units=this%flux_hist2d_var(jj)%units)
    else
      this%hist_fluxes_2d(begc:endc,1:betr_nlevtrc_soil, jj) = spval
      data2dptr => this%hist_fluxes_2d(begc:endc,1:betr_nlevtrc_soil, jj)
      call hist_addfld2d (fname=this%flux_hist2d_var(jj)%varname, &
           units=this%flux_hist1d_var(jj)%units, type2d="levtrc",  &
           avgflag=this%flux_hist2d_var(jj)%avg_flag, &
           long_name=this%flux_hist2d_var(jj)%long_name,  ptr_col=data2dptr, &
           default=this%flux_hist2d_var(jj)%use_default)
    endif
  enddo

  do jj = 1, num_flux1d

    if(betr_offline)then
      if(trim(this%flux_hist1d_var(jj)%use_default)/='inactive')&
      call hist_def_fld1d (ncid, varname=this%flux_hist1d_var(jj)%varname, &
        nf90_type=ncd_float, &
        dim1name="ncol", long_name=this%flux_hist1d_var(jj)%long_name,&
        units=this%flux_hist1d_var(jj)%units)
    else
      this%hist_fluxes_1d(begc:endc,jj) = spval
      data1dptr => this%hist_fluxes_1d(begc:endc, jj)
      call hist_addfld1d (fname=this%flux_hist1d_var(jj)%varname, &
        units=this%flux_hist1d_var(jj)%units,  &
        avgflag=this%flux_hist1d_var(jj)%avg_flag, &
        long_name=this%flux_hist1d_var(jj)%long_name, &
        ptr_col=data1dptr, default=this%flux_hist1d_var(jj)%use_default)
    endif

  enddo

  end subroutine hist_create_fluxes

  !------------------------------------------------------------------------
  subroutine hist_create_states(this, bounds, betr_nlevtrc_soil, num_state1d, num_state2d, ncid)
  !
  !create history file for betr states variables
  use histFileMod   , only: hist_addfld1d, hist_addfld2d
  use bhistFileMod        , only : hist_def_fld2d , hist_def_fld1d
  use bncdio_pio , only : file_desc_t, ncd_float
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer, intent(in) :: betr_nlevtrc_soil
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  type(file_desc_t) ,   optional,  intent(inout)   :: ncid
  !local variables
  integer :: begc, endc
  integer :: jj

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  character(len=*), parameter :: subname = 'hist_create_states'


  if(betr_offline .and. (.not. present(ncid)))then
    call endrun(msg="ncid not defined in "//subname//errmsg(mod_filename, __LINE__))
  endif

  begc = bounds%begc; endc = bounds%endc

  do jj = 1, num_state2d
    !read namelist

    if(betr_offline)then
      if(trim(this%state_hist2d_var(jj)%use_default)/='inactive') &
        call hist_def_fld2d (ncid=ncid, varname=this%state_hist2d_var(jj)%varname, &
          nf90_type=ncd_float, dim1name = "ncol",&
          dim2name="levtrc", long_name=this%state_hist2d_var(jj)%long_name,&
          units=this%state_hist2d_var(jj)%units)
    else
      this%hist_states_2d(begc:endc,1:betr_nlevtrc_soil, jj) = spval
      data2dptr => this%hist_states_2d(begc:endc,1:betr_nlevtrc_soil, jj)
      call hist_addfld2d (fname=this%state_hist2d_var(jj)%varname, &
           units=this%state_hist2d_var(jj)%units, type2d="levtrc",  &
           avgflag=this%state_hist2d_var(jj)%avg_flag, &
           long_name=this%state_hist2d_var(jj)%long_name, &
           ptr_col=data2dptr, default=this%state_hist2d_var(jj)%use_default)
    endif
  enddo

  do jj = 1, num_state1d

    if(betr_offline)then
      if(trim(this%state_hist1d_var(jj)%use_default)/='inactive') &
        call hist_def_fld1d (ncid, varname=this%state_hist1d_var(jj)%varname,  nf90_type=ncd_float, &
          dim1name="ncol", long_name=this%state_hist1d_var(jj)%long_name, units=this%state_hist1d_var(jj)%units)
    else
      this%hist_states_1d(begc:endc,jj) = spval
      data1dptr => this%hist_states_1d(begc:endc,jj)
      call hist_addfld1d (fname=this%state_hist1d_var(jj)%varname, &
          units=this%state_hist1d_var(jj)%units,      &
          avgflag=this%state_hist1d_var(jj)%avg_flag, &
          long_name=this%state_hist1d_var(jj)%long_name, &
          ptr_col=data1dptr, default=this%state_hist1d_var(jj)%use_default)
    endif
  enddo

  end subroutine hist_create_states
  !------------------------------------------------------------------------

  subroutine hist_output_fluxes(this,  ncid, record, bounds, numf, filter, &
     betr_nlevtrc_soil, num_flux1d, num_flux2d)
  !
  !DESCRIPTION
  !create history file for betr fluxes
  !
  use histFileMod         , only: hist_addfld1d, hist_addfld2d
  use bncdio_pio           , only :   file_desc_t, ncd_putvar
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: betr_nlevtrc_soil
  integer, intent(in) :: record
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)
  type(betr_bounds_type)           , intent(in)    :: bounds               ! bounds
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d
  type(file_desc_t) ,     intent(inout)   :: ncid
  !local variables
  integer :: jj, begc, endc, jl, fc, c

  character(len=*), parameter :: subname = 'hist_output_fluxes'

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays


  begc=bounds%begc; endc=bounds%endc
  this%hist_naccum = this%hist_naccum+1._r8

  do jj = 1, num_flux2d
    do jl = 1, betr_nlevtrc_soil
      do fc = 1, numf
        c =filter(fc)
        this%hist_fluxes_2d_accum(c,jl,jj)=this%hist_fluxes_2d_accum(c,jl,jj)/this%hist_naccum
      enddo
    enddo
    if(trim(this%flux_hist2d_var(jj)%use_default)/='inactive') then
      data2dptr => this%hist_fluxes_2d_accum(begc:endc,1:betr_nlevtrc_soil, jj)
      call ncd_putvar(ncid, this%flux_hist2d_var(jj)%varname, record, data2dptr)
    endif
  enddo

  do jj = 1, num_flux1d
    do fc = 1, numf
      c =filter(fc)
      this%hist_fluxes_1d_accum(c,jj) = this%hist_fluxes_1d_accum(c,jj)/this%hist_naccum
    enddo
    if(trim(this%flux_hist1d_var(jj)%use_default)/='inactive')then
      data1dptr => this%hist_fluxes_1d_accum(begc:endc, jj)
      call ncd_putvar(ncid,this%flux_hist1d_var(jj)%varname, record, data1dptr)
    endif
  enddo

  this%hist_naccum = 0._r8
  do jj = 1, num_flux2d
    do jl = 1, betr_nlevtrc_soil
      do fc = 1, numf
        c =filter(fc)
        this%hist_fluxes_2d_accum(c,jl,jj) = 0._r8
      enddo
    enddo
  enddo

  do jj = 1, num_flux1d
    do fc = 1, numf
      c =filter(fc)
      this%hist_fluxes_1d_accum(c,jj) = 0._r8
    enddo
  enddo
  end subroutine hist_output_fluxes

  !------------------------------------------------------------------------
  subroutine hist_flux_accum(this, bounds, numf, filter, &
     betr_nlevtrc_soil, num_flux1d, num_flux2d)
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  type(betr_bounds_type)      , intent(in)    :: bounds               ! bounds
  integer                     , intent(in)    :: betr_nlevtrc_soil
  integer                     , intent(in)    :: numf
  integer                     , intent(in)    :: filter(:)
  integer                     , intent(in)    :: num_flux1d
  integer                     , intent(in)    :: num_flux2d

  integer :: jj, jl, fc, c

  do jj = 1, num_flux2d
    do jl = 1, betr_nlevtrc_soil
      do fc = 1, numf
        c =filter(fc)
        this%hist_fluxes_2d_accum(c, jl, jj) = &
          this%hist_fluxes_2d_accum(c, jl, jj) + this%hist_fluxes_2d(c, jl, jj)
      enddo
    enddo
  enddo

  do jj = 1, num_flux1d
    do fc = 1, numf
      c =filter(fc)
      this%hist_fluxes_1d_accum(c, jj) = this%hist_fluxes_1d_accum(c, jj) + &
         this%hist_fluxes_1d(c, jj)
    enddo
  enddo

  this%hist_naccum=this%hist_naccum + 1._r8
  end subroutine hist_flux_accum

  !------------------------------------------------------------------------
  subroutine hist_output_states(this,  ncid,  record, bounds, numf, filter, &
     betr_nlevtrc_soil, num_state1d, num_state2d)
  !
  !create history file for betr states variables
  use histFileMod   , only: hist_addfld1d, hist_addfld2d
  use bncdio_pio , only : file_desc_t, ncd_putvar
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: record
  type(betr_bounds_type)           , intent(in)    :: bounds               ! bounds
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)
  integer, intent(in) :: betr_nlevtrc_soil
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  type(file_desc_t) ,     intent(inout)   :: ncid
  !local variables
  integer :: begc, endc
  integer :: jj

  real(r8), pointer :: data2dptr(:,:) ! temp. pointers for slicing larger arrays
  real(r8), pointer :: data1dptr(:)   ! temp. pointers for slicing larger arrays

  character(len=*), parameter :: subname = 'hist_output_states'

  begc = bounds%begc; endc = bounds%endc

  do jj = 1, num_state2d
    if(trim(this%state_hist2d_var(jj)%use_default)/='inactive') then
      data2dptr => this%hist_states_2d(begc:endc,1:betr_nlevtrc_soil, jj)
      call ncd_putvar(ncid,this%state_hist2d_var(jj)%varname, record, data2dptr)
    endif
  enddo

  do jj = 1, num_state1d
    if(trim(this%state_hist1d_var(jj)%use_default)/='inactive') then
      data1dptr => this%hist_states_1d(begc:endc,jj)
      call ncd_putvar(ncid,this%state_hist1d_var(jj)%varname, record, data1dptr)
    endif
  enddo

  end subroutine hist_output_states
  !------------------------------------------------------------------------
  subroutine BeTRSimulationCreateHistory(this, bounds, betr_nlevtrc_soil,&
     num_state1d, num_state2d, num_flux1d, num_flux2d)
  !
  !links the variable to output
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)           , intent(in)    :: bounds               ! bounds
  integer           , intent(in) :: betr_nlevtrc_soil
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d

  call this%hist_create_states(bounds, betr_nlevtrc_soil, num_state1d, num_state2d)

  call this%hist_create_fluxes(bounds, betr_nlevtrc_soil, num_flux1d, num_flux2d)

  end subroutine BeTRSimulationCreateHistory
  !------------------------------------------------------------------------
  subroutine BeTRSimulationRetrieveHistoryState(this, numf, filter)
  !
  !DESCRIPTION
  !retrieve state variables for history writing
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !TEMPORARY VARIABLES
  integer :: fc, c
  type(betr_bounds_type)     :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)

  do fc = 1, numf
    c = filter(fc)
    if(.not. this%active_col(c))cycle
    call this%betr(c)%HistRetrieveState(betr_bounds, 1, betr_bounds%ubj, &
       this%num_hist_state1d, this%num_hist_state2d,&
       this%hist_states_1d(c:c,:), this%hist_states_2d(c:c,1:betr_bounds%ubj,:))
  enddo

  end subroutine BeTRSimulationRetrieveHistoryState
  !------------------------------------------------------------------------
  subroutine BeTRSimulationRetrieveHistoryFlux(this, numf, filter)
  use tracer_varcon  , only :  betr_nlevtrc_soil
  implicit none
  !ARGUMENTS
  class(betr_simulation_type) , intent(inout) :: this
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !TEMPORARY VARIABLES
  integer :: fc, c
  type(betr_bounds_type)     :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)
  this%hist_fluxes_1d(:,:)=spval
  this%hist_fluxes_2d(:,:,:)=spval
  do fc = 1, numf
    c = filter(fc)
    if(.not. this%active_col(c))cycle
    call this%betr(c)%HistRetrieveFlux(betr_bounds, 1, betr_bounds%ubj, &
       this%num_hist_flux1d,this%num_hist_flux2d, &
       this%hist_fluxes_1d(c:c,:),this%hist_fluxes_2d(c:c,1:betr_bounds%ubj,:))
  enddo

  end subroutine BeTRSimulationRetrieveHistoryFlux
  !------------------------------------------------------------------------
  subroutine BeTRSimulationRestartOffline(this, bounds, ncid, numf, filter, flag)
  !DESCRIPTION
  !create or read restart file
  use restUtilMod    , only : restartvar
  use bncdio_pio      , only : file_desc_t,ncd_double
  use bncdio_pio      , only : ncd_defvar, ncd_defdim
  use bncdio_pio      , only : ncd_enddef, ncd_putvar
  use bncdio_pio      , only : ncd_getvar, ncd_int
  implicit none
  ! !ARGUMENTS:
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)    , intent(in)    :: bounds
  type(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
  character(len=*)     , intent(in)    :: flag ! 'read' or 'write'
  integer, intent(in) :: numf
  integer, intent(in) :: filter(:)

  !local variables

  integer :: c, jj, fc
  character(len=255),allocatable :: rest_varname_1d(:)
  character(len=255),allocatable :: rest_varname_2d(:)
  logical :: readvar      ! determine if variable is on initial file
  real(r8), pointer :: ptr1d(:)
  real(r8), pointer :: ptr2d(:,:)
  type(betr_bounds_type)     :: betr_bounds
  integer :: recordDimID

  allocate(rest_varname_1d(this%num_rest_state1d)); rest_varname_1d=''
  allocate(rest_varname_2d(this%num_rest_state2d)); rest_varname_2d=''

  c = bounds%begc
  call this%betr(c)%get_restartvar_info(this%num_rest_state1d, &
    this%num_rest_state2d,rest_varname_1d, rest_varname_2d)

  !assign initial conditions
  call this%BeTRSetBounds(betr_bounds)

  print*,'offline restart ', flag
  if(flag=='define')then
    ! print*,'define restart file'
    ! define the dimensions
    !the temporal dimension is infinite

    !number of vertical layers
    call ncd_defdim(ncid, 'levtrc', betr_nlevsoi, recordDimID)

    !number of columns
    call ncd_defdim(ncid, 'column', numf, recordDimID)

    call ncd_defvar(ncid, 'hist_naccum',ncd_double, long_name='', &
         units = '',  missing_value=spval, fill_value=spval)

    do jj = 1, this%num_rest_state1d
        !x print*,jj,trim(rest_varname_1d(jj))
      call ncd_defvar(ncid, trim(rest_varname_1d(jj)),ncd_double,dim1name='column',  &
          long_name='', units = '',  missing_value=spval, fill_value=spval)
    enddo

    do jj =1, this%num_rest_state2d
      !x print*,jj,trim(rest_varname_2d(jj))
      call ncd_defvar(ncid, trim(rest_varname_2d(jj)),ncd_double,dim1name='column',  &
        dim2name='levtrc', long_name='', units = '',  missing_value=spval, fill_value=spval)
    enddo

    do jj = 1, this%num_hist_flux1d
       call ncd_defvar(ncid, trim(this%flux_hist1d_var(jj)%varname)//'_accum',ncd_double,dim1name='column',  &
           long_name='', units = '',  missing_value=spval, fill_value=spval)
    enddo

    do jj = 1, this%num_hist_flux2d
      call ncd_defvar(ncid, trim(this%flux_hist2d_var(jj)%varname)//'_accum',ncd_double,dim1name='column',  &
        dim2name='levtrc', long_name='', units = '',  missing_value=spval, fill_value=spval)
    enddo
    !define time information
    call ncd_defvar(ncid, 'tod',ncd_double,long_name='time of day',units='',missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid, 'toy',ncd_double,long_name='time of year',units='',missing_value=spval, fill_value=spval)
    call ncd_defvar(ncid, 'dow',ncd_int, long_name='day of week', imissing_value=ispval, ifill_value=ispval)
    call ncd_defvar(ncid, 'dom',ncd_int, long_name='day of month', imissing_value=ispval, ifill_value=ispval)
    call ncd_defvar(ncid, 'doy',ncd_int, long_name='day of year', imissing_value=ispval, ifill_value=ispval)
    call ncd_defvar(ncid, 'moy',ncd_int, long_name='month of year', imissing_value=ispval, ifill_value=ispval)
    call ncd_defvar(ncid, 'cyears',ncd_int, long_name='cumulative years', imissing_value=ispval, ifill_value=ispval)
    call ncd_defvar(ncid, 'cdays',ncd_int, long_name='cumulative days', imissing_value=ispval, ifill_value=ispval)
    call ncd_defvar(ncid, 'tstep',ncd_int, long_name='steps of the year', imissing_value=ispval, ifill_value=ispval)

    call ncd_enddef(ncid)

  elseif(flag=='write')then

    do fc = 1, numf
      c = filter(fc)
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
        this%num_rest_state1d, this%num_rest_state2d, &
        this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo

    print*,'write restart file'
    do jj = 1, this%num_rest_state1d
       ptr1d => this%rest_states_1d(:, jj)
       call ncd_putvar(ncid, trim(rest_varname_1d(jj)), 1, ptr1d)
    enddo

    do jj = 1, this%num_rest_state2d
      ptr2d => this%rest_states_2d(:, :, jj)
      call ncd_putvar(ncid, trim(rest_varname_2d(jj)), 1, ptr2d)
    enddo

    do jj = 1, this%num_hist_flux1d
       ptr1d => this%hist_fluxes_1d_accum(:, jj)
       call ncd_putvar(ncid, trim(this%flux_hist1d_var(jj)%varname)//'_accum', 1, ptr1d)
    enddo

    do jj = 1, this%num_hist_flux2d
      ptr2d => this%hist_fluxes_2d_accum(:,:,jj)
      call ncd_putvar(ncid, trim(this%flux_hist2d_var(jj)%varname)//'_accum', 1, ptr2d)
    enddo

    call ncd_putvar(ncid,'hist_naccum',this%hist_naccum)

    call ncd_putvar(ncid,'tod',this%betr_time%tod)
    call ncd_putvar(ncid,'toy',this%betr_time%toy)
    call ncd_putvar(ncid,'dow',this%betr_time%dow)
    call ncd_putvar(ncid,'dom',this%betr_time%dom)
    call ncd_putvar(ncid,'doy',this%betr_time%doy)
    call ncd_putvar(ncid,'moy',this%betr_time%moy)
    call ncd_putvar(ncid,'cyears',this%betr_time%cyears)
    call ncd_putvar(ncid,'cdays',this%betr_time%cdays)
    call ncd_putvar(ncid,'tstep',this%betr_time%tstep)

  elseif(flag=='read')then
      ! print*,'read restart file'
    do jj = 1, this%num_rest_state1d
       ptr1d => this%rest_states_1d(:, jj)
       call ncd_getvar(ncid, trim(rest_varname_1d(jj)), ptr1d)
    enddo

    do jj = 1, this%num_rest_state2d
      ptr2d => this%rest_states_2d(:, :, jj)
      call ncd_getvar(ncid, trim(rest_varname_2d(jj)), ptr2d)
    enddo

    do jj = 1, this%num_hist_flux1d
       ptr1d => this%hist_fluxes_1d_accum(:, jj)
       call ncd_getvar(ncid, trim(this%flux_hist1d_var(jj)%varname)//'_accum', ptr1d)
    enddo

    do jj = 1, this%num_hist_flux2d
      ptr2d => this%hist_fluxes_2d_accum(:,:,jj)
      call ncd_getvar(ncid, trim(this%flux_hist2d_var(jj)%varname)//'_accum', ptr2d)
    enddo

    call ncd_getvar(ncid,'hist_naccum',this%hist_naccum)

    ! print*,'assign values to state variables',flag
    do fc = 1, numf
      c = filter(fc)
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
      this%num_rest_state1d,this%num_rest_state2d, &
      this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo

    call ncd_getvar(ncid,'tod',this%betr_time%tod)
    call ncd_getvar(ncid,'toy',this%betr_time%toy)
    call ncd_getvar(ncid,'dow',this%betr_time%dow)
    call ncd_getvar(ncid,'dom',this%betr_time%dom)
    call ncd_getvar(ncid,'doy',this%betr_time%doy)
    call ncd_getvar(ncid,'moy',this%betr_time%moy)
    call ncd_getvar(ncid,'cyears',this%betr_time%cyears)
    call ncd_getvar(ncid,'cdays',this%betr_time%cdays)
    call ncd_getvar(ncid,'tstep',this%betr_time%tstep)
  endif

  deallocate(rest_varname_1d)
  deallocate(rest_varname_2d)

  end subroutine BeTRSimulationRestartOffline

  !------------------------------------------------------------------------
  subroutine BeTRSimulationRestart(this, bounds, ncid, flag)
  !DESCRIPTION
  !create or read restart file
  use restUtilMod    , only : restartvar
  use ncdio_pio      , only : file_desc_t,ncd_double, ncd_int
  use elm_varctl     , only : spinup_state
  use clm_time_manager, only : get_nstep
  use betr_ctrl      , only : exit_spinup, enter_spinup,betr_spinup_state
  use tracer_varcon  , only : reaction_method
  implicit none
  ! !ARGUMENTS:
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type)    , intent(in)    :: bounds
  type(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
  character(len=*)     , intent(in)    :: flag ! 'read' or 'write'

  !local variables
  integer :: nrest_1d, nrest_2d
  integer :: c, jj, fc
  character(len=255), allocatable :: rest_varname_1d(:)
  character(len=255), allocatable :: rest_varname_2d(:)
  logical :: readvar      ! determine if variable is on initial file
  real(r8), pointer :: ptr1d(:)
  real(r8), pointer :: ptr2d(:,:)
  type(betr_bounds_type)     :: betr_bounds
  integer :: recordDimID
  integer  :: idata
  integer  :: restart_file_spinup_state
  integer  :: c_l

  c_l = 1
  c = bounds%begc
  restart_file_spinup_state =0
  allocate(rest_varname_1d(this%num_rest_state1d)); rest_varname_1d=''
  allocate(rest_varname_2d(this%num_rest_state2d)); rest_varname_2d=''

  c = bounds%begc
  call this%betr(c)%get_restartvar_info(this%num_rest_state1d, &
    this%num_rest_state2d,rest_varname_1d, rest_varname_2d)

  call this%BeTRSetBounds(betr_bounds)
  if(trim(flag)=='write')then
    if(this%do_soibgc())then
      idata = spinup_state
      do c = bounds%begc, bounds%endc
        this%scalaravg_col(c) = this%biophys_forc(c)%scalaravg_col(c_l)
        this%dom_scalar_col(c) = this%biophys_forc(c)%dom_scalar_col(c_l)
      enddo
    endif
    do c = bounds%begc, bounds%endc
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
        this%num_rest_state1d,this%num_rest_state2d, &
        this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo
  endif

  if(this%do_soibgc())then
    call restartvar(ncid=ncid, flag=flag, varname='spinscalar', xtype=ncd_double, &
         dim1name='column', long_name='', units='', &
         interpinic_flag = 'interp', readvar=readvar, data=this%scalaravg_col)

     call restartvar(ncid=ncid, flag=flag, varname='domspinscalar', xtype=ncd_double, &
         dim1name='column', long_name='', units='', &
         interpinic_flag = 'interp', readvar=readvar, data=this%dom_scalar_col)

    call restartvar(ncid=ncid, flag=flag, varname='betr_spinup_state', xtype=ncd_int,  &
           long_name='Spinup state of betr model that wrote this restart file: ' &
           // ' 0,1,2=not ready for spinup scalar, 3 = apply spinup scalar', units='', &
           interpinic_flag='copy', readvar=readvar,  data=betr_spinup_state)

    call restartvar(ncid=ncid, flag=flag, varname='spinup_count', xtype=ncd_int,  &
           long_name='Spinup count of the model that wrote this restart file: ' &
           // ' 0 <=2 skip mass bal check, 3 = do mass bal check', units='', &
           interpinic_flag='copy', readvar=readvar,  data=this%spinup_count)

    if(trim(flag)=='read')then
      call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
             long_name='Spinup state of the model that wrote this restart file: ' &
             // ' 0 = normal model mode, 1 = AD spinup', units='', &
             interpinic_flag='copy', readvar=readvar,  data=idata)
      if (readvar) then
        restart_file_spinup_state = idata
      else
        restart_file_spinup_state = spinup_state
      endif

    endif

  endif
  do jj = 1, this%num_rest_state1d
    ptr1d => this%rest_states_1d(:, jj)
    call restartvar(ncid=ncid, flag=flag, varname=trim(rest_varname_1d(jj)), &
       xtype=ncd_double,  dim1name='column', long_name='',  units='', &
       interpinic_flag='interp' , readvar=readvar, data=ptr1d)
  enddo
  do jj = 1, this%num_rest_state2d
    ptr2d => this%rest_states_2d(:, :, jj)
    call restartvar(ncid=ncid, flag=flag, varname=trim(rest_varname_2d(jj)), xtype=ncd_double,  &
      dim1name='column',dim2name='levtrc', switchdim=.true., &
      long_name='',  units='', interpinic_flag='interp',readvar=readvar, data=ptr2d)
  enddo
  if(trim(flag)=='read')then

    !assign initial conditions
    do c = bounds%begc, bounds%endc
      call this%betr(c)%set_restvar(betr_bounds, 1, betr_nlevtrc_soil, &
        this%num_rest_state1d,this%num_rest_state2d, &
        this%rest_states_1d(c:c,:), this%rest_states_2d(c:c,:,:), flag)
    enddo

    if(this%do_soibgc())then
      exit_spinup = (spinup_state == 0 .and. restart_file_spinup_state == 1 )
      enter_spinup = (spinup_state == 1 .and. restart_file_spinup_state == 0)
      if(get_nstep() >= 2)then
        exit_spinup = .false.; enter_spinup=.false.
      endif
      do c = bounds%begc, bounds%endc
        this%biophys_forc(c)%scalaravg_col(c_l) = max(this%scalaravg_col(c),0.01_r8)
        this%biophys_forc(c)%dom_scalar_col(c_l)= this%dom_scalar_col(c)
      enddo
    endif
  endif

  deallocate(rest_varname_1d)
  deallocate(rest_varname_2d)
  end subroutine BeTRSimulationRestart

  !------------------------------------------------------------------------

  subroutine BeTRSimulationSetSpinup(this, bounds)
  !
  ! set spinup for betr bgc runs
  use betr_ctrl      , only : exit_spinup, enter_spinup,betr_spinup_state
  use ApplicationsFactory, only : AppSetSpinup
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type), intent(in) :: bounds
  type(betr_bounds_type)     :: betr_bounds
  integer :: c

  if(exit_spinup .or. enter_spinup)then
     call AppSetSpinup()
     call this%BeTRSetBounds(betr_bounds)
     do c = bounds%begc, bounds%endc
       if(.not. this%active_col(c))cycle
       call this%betr(c)%set_bgc_spinup(betr_bounds, 1,  betr_nlevtrc_soil, this%biophys_forc(c))
     enddo
  endif
  if(exit_spinup)betr_spinup_state=0

  end subroutine BeTRSimulationSetSpinup

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetcps(this, bounds, col, pft)
  !
  !DESCRIPTION
  ! set up columns
  !USES
  use decompMod             , only : bounds_type
  use pftvarcon             , only : noveg, crop
  use tracer_varcon         , only : betr_nlevsoi
  !ARGUMENTS
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  type(bounds_type), intent(in) :: bounds
  type(column_type), intent(in) :: col
  type(patch_type), optional, intent(in) :: pft
  integer :: c, p, pi, pp
  integer :: c_l

  c_l=1
  do c = bounds%begc, bounds%endc
    this%betr_col(c)%snl(c_l) = col%snl(c)
    if(col%snl(c)<0)this%betr_col(c)%dz_snow(c_l, col%snl(c)+1:0) = col%dz(c,col%snl(c)+1:0)
    this%betr_col(c)%zi(c_l,0:betr_nlevsoi)= col%zi(c,0:betr_nlevsoi)
    this%betr_col(c)%dz(c_l,1:betr_nlevsoi)= col%dz(c,1:betr_nlevsoi)
    this%betr_col(c)%z(c_l,1:betr_nlevsoi)= col%z(c,1:betr_nlevsoi)

    this%biophys_forc(c)%zi(c_l,0:betr_nlevsoi)= col%zi(c,0:betr_nlevsoi)
    this%biophys_forc(c)%dz(c_l,1:betr_nlevsoi)= col%dz(c,1:betr_nlevsoi)
    this%biophys_forc(c)%z(c_l,1:betr_nlevsoi)= col%z(c,1:betr_nlevsoi)

    this%betr_col(c)%pfti(1)= col%pfti(c)
    this%betr_col(c)%pftf(1)= col%pftf(c)
    this%betr_col(c)%npfts(1)= col%npfts(c)

    if(present(pft))then
      this%betr_pft(c)%column(:)=1
      this%betr_pft(c)%npfts = 0
      pp = 0
      do pi = 1, betr_maxpatch_pft
        if (pi <= col%npfts(c)) then
          p = col%pfti(c) + pi - 1
          if (pft%active(p) .and. (pft%itype(p) /= noveg)) then
            pp = pp + 1
            this%betr_pft(c)%wtcol(pp) = pft%wtcol(p)
            this%betr_pft(c)%itype(pp) = pft%itype(p)
            this%betr_pft(c)%crop(pp) = crop(pi)         !the crop looks weird here, jyt
          endif
        endif
      enddo
      this%betr_pft(c)%npfts = pp
    endif
  enddo


  end subroutine BeTRSimulationSetcps

  !------------------------------------------------------------------------
  subroutine BeTRSimulationSetBounds(this, betr_bounds)
  !
  !DESCRIPTION
  !set betr_bounds
  !
  use betr_varcon    , only : betr_maxpatch_pft
  use tracer_varcon  , only : betr_nlevsoi,betr_nlevsno
  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  type(betr_bounds_type), intent(out)  :: betr_bounds
  !the following will be adpated to simulate lake and wetlands, by allowing
  !lbj to be negative
  betr_bounds%lbj  = 1; betr_bounds%ubj  = betr_nlevsoi
  betr_bounds%begp = 1; betr_bounds%endp = betr_maxpatch_pft
  betr_bounds%begc = 1; betr_bounds%endc = 1
  betr_bounds%begl = 1; betr_bounds%endl = 1
  betr_bounds%begg = 1; betr_bounds%endg = 1
  betr_bounds%nlevsno=betr_nlevsno
  end subroutine BeTRSimulationSetBounds

  !------------------------------------------------------------------------
  function do_soibgc(this)result(yesno)

  implicit none
  class(betr_simulation_type) , intent(inout) :: this

  logical :: yesno
  yesno = this%active_soibgc
  return
  end function do_soibgc
  !------------------------------------------------------------------------
  function do_regress_test(this) result(yesno)

  implicit none
  class(betr_simulation_type) , intent(inout) :: this
  logical :: yesno

  yesno = this%regression%write_regression_output
  end function do_regress_test


end module BeTRSimulation
