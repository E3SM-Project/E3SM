module BeTRSimulationALM
  !
  ! !DESCRIPTION:
  !  API for using BeTR in ALM
  !
  ! !USES:
  !
#include "shr_assert.h"
  use abortutils          , only : endrun
  use elm_varctl          , only : iulog,use_cn
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use BeTRSimulation      , only : betr_simulation_type
  use decompMod           , only : bounds_type
  use BeTRSimulation      , only : betr_simulation_type
  use BeTR_TimeMod        , only : betr_time_type
  use tracer_varcon       , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use betr_decompMod      , only : betr_bounds_type
  use betr_varcon         , only : betr_maxpatch_pft
  use pftvarcon           , only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
#if (defined SBETR)
  use PatchType      , only : patch_type
  use ColumnType     , only : column_type
  use LandunitType   , only : landunit_type
  use CNCarbonStateType  , only : carbonstate_type
  use CNNitrogenStateType, only : nitrogenstate_type
  use pfNitrogenStateType, only : pf_nitrogenstate_type
  use PhosphorusStateType, only : phosphorusstate_type
  use CNCarbonFluxType   , only : carbonflux_type
  use PfCarbonFluxType  , only : pf_carbonflux_type
  use CNNitrogenFluxType , only : nitrogenflux_type
  use PfNitrogenFluxType , only : pf_nitrogenflux_type
  use PhosphorusFluxType , only : phosphorusflux_type
  use PfPhosphorusFluxType , only : pf_phosphorusflux_type
  use WaterStateType  , only : waterstate_type
  use WaterfluxType     , only : waterflux_type
  use TemperatureType   , only : temperature_type
  use PfTemperatureType   , only : pf_temperature_type
  use PfWaterfluxType     , only : pf_waterflux_type
#else
  use VegetationType      , only : patch_type => vegetation_physical_properties
  use ColumnType          , only : column_type => column_physical_properties
  use LandunitType        , only : landunit_type => landunit_physical_properties
  use ColumnDataType      , only : carbonflux_type => column_carbon_flux
  use ColumnDataType      , only : nitrogenflux_type => column_nitrogen_flux
  use ColumnDataType      , only : phosphorusflux_type => column_phosphorus_flux
  use ColumnDataType      , only : carbonstate_type=> column_carbon_state
  use ColumnDataType      , only : nitrogenstate_type => column_nitrogen_state
  use ColumnDataType      , only : phosphorusstate_type => column_phosphorus_state
  use ColumnDataType      , only : waterstate_type => column_water_state
  use ColumnDataType      , only : waterflux_type => column_water_flux
  use ColumnDataType      , only : temperature_type=> column_energy_state
  use VegetationDataType  , only : vegetation_carbon_state, pf_carbonflux_type => vegetation_carbon_flux
  use VegetationDataType  , only : pf_nitrogenstate_type=>vegetation_nitrogen_state, pf_nitrogenflux_type=>vegetation_nitrogen_flux
  use VegetationDataType  , only : vegetation_phosphorus_state, pf_phosphorusflux_type=> vegetation_phosphorus_flux
  use VegetationDataType, only : pf_temperature_type => vegetation_energy_state
  use VegetationDataType, only : pf_waterflux_type => vegetation_water_flux
#endif
  use calibrationType, only : calibration_type
  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_alm_type
   private
     type(calibration_type) :: calibration_vpars
   contains
     procedure :: InitOnline                        => ALMInit
     procedure :: Init                              => ALMInitOffline
     procedure, public :: StepWithoutDrainage       => ALMStepWithoutDrainage
     procedure, public :: StepWithDrainage          => ALMStepWithDrainage
     procedure, public :: SetBiophysForcing         => ALMSetBiophysForcing
     !unique subroutines
     procedure, public :: CalcDewSubFlux            => ALMCalcDewSubFlux
     procedure, public :: CalcSmpL                  => ALMCalcSmpL
     procedure, public :: PlantSoilBGCSend          => ALMBetrPlantSoilBGCSend
     procedure, public :: RetriveBGCInput           => ALMBeTRRetriveBGCInput
     procedure, public :: PlantSoilBGCRecv          => ALMBetrPlantSoilBGCRecv
     procedure, public :: DiagnoseLnd2atm           => ALMDiagnoseLnd2atm
     procedure, public :: set_active                => ALMset_active
     procedure, public :: OutLoopSoilBGC            => ALMOutLoopSoilBGC
     procedure, public :: EnterOutLoopBGC           => ALMEnterOutLoopBGC
     procedure, public :: ExitOutLoopBGC            => ALMExitOutLoopBGC
     procedure, private:: set_transient_kinetics_par
     procedure, private:: set_vegpar_calibration
     procedure, public :: set_iP_prof
     procedure, public :: skip_balcheck
     procedure, public :: checkpmassyes
     procedure, public :: do_bgc_type
  end type betr_simulation_alm_type

  public :: create_betr_simulation_alm

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_alm() result(simulation)
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_simulation_alm_type), pointer :: simulation

    allocate(simulation)

  end function create_betr_simulation_alm
!-------------------------------------------------------------------------------

  function checkpmassyes(this)result(yesno)
  use tracer_varcon, only : fix_ip
  implicit none
  class(betr_simulation_alm_type)          , intent(inout) :: this

  logical :: yesno

  yesno = .not. fix_ip

  end function checkpmassyes
!-------------------------------------------------------------------------------
  subroutine Set_iP_prof(this, bounds)
  !
  !set initial inorganic P profile
  use decompMod       , only : bounds_type
  use elm_varpar      , only : nlevtrc_soil
  implicit none
  class(betr_simulation_alm_type)          , intent(inout) :: this
  type(bounds_type), intent(in) :: bounds

  integer :: c
  type(betr_bounds_type)   :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)

  betr_nlevtrc_soil = nlevtrc_soil

  do c = bounds%begc, bounds%endc
    if(.not. this%active_col(c))cycle
    call this%betr(c)%Set_iP_prof(betr_bounds, 1,  betr_nlevtrc_soil, this%biophys_forc(c))
  enddo

  end subroutine Set_iP_prof

!-------------------------------------------------------------------------------
  function skip_balcheck(this)result(stats)
  use betr_ctrl         , only : betr_spinup_state
  implicit none
  class(betr_simulation_alm_type)          , intent(inout) :: this

  logical :: stats

  stats = (this%spinup_count < 3) .and. (betr_spinup_state==3)
  return
  end function skip_balcheck
!-------------------------------------------------------------------------------

  subroutine ALMInit(this, bounds, lun, col, pft, waterstate, namelist_buffer, masterproc)
    !DESCRIPTION
    !Initialize BeTR for ALM
    !
    !USES
    !data types from alm
    use landunit_varcon , only : istcrop, istice, istsoil
    use elm_varpar      , only : nlevsno, nlevsoi, nlevtrc_soil
    !betr types
    use betr_constants      , only : betr_filename_length
    use betr_constants      , only : betr_namelist_buffer_size
    use BeTR_pftvarconType  , only : betr_pftvarcon
    use BeTR_landvarconType , only : betr_landvarcon
    use BeTR_decompMod      , only : betr_bounds_type
    use tracer_varcon       , only : do_bgc_calibration
    implicit none
    class(betr_simulation_alm_type)          , intent(inout) :: this
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(inout) :: col
    type(patch_type)                         , intent(in) :: pft
    character(len=*)                         , intent(in)    :: namelist_buffer
    logical,                      optional   , intent(in) :: masterproc
    type(waterstate_type)                    , intent(inout) :: waterstate

    !grid size
    betr_nlevsoi      = nlevsoi
    betr_nlevsno      = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    betr_landvarcon%istsoil            = istsoil
    betr_landvarcon%istcrop            = istcrop
    betr_landvarcon%istice             = istice
    this%spinup_count = 0
    ! now call the base simulation init to continue initialization
    if(present(masterproc))then
      call this%BeTRInit(bounds, lun, col, pft, waterstate, namelist_buffer, masterproc=masterproc)
    else
      call this%BeTRInit(bounds, lun, col, pft, waterstate, namelist_buffer)
    endif

    if(do_bgc_calibration)call this%calibration_vpars%Init(bounds, betr_maxpatch_pft)
  end subroutine ALMInit
!-------------------------------------------------------------------------------

  subroutine ALMInitOffline(this, bounds, lun, col, pft, waterstate, namelist_buffer, base_filename, case_id)
    !DESCRIPTION
    !Initialize BeTR for ALM
    !
    !USES
    !data types from alm
    use landunit_varcon , only : istcrop, istice, istsoil
    use elm_varpar      , only : nlevsno, nlevsoi, nlevtrc_soil
    !betr types
    use betr_constants      , only : betr_filename_length
    use betr_constants      , only : betr_namelist_buffer_size
    use BeTR_pftvarconType  , only : betr_pftvarcon
    use BeTR_landvarconType , only : betr_landvarcon
    use BeTR_decompMod      , only : betr_bounds_type

    implicit none
    class(betr_simulation_alm_type)          , intent(inout) :: this
    character(len=*)                         , intent(in)    :: namelist_buffer
    character(len=*)                         , intent(in)    :: base_filename
    character(len=*)                         , intent(in)    :: case_id
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(inout) :: col
    type(patch_type)                         , intent(in) :: pft
    type(waterstate_type)                    , intent(inout) :: waterstate

    !grid size
    betr_nlevsoi      = nlevsoi
    betr_nlevsno      = nlevsno
    betr_nlevtrc_soil = nlevtrc_soil

    betr_pftvarcon%nc3_arctic_grass    = nc3_arctic_grass
    betr_pftvarcon%nc3_nonarctic_grass = nc3_nonarctic_grass
    betr_pftvarcon%nc4_grass           = nc4_grass
    betr_pftvarcon%noveg               = noveg

    betr_landvarcon%istsoil            = istsoil
    betr_landvarcon%istcrop            = istcrop
    betr_landvarcon%istice             = istice

    ! now call the base simulation init to continue initialization
    call this%BeTRInit(bounds, lun, col, pft, waterstate, namelist_buffer, base_filename)

  end subroutine ALMInitOffline
!-------------------------------------------------------------------------------
  subroutine ALMStepWithoutDrainage(this, bounds,  col, pft)
   !DESCRIPTION
   !march one time step without doing drainage
   !
   !USES
    use elm_varpar        , only : nlevsno, nlevsoi, nlevtrc_soil
    use elm_varctl        , only : spinup_state
    use tracer_varcon     , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil, AA_spinup_on
    use betr_ctrl         , only : betr_spinup_state, enter_spinup
    use MathfuncMod       , only : num2str
    use betr_varcon       , only : kyr_spinup
    use clm_time_manager  , only : get_curr_date,is_end_curr_day,is_beg_curr_day,get_nstep
    implicit none
    ! !ARGUMENTS :
    class(betr_simulation_alm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds ! bounds
    type(column_type)               , intent(in)    :: col ! column type
    type(patch_type)                , intent(in)    :: pft
    !TEMPORARY VARIABLES
    type(betr_bounds_type)     :: betr_bounds
    integer :: c, c_l, begc_l, endc_l
    integer :: year, mon, day, sec

    call get_curr_date(year, mon, day, sec)
    c_l=1
    if(this%do_soibgc())then
      if(spinup_state==1)then
        do c = bounds%begc, bounds%endc
          this%biophys_forc(c)%dom_scalar_col(c_l)=this%dom_scalar_col(c)
        enddo
      else
        betr_spinup_state=0
      endif
    endif
    call this%bsimstatus%reset()

    !pass necessary data for correct subroutine call
    betr_nlevsoi       = nlevsoi
    betr_nlevsno       = nlevsno
    betr_nlevtrc_soil  = nlevtrc_soil

    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col, pft)

    c_l = 1; begc_l = betr_bounds%begc; endc_l=betr_bounds%endc;

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      !this%betr(c)%tracers%debug=col%debug_flag(c)
      call this%biophys_forc(c)%frac_normalize(this%betr_pft(c)%npfts, 1, betr_nlevtrc_soil)

      if(this%betr(c)%tracers%debug)call this%betr(c)%debug_info(betr_bounds, this%betr_col(c), &
         this%num_surfc, this%filter_soilc, 'bef w/o drain',this%bstatus(c))

      call this%betr(c)%step_without_drainage(this%betr_time, betr_bounds, this%betr_col(c), &
         this%betr_pft(c), this%num_surfc, this%filter_soilc, this%num_soilp, this%filter_soilp, &
         this%biophys_forc(c), this%biogeo_flux(c), this%biogeo_state(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err(),c)
        exit
      endif

      if(this%betr(c)%tracers%debug)call this%betr(c)%debug_info(betr_bounds, this%betr_col(c),&
          this%num_surfc, this%filter_soilc, 'aft w/o drain',this%bstatus(c))
    enddo
!    print*,'out without drainage'
    if(this%bsimstatus%check_status())then
      call endrun(msg=this%bsimstatus%print_msg())
    endif
    if(betr_spinup_state>0)then
      !the following needs double check for whether to keep or remove it.
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle
        this%dom_scalar_col(c) = this%biophys_forc(c)%dom_scalar_col(c_l)
      enddo
    endif

  end subroutine ALMStepWithoutDrainage

  !---------------------------------------------------------------------------------
  subroutine ALMset_active(this,bounds,col)

  !
  !DESCRIPTION
  !activate columuns that are active in alm
  implicit none
  ! !ARGUMENTS:
  class(betr_simulation_alm_type) , intent(inout) :: this
  type(bounds_type)               , intent(in)    :: bounds
  type(column_type)               , intent(in)    :: col ! column type

  integer :: c
  do c = bounds%begc, bounds%endc
    this%active_col(c) = (this%active_col(c) .and. col%active(c))
  enddo
  end subroutine ALMset_active
  !---------------------------------------------------------------------------------
  subroutine ALMBeTRRetriveBGCInput(this, num_surfc, filter_soilc, carbonflux_vars, nitrogenflux_vars, phosphorusflux_vars)


  implicit none
  class(betr_simulation_alm_type) , intent(inout) :: this
  integer           , intent(in)  :: num_surfc
  integer           , intent(in)  :: filter_soilc(:)
  type(carbonflux_type), intent(inout):: carbonflux_vars
  type(nitrogenflux_type), intent(inout):: nitrogenflux_vars
  type(phosphorusflux_type), intent(inout):: phosphorusflux_vars

  integer :: j, c, fc, c_l
  type(betr_bounds_type) :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)

  c_l=1

!  do j = betr_bounds%lbj, betr_bounds%ubj
!    do fc = 1, num_soilc
!      c = filter_soilc(fc)

!      carbonflux_vars%cflx_input_litr_met_vr_col(c,j)=this%biophys_forc(c)%c12flx%cflx_input_litr_met_vr_col(c_l,j)

!      carbonflux_vars%cflx_input_litr_cel_vr_col(c,j)=this%biophys_forc(c)%c12flx%cflx_input_litr_cel_vr_col(c_l,j)

!      carbonflux_vars%cflx_input_litr_lig_vr_col(c,j)=this%biophys_forc(c)%c12flx%cflx_input_litr_lig_vr_col(c_l,j)

!      carbonflux_vars%cflx_input_litr_cwd_vr_col(c,j)=this%biophys_forc(c)%c12flx%cflx_input_litr_cwd_vr_col(c_l,j)

!      nitrogenflux_vars%nflx_input_litr_met_vr_col(c,j)=this%biophys_forc(c)%n14flx%nflx_input_litr_met_vr_col(c_l,j)

!      nitrogenflux_vars%nflx_input_litr_cel_vr_col(c,j)=this%biophys_forc(c)%n14flx%nflx_input_litr_cel_vr_col(c_l,j)

!      nitrogenflux_vars%nflx_input_litr_lig_vr_col(c,j)=this%biophys_forc(c)%n14flx%nflx_input_litr_lig_vr_col(c_l,j)

!      nitrogenflux_vars%nflx_input_litr_cwd_vr_col(c,j)=this%biophys_forc(c)%n14flx%nflx_input_litr_cwd_vr_col(c_l,j)

!      phosphorusflux_vars%pflx_input_litr_met_vr_col(c,j)=this%biophys_forc(c)%p31flx%pflx_input_litr_met_vr_col(c_l,j)

!      phosphorusflux_vars%pflx_input_litr_cel_vr_col(c,j)=this%biophys_forc(c)%p31flx%pflx_input_litr_cel_vr_col(c_l,j)

!      phosphorusflux_vars%pflx_input_litr_lig_vr_col(c,j)=this%biophys_forc(c)%p31flx%pflx_input_litr_lig_vr_col(c_l,j)

!      phosphorusflux_vars%pflx_input_litr_cwd_vr_col(c,j)=this%biophys_forc(c)%p31flx%pflx_input_litr_cwd_vr_col(c_l,j)

!      nitrogenflux_vars%nflx_minn_input_nh4_vr_col(c,j)=this%biophys_forc(c)%n14flx%nflx_minn_input_nh4_vr_col(c_l,j) + &
!        this%biophys_forc(c)%n14flx%nflx_minn_nh4_fix_nomic_vr_col(c_l,j)

!      nitrogenflux_vars%nflx_minn_input_no3_vr_col(c,j)=this%biophys_forc(c)%n14flx%nflx_minn_input_no3_vr_col(c_l,j)

!      phosphorusflux_vars%pflx_minp_input_po4_vr_col(c,j)=this%biophys_forc(c)%p31flx%pflx_minp_input_po4_vr_col(c_l,j) + &
!        this%biophys_forc(c)%p31flx%pflx_minp_weathering_po4_vr_col(c_l,j)
!    enddo
!  enddo
  end subroutine ALMBeTRRetriveBGCInput
  !---------------------------------------------------------------------------------
  subroutine ALMStepWithDrainage(this, bounds,  col)
   !
   !DESCRIPTION
   !interface for using diagnose land fluxes to atm and river copmonents
   !
   !USES
    use MathfuncMod   , only : safe_div
    use lnd2atmType    , only : lnd2atm_type
    use elm_varpar     , only : nlevsno, nlevsoi, nlevtrc_soil
    use tracer_varcon  , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    implicit none
    !ARGUMENTS
    class(betr_simulation_alm_type) , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    type(column_type)           , intent(in)    :: col ! column type

    !temporary variables
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer                :: c, c_l, begc_l, endc_l

    call this%bsimstatus%reset()

    betr_nlevsoi       = nlevsoi
    betr_nlevsno       = nlevsno
    betr_nlevtrc_soil  = nlevtrc_soil

    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col)

    c_l = 1; begc_l = betr_bounds%begc; endc_l=betr_bounds%endc;

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle

      call this%betr(c)%step_with_drainage(betr_bounds,      &
         this%betr_col(c),this%num_surfc, this%filter_soilc, this%jtops, &
         this%biogeo_flux(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif

    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=this%bsimstatus%print_msg())

  end subroutine ALMStepWithDrainage

  !---------------------------------------------------------------------------------
  subroutine ALMDiagnoseLnd2atm(this, bounds,  col, lnd2atm_vars)
   !DESCRIPTION
   ! march one step with drainage
   !
   !USES
    use subgridAveMod  , only : c2g
    use elm_varpar     , only : nlevsno, nlevsoi, nlevtrc_soil
    use lnd2atmType    , only : lnd2atm_type
    use betr_decompMod , only : betr_bounds_type
    use tracer_varcon  , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    use tracer_varcon  , only : reaction_method
    implicit none
    ! !ARGUMENTS:
    class(betr_simulation_alm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    type(column_type)               , intent(in)    :: col ! column type
    type(lnd2atm_type)              , intent(inout) :: lnd2atm_vars

    !temporary variables
    type(betr_bounds_type) :: betr_bounds
    integer                :: lbj, ubj ! lower and upper bounds, make sure they are > 0
    integer                :: c, c_l
    real(r8)  :: qflx_rofliq_qsur_doc_col(bounds%begc:bounds%endc)
    real(r8)  :: qflx_rofliq_qsur_dic_col(bounds%begc:bounds%endc)
    real(r8)  :: qflx_rofliq_qsub_doc_col(bounds%begc:bounds%endc)
    real(r8)  :: qflx_rofliq_qsub_dic_col(bounds%begc:bounds%endc)

    associate(  &
     begc => bounds%begc, &
     endc => bounds%endc, &
     begg => bounds%begg, &
     endg => bounds%endg, &
     qflx_rofliq_qsur_doc_grc  => lnd2atm_vars%qflx_rofliq_qsur_doc_grc, &
     qflx_rofliq_qsur_dic_grc  => lnd2atm_vars%qflx_rofliq_qsur_dic_grc, &
     qflx_rofliq_qsub_doc_grc  => lnd2atm_vars%qflx_rofliq_qsub_doc_grc, &
     qflx_rofliq_qsub_dic_grc  => lnd2atm_vars%qflx_rofliq_qsub_dic_grc  &
    )

    call this%BeTRSetBounds(betr_bounds)

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
         call this%betr(c)%diagnoselnd2atm(betr_bounds,      &
           this%num_surfc, this%filter_soilc, this%biogeo_flux(c))
    enddo

    if(index(reaction_method,'mosart')/=0)then

      c_l = 1
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle
        qflx_rofliq_qsur_doc_col(c)=this%biogeo_flux(c)%c12flux_vars%som_c_runoff_col(c_l)
        qflx_rofliq_qsur_dic_col(c)=this%biogeo_flux(c)%qflx_rofliq_qsur_dic_col(c_l)
        qflx_rofliq_qsub_doc_col(c)=this%biogeo_flux(c)%c12flux_vars%som_c_qdrain_col(c_l)
        qflx_rofliq_qsub_dic_col(c)=this%biogeo_flux(c)%qflx_rofliq_qsub_dic_col(c_l)
      enddo

      call c2g( bounds, &
           qflx_rofliq_qsur_doc_col(begc:endc), qflx_rofliq_qsur_doc_grc(begg:endg),     &
           c2l_scale_type= 'unity', l2g_scale_type='unity' )

      call c2g( bounds, &
           qflx_rofliq_qsur_dic_col(begc:endc), qflx_rofliq_qsur_dic_grc(begg:endg),     &
           c2l_scale_type= 'unity', l2g_scale_type='unity' )

      call c2g( bounds, &
           qflx_rofliq_qsub_doc_grc(begc:endc), qflx_rofliq_qsub_doc_grc(begg:endg),     &
           c2l_scale_type= 'unity', l2g_scale_type='unity' )

      call c2g( bounds, &
           qflx_rofliq_qsub_doc_grc(begc:endc), qflx_rofliq_qsub_doc_grc(begg:endg),     &
           c2l_scale_type= 'unity', l2g_scale_type='unity' )
    endif
    end associate
  end subroutine ALMDiagnoseLnd2atm

  !------------------------------------------------------------------------
  subroutine ALMBetrPlantSoilBGCSend(this, bounds, col, pft, num_surfc,  filter_soilc, cnstate_vars, &
    carbonstate_vars, carbonflux_vars,  c13state_vars, c13_cflx_vars, &
    c14state_vars, c14_cflx_vars, nitrogenstate_vars, nitrogenflux_vars, &
    phosphorusstate_vars, phosphorusflux_vars, &
    PlantMicKinetics_vars)

  !read in biogeochemical fluxes from alm for soil bgc modeling
  !these are C, N and P fluxes from root input, surface litter input
  !atmospheric deposition, fire (negative), and fertilization
  !Because of possible harvest activity that is
  !related to dynamic land use, input profiles are computed in alm.
  !

  use CNStateType, only : cnstate_type
  use elm_varpar, only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use CNStateType        , only : cnstate_type
  use elm_varpar         , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use PlantMicKineticsMod, only : PlantMicKinetics_type
  use mathfuncMod        , only : apvb,bisnan
  use tracer_varcon      , only : use_c13_betr, use_c14_betr, do_bgc_calibration
  use tracer_varcon      , only : betr_nlevsoi
  implicit none
  class(betr_simulation_alm_type), intent(inout)  :: this
  type(bounds_type) , intent(in)  :: bounds
  type(column_type) , intent(in)  :: col ! column type
  type(patch_type)  , intent(in)  :: pft ! pft type
  integer           , intent(in)  :: num_surfc
  integer           , intent(in)  :: filter_soilc(:)
  type(cnstate_type), intent(in)  :: cnstate_vars
  type(carbonstate_type), intent(in) :: carbonstate_vars
  type(carbonflux_type), intent(in)  :: carbonflux_vars
  type(carbonstate_type), intent(in) :: c13state_vars
  type(carbonflux_type), intent(in) :: c13_cflx_vars
  type(carbonstate_type), intent(in) :: c14state_vars
  type(carbonflux_type), intent(in):: c14_cflx_vars
  type(nitrogenstate_type),intent(in):: nitrogenstate_vars
  type(nitrogenflux_type), intent(in):: nitrogenflux_vars
  type(phosphorusstate_type), intent(in):: phosphorusstate_vars
  type(phosphorusflux_type), intent(in):: phosphorusflux_vars
  type(PlantMicKinetics_type), intent(in) :: PlantMicKinetics_vars

  !temporary variables
  type(betr_bounds_type) :: betr_bounds
  integer :: c, fc, j, c_l, begc_l, endc_l, kk
  real(r8) :: fport
  real(r8) :: ndep_prof_loc(1:betr_nlevsoi)

  associate(                                           &
    ndep_prof     => cnstate_vars%ndep_prof_col     ,  &
    pdep_prof     => cnstate_vars%pdep_prof_col     ,  &
    nfixation_prof=> cnstate_vars%nfixation_prof_col,  &
    biochem_pmin_vr=> phosphorusflux_vars%biochem_pmin_vr, &
    frac_loss_lit_to_fire_col=> cnstate_vars%frac_loss_lit_to_fire_col, &
    frac_loss_cwd_to_fire_col=> cnstate_vars%frac_loss_cwd_to_fire_col  &
  )
  call this%BeTRSetBounds(betr_bounds)

  !reset and prepare for retrieval
  do fc = 1, num_surfc
    c = filter_soilc(fc)
    call this%biogeo_flux(c)%reset(value_column=0._r8, active_soibgc=this%do_soibgc())
    call this%biogeo_state(c)%reset(value_column=0._r8, active_soibgc=this%do_soibgc())
  enddo
  if(.not. this%do_soibgc())return
  if(this%do_bgc_type('type0_bgc'))return
  c_l=1
  if(this%do_bgc_type('type1_bgc'))then
    !transfer state variables from elm to betr
    associate(                                                             &
    c12_decomp_cpools_vr_col =>   carbonstate_vars%decomp_cpools_vr   , &
    decomp_npools_vr_col =>   nitrogenstate_vars%decomp_npools_vr    , &
    decomp_ppools_vr_col =>   phosphorusstate_vars%decomp_ppools_vr    , &
    smin_no3_vr         => nitrogenstate_vars%smin_no3_vr         , &
    smin_nh4_vr         => nitrogenstate_vars%smin_nh4_vr         , &
    solutionp_vr        => phosphorusstate_vars%solutionp_vr          &
    )
    do j = 1,betr_nlevtrc_soil
      do fc = 1, num_surfc
        c = filter_soilc(fc)
        this%biophys_forc(c)%n14flx%in_sminn_no3_vr_col(c_l,j) = smin_no3_vr(c,j)
        this%biophys_forc(c)%n14flx%in_sminn_nh4_vr_col(c_l,j) = smin_nh4_vr(c,j)
        this%biophys_forc(c)%p31flx%in_sminp_vr_col(c_l,j) = solutionp_vr(c,j)
      enddo
    enddo

    do kk = 1, 7
      do j = 1,betr_nlevtrc_soil
        do fc = 1, num_surfc
          c = filter_soilc(fc)
          this%biophys_forc(c)%c12flx%in_decomp_cpools_vr_col(c_l,j,kk)=c12_decomp_cpools_vr_col(c,j,kk)
          this%biophys_forc(c)%n14flx%in_decomp_npools_vr_col(c_l,j,kk)=decomp_npools_vr_col(c,j,kk)
          this%biophys_forc(c)%p31flx%in_decomp_ppools_vr_col(c_l,j,kk)=decomp_ppools_vr_col(c,j,kk)
        enddo
      enddo
    enddo

    end associate
    do fc = 1, num_surfc
      c = filter_soilc(fc)
      call this%betr(c)%reset_biostates(betr_bounds, this%num_surfc, this%filter_soilc, &
        this%biophys_forc(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err(),c)
        exit
      endif
    enddo
    return
  endif
  !Note for improvement:
  !The following is a work around of the nutrient deposition problem.
  !Due to the lack of information on dry and wet deposition partitioning,
  !all deposition is distributed in the top 2 soil layers. This treatment
  !is overestimating incoming nutrient flux to soil in cold regions.
  !Nov, 2018
  ndep_prof_loc(:) = 0._r8
  ndep_prof_loc(1) = col%dz(bounds%begc,1)/(col%dz(bounds%begc,1)+col%dz(bounds%begc,2))
  ndep_prof_loc(2) = 1._r8-ndep_prof_loc(1)
  ndep_prof_loc(1) = ndep_prof_loc(1)/col%dz(bounds%begc,1)
  ndep_prof_loc(2) = ndep_prof_loc(2)/col%dz(bounds%begc,2)

  !set kinetic parameters
  call this%set_transient_kinetics_par(betr_bounds, col, pft, num_surfc, filter_soilc, PlantMicKinetics_vars)

  if(do_bgc_calibration)then
    call this%set_vegpar_calibration(betr_bounds, col, pft, num_surfc, filter_soilc, this%calibration_vpars)
  endif

  !set biophysical forcing
  c_l = 1
  do fc = 1, num_surfc
    c = filter_soilc(fc)
    call this%biophys_forc(c)%reset(value_column=0._r8)
    this%biophys_forc(c)%annsum_counter_col(c_l) = cnstate_vars%annsum_counter_col(c)
    this%biophys_forc(c)%isoilorder(c_l) = 1                 !this needs update
!    this%biophys_forc(c)%lithotype_col(c_l) = cnstate_vars%lithoclass_col(c)
    this%biophys_forc(c)%frac_loss_lit_to_fire_col(c_l) = frac_loss_lit_to_fire_col(c)
    this%biophys_forc(c)%frac_loss_cwd_to_fire_col(c_l) = frac_loss_cwd_to_fire_col(c)
    this%biophys_forc(c)%biochem_pmin_vr(c_l,1:betr_nlevsoi)= biochem_pmin_vr(c,1:betr_nlevsoi)
    call this%biophys_forc(c)%c12flx%reset(value_column=0._r8)
    call this%biophys_forc(c)%n14flx%reset(value_column=0._r8)
    call this%biophys_forc(c)%p31flx%reset(value_column=0._r8)

    if(use_c13_betr)then
      call this%biophys_forc(c)%c13flx%reset(value_column=0._r8)
    endif

    if(use_c14_betr)then
      call this%biophys_forc(c)%c14flx%reset(value_column=0._r8)
    endif
  enddo

  !sum up carbon input profiles
  do j = 1, betr_bounds%ubj
    do fc = 1, num_surfc
      c = filter_soilc(fc)
      this%biophys_forc(c)%pweath_prof_col(c_l,j) = pdep_prof(c,j)
      this%biophys_forc(c)%c12flx%rt_vr_col(c_l,j) = carbonflux_vars%rr_vr(c,j)

      !!------------------------------------------------------------------------
      !carbon input
      !metabolic carbon
      call apvb(this%biophys_forc(c)%c12flx%cflx_input_litr_met_vr_col(c_l,j), & !
         (/carbonflux_vars%phenology_c_to_litr_met_c(c,j)       , & !phenology
         carbonflux_vars%dwt_frootc_to_litr_met_c(c,j)          , & !dynamic land cover
         carbonflux_vars%gap_mortality_c_to_litr_met_c(c,j)     , & !gap mortality
         carbonflux_vars%harvest_c_to_litr_met_c(c,j)           , & !harvest
         carbonflux_vars%m_c_to_litr_met_fire(c,j)/))              ! fire mortality

      !cellulose carbon
      call apvb(this%biophys_forc(c)%c12flx%cflx_input_litr_cel_vr_col(c_l,j)   , &
         (/carbonflux_vars%phenology_c_to_litr_cel_c(c,j)          , &  !phenology
         carbonflux_vars%dwt_frootc_to_litr_cel_c(c,j)             , &  !dynamic land cover
         carbonflux_vars%gap_mortality_c_to_litr_cel_c(c,j)        , & !gap mortality
         carbonflux_vars%harvest_c_to_litr_cel_c(c,j)              , & !harvest
         carbonflux_vars%m_c_to_litr_cel_fire(c,j)/))              ! fire mortality

      !lignin carbon
      call apvb(this%biophys_forc(c)%c12flx%cflx_input_litr_lig_vr_col(c_l,j) , &
         (/carbonflux_vars%phenology_c_to_litr_lig_c(c,j)        , & !phenology
         carbonflux_vars%dwt_frootc_to_litr_lig_c(c,j)           , & !dynamic land cover
         carbonflux_vars%gap_mortality_c_to_litr_lig_c(c,j)      , & !gap mortality
         carbonflux_vars%harvest_c_to_litr_lig_c(c,j)            , & !harvest
         carbonflux_vars%m_c_to_litr_lig_fire(c,j)/))                ! fire mortality

      !cwd carbon
      call apvb(this%biophys_forc(c)%c12flx%cflx_input_litr_cwd_vr_col(c_l,j) , &
        (/carbonflux_vars%dwt_livecrootc_to_cwdc(c,j)            , &
        carbonflux_vars%dwt_deadcrootc_to_cwdc(c,j)              , &
        carbonflux_vars%gap_mortality_c_to_cwdc(c,j)             , &
        carbonflux_vars%harvest_c_to_cwdc(c,j)                   , &
        carbonflux_vars%fire_mortality_c_to_cwdc(c,j)/))


      !!------------------------------------------------------------------------
      if(use_c13_betr)then
        !metabolic carbon
        call apvb(this%biophys_forc(c)%c13flx%cflx_input_litr_met_vr_col(c_l,j), & !
           (/carbonflux_vars%phenology_c_to_litr_met_c(c,j)       , & !phenology
           carbonflux_vars%dwt_frootc_to_litr_met_c(c,j)          , & !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_met_c(c,j)     , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_met_c(c,j)           , & !harvest
           carbonflux_vars%m_c_to_litr_met_fire(c,j)/))              ! fire mortality

        !cellulose carbon
        call apvb(this%biophys_forc(c)%c13flx%cflx_input_litr_cel_vr_col(c_l,j)   , &
           (/carbonflux_vars%phenology_c_to_litr_cel_c(c,j)          , &  !phenology
           carbonflux_vars%dwt_frootc_to_litr_cel_c(c,j)             , &  !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_cel_c(c,j)        , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_cel_c(c,j)              , & !harvest
           carbonflux_vars%m_c_to_litr_cel_fire(c,j)/))              ! fire mortality

        !lignin carbon
        call apvb(this%biophys_forc(c)%c13flx%cflx_input_litr_lig_vr_col(c_l,j) , &
           (/carbonflux_vars%phenology_c_to_litr_lig_c(c,j)        , & !phenology
           carbonflux_vars%dwt_frootc_to_litr_lig_c(c,j)           , & !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_lig_c(c,j)      , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_lig_c(c,j)            , & !harvest
           carbonflux_vars%m_c_to_litr_lig_fire(c,j)/))                ! fire mortality

        !cwd carbon
        call apvb(this%biophys_forc(c)%c13flx%cflx_input_litr_cwd_vr_col(c_l,j) , &
          (/carbonflux_vars%dwt_livecrootc_to_cwdc(c,j)            , &
          carbonflux_vars%dwt_deadcrootc_to_cwdc(c,j)              , &
          carbonflux_vars%gap_mortality_c_to_cwdc(c,j)             , &
          carbonflux_vars%harvest_c_to_cwdc(c,j)                   , &
          carbonflux_vars%fire_mortality_c_to_cwdc(c,j)/))

      endif
      if(use_c14_betr)then
        !metabolic carbon
        call apvb(this%biophys_forc(c)%c14flx%cflx_input_litr_met_vr_col(c_l,j), & !
           (/carbonflux_vars%phenology_c_to_litr_met_c(c,j)       , & !phenology
           carbonflux_vars%dwt_frootc_to_litr_met_c(c,j)          , & !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_met_c(c,j)     , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_met_c(c,j)           , & !harvest
           carbonflux_vars%m_c_to_litr_met_fire(c,j)/))              ! fire mortality

        !cellulose carbon
        call apvb(this%biophys_forc(c)%c14flx%cflx_input_litr_cel_vr_col(c_l,j)   , &
           (/carbonflux_vars%phenology_c_to_litr_cel_c(c,j)          , &  !phenology
           carbonflux_vars%dwt_frootc_to_litr_cel_c(c,j)             , &  !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_cel_c(c,j)        , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_cel_c(c,j)              , & !harvest
           carbonflux_vars%m_c_to_litr_cel_fire(c,j)/))              ! fire mortality

        !lignin carbon
        call apvb(this%biophys_forc(c)%c14flx%cflx_input_litr_lig_vr_col(c_l,j) , &
           (/carbonflux_vars%phenology_c_to_litr_lig_c(c,j)        , & !phenology
           carbonflux_vars%dwt_frootc_to_litr_lig_c(c,j)           , & !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_lig_c(c,j)      , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_lig_c(c,j)            , & !harvest
           carbonflux_vars%m_c_to_litr_lig_fire(c,j)/))                ! fire mortality

        !cwd carbon
        call apvb(this%biophys_forc(c)%c14flx%cflx_input_litr_cwd_vr_col(c_l,j) , &
          (/carbonflux_vars%dwt_livecrootc_to_cwdc(c,j)            , &
          carbonflux_vars%dwt_deadcrootc_to_cwdc(c,j)              , &
          carbonflux_vars%gap_mortality_c_to_cwdc(c,j)             , &
          carbonflux_vars%harvest_c_to_cwdc(c,j)                   , &
          carbonflux_vars%fire_mortality_c_to_cwdc(c,j)/))

      endif

      !nitrogen input
      !metabolic nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_input_litr_met_vr_col(c_l,j) , &
         (/nitrogenflux_vars%phenology_n_to_litr_met_n(c,j)      , & !phenology
         nitrogenflux_vars%dwt_frootn_to_litr_met_n(c,j)         , & !dynamic land cover
         nitrogenflux_vars%gap_mortality_n_to_litr_met_n(c,j)    , & !gap mortality
         nitrogenflux_vars%harvest_n_to_litr_met_n(c,j)          , & !harvest
         nitrogenflux_vars%m_n_to_litr_met_fire(c,j)/))              ! fire mortality

      !cellulose nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_input_litr_cel_vr_col(c_l,j), &
         (/nitrogenflux_vars%phenology_n_to_litr_cel_n(c,j)     , & !phenology
         nitrogenflux_vars%dwt_frootn_to_litr_cel_n(c,j)        , & !dynamic land cover
         nitrogenflux_vars%gap_mortality_n_to_litr_cel_n(c,j)   , & !gap mortality
         nitrogenflux_vars%harvest_n_to_litr_cel_n(c,j)         , & !harvest
         nitrogenflux_vars%m_n_to_litr_cel_fire(c,j)/))             ! fire mortality

      !lignin nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_input_litr_lig_vr_col(c_l,j) , &
         (/nitrogenflux_vars%phenology_n_to_litr_lig_n(c,j)      , &  !phenology
         nitrogenflux_vars%dwt_frootn_to_litr_lig_n(c,j)         , &   !dynamic land cover
         nitrogenflux_vars%gap_mortality_n_to_litr_lig_n(c,j)    , & !gap mortality
         nitrogenflux_vars%harvest_n_to_litr_lig_n(c,j)          , & !harvest
         nitrogenflux_vars%m_n_to_litr_lig_fire(c,j)/))              ! fire mortality

      !cwd nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_input_litr_cwd_vr_col(c_l,j) , &
        (/nitrogenflux_vars%dwt_livecrootn_to_cwdn(c,j)          , &
        nitrogenflux_vars%dwt_deadcrootn_to_cwdn(c,j)            , &
        nitrogenflux_vars%gap_mortality_n_to_cwdn(c,j)           , &
        nitrogenflux_vars%harvest_n_to_cwdn(c,j)                 , &
        nitrogenflux_vars%fire_mortality_n_to_cwdn(c,j)/))

      !phosphorus input
      !metabolic phosphorus
      call apvb(this%biophys_forc(c)%p31flx%pflx_input_litr_met_vr_col(c_l,j) , &
         (/phosphorusflux_vars%phenology_p_to_litr_met_p(c,j)    , & !phenology
         phosphorusflux_vars%dwt_frootp_to_litr_met_p(c,j)       , & !dynamic land cover
         phosphorusflux_vars%gap_mortality_p_to_litr_met_p(c,j)  , & !gap mortality
         phosphorusflux_vars%harvest_p_to_litr_met_p(c,j)        , & !harvest
         phosphorusflux_vars%m_p_to_litr_met_fire(c,j)/))            ! fire mortality

      !cellulose phosphorus
      call apvb(this%biophys_forc(c)%p31flx%pflx_input_litr_cel_vr_col(c_l,j) , &
         (/phosphorusflux_vars%phenology_p_to_litr_cel_p(c,j)    , & !phenology
         phosphorusflux_vars%dwt_frootp_to_litr_cel_p(c,j)       , & !dynamic land cover
         phosphorusflux_vars%gap_mortality_p_to_litr_cel_p(c,j)  , & !gap mortality
         phosphorusflux_vars%harvest_p_to_litr_cel_p(c,j)        , & !harvest
         phosphorusflux_vars%m_p_to_litr_cel_fire(c,j)/))            ! fire mortality

      !lignin phosphorus
      call apvb(this%biophys_forc(c)%p31flx%pflx_input_litr_lig_vr_col(c_l,j) , &
         (/phosphorusflux_vars%phenology_p_to_litr_lig_p(c,j)    , & !phenology
         phosphorusflux_vars%dwt_frootp_to_litr_lig_p(c,j)       , & !dynamic land cover
         phosphorusflux_vars%gap_mortality_p_to_litr_lig_p(c,j)  , & !gap mortality
         phosphorusflux_vars%harvest_p_to_litr_lig_p(c,j)        , & !harvest
         phosphorusflux_vars%m_p_to_litr_lig_fire(c,j)/))            ! fire mortality

      !cwd phosphorus
      call apvb(this%biophys_forc(c)%p31flx%pflx_input_litr_cwd_vr_col(c_l,j) , &
        (/phosphorusflux_vars%dwt_livecrootp_to_cwdp(c,j) , &
        phosphorusflux_vars%dwt_deadcrootp_to_cwdp(c,j)   , &
        phosphorusflux_vars%gap_mortality_p_to_cwdp(c,j)  , &
        phosphorusflux_vars%harvest_p_to_cwdp(c,j)        , &
        phosphorusflux_vars%fire_mortality_p_to_cwdp(c,j)/))

      !mineral nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_minn_input_nh4_vr_col(c_l,j) , &
         (/nitrogenflux_vars%ndep_to_sminn_nh3(c)                    , &
         nitrogenflux_vars%fert_to_sminn(c)/),  ndep_prof_loc(j))

      call apvb(this%biophys_forc(c)%n14flx%nflx_minn_input_no3_vr_col(c_l,j) , &
         (/nitrogenflux_vars%ndep_to_sminn_no3(c)/),  ndep_prof_loc(j))

      !the following could be commented out if a fixation model is done in betr
      call apvb(this%biophys_forc(c)%n14flx%nflx_minn_nh4_fix_nomic_vr_col(c_l,j) , &
         (/nitrogenflux_vars%nfix_to_sminn(c)                        , &
         nitrogenflux_vars%soyfixn_to_sminn(c)/), nfixation_prof(c,j))

      !mineral phosphorus, the deposition is assumed to be of primary form
      call apvb(this%biophys_forc(c)%p31flx%pflx_minp_input_po4_vr_col(c_l,j) , &
         (/phosphorusflux_vars%fert_p_to_sminp(c)/),   ndep_prof_loc(j))

    enddo
  enddo

  end associate
  !pull in all state variables and update tracers
  end subroutine ALMBetrPlantSoilBGCSend

  !------------------------------------------------------------------------

  subroutine ALMBetrPlantSoilBGCRecv(this, bounds, col, pft, num_surfc,  filter_soilc,&
   c12state_vars, c12flux_vars, pf_c12flux_vars, c13state_vars, c13flux_vars, &
   c14state_vars, c14flux_vars, n14state_vars, pf_n14state_vars, n14flux_vars, &
   pf_n14flux_vars, p31state_vars, p31flux_vars, pf_p31flux_vars)

  !this returns the flux back to ALM after doing soil BGC
  !this specifically returns plant nutrient yield
  use clm_time_manager    , only : get_nstep

  use tracer_varcon       , only : use_c13_betr, use_c14_betr
  use MathfuncMod         , only : safe_div
  use tracer_varcon       , only : reaction_method

  implicit none
  class(betr_simulation_alm_type), intent(inout)  :: this
  type(bounds_type) , intent(in)  :: bounds
  type(patch_type)            , intent(in) :: pft
  type(column_type)           , intent(in)    :: col ! column type
  integer           , intent(in)  :: num_surfc
  integer           , intent(in)  :: filter_soilc(:)
  type(carbonstate_type), intent(inout) :: c12state_vars
  type(carbonstate_type), intent(inout) :: c13state_vars
  type(carbonstate_type), intent(inout) :: c14state_vars
  type(nitrogenstate_type), intent(inout) :: n14state_vars
  type(phosphorusstate_type), intent(inout) :: p31state_vars
  type(carbonflux_type)  , intent(inout):: c12flux_vars    !return carbon fluxes through DON?
  type(carbonflux_type)  , intent(inout):: c13flux_vars    !return carbon fluxes through DON?
  type(carbonflux_type)  , intent(inout):: c14flux_vars    !return carbon fluxes through DON?
  type(nitrogenflux_type), intent(inout):: n14flux_vars
  type(phosphorusflux_type), intent(inout):: p31flux_vars
  type(pf_carbonflux_type), intent(inout) :: pf_c12flux_vars
  type(pf_nitrogenstate_type), intent(inout) :: pf_n14state_vars
  type(pf_nitrogenflux_type), intent(inout) :: pf_n14flux_vars
  type(pf_phosphorusflux_type), intent(inout) :: pf_p31flux_vars

  integer :: c, fc, p, pi, c_l, j, kk

    !TEMPORARY VARIABLES
  type(betr_bounds_type)     :: betr_bounds
  integer :: begc_l, endc_l

  !summarize the fluxes and state variables
  c_l = 1
  call this%BeTRSetBounds(betr_bounds)
  begc_l = betr_bounds%begc; endc_l=betr_bounds%endc;

  !retrieve and return
  do fc =1, num_surfc
    c = filter_soilc(fc)
    call this%betr(c)%retrieve_biostates(betr_bounds,  1, betr_nlevsoi, &
       this%num_surfc, this%filter_soilc,&
       this%jtops, this%biogeo_state(c),this%bstatus(c))

    if(this%bstatus(c)%check_status())then
      call this%bsimstatus%setcol(c)
      call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
      exit
    endif
    call this%betr(c)%retrieve_biofluxes(this%num_surfc, this%filter_soilc, &
      this%num_soilp, this%filter_soilp,  this%biogeo_flux(c))

    call this%biogeo_state(c)%summary(betr_bounds, 1, betr_nlevtrc_soil,&
         this%betr_col(c)%dz(begc_l:endc_l,1:betr_nlevtrc_soil), &
         this%betr_col(c)%zi(begc_l:endc_l,1:betr_nlevtrc_soil),this%do_soibgc())

    call this%biogeo_flux(c)%summary(betr_bounds, 1, betr_nlevtrc_soil, &
         this%betr_col(c)%dz(begc_l:endc_l,1:betr_nlevtrc_soil),this%do_soibgc())
  enddo

  if(.not. this%do_soibgc())return

  !retrieve plant nutrient uptake from biogeo_flux
  if (this%do_bgc_type('type2_bgc')) then
    do fc = 1, num_surfc
      c = filter_soilc(fc)
      pi = 0
      do p = col%pfti(c), col%pftf(c)
        if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
          pi = pi + 1
          pf_n14flux_vars%smin_nh4_to_plant(p) = this%biogeo_flux(c)%n14flux_vars%smin_nh4_to_plant_patch(pi)
          pf_n14flux_vars%smin_no3_to_plant(p) = this%biogeo_flux(c)%n14flux_vars%smin_no3_to_plant_patch(pi)
          pf_p31flux_vars%sminp_to_plant(p)  = this%biogeo_flux(c)%p31flux_vars%sminp_to_plant_patch(pi)
          pf_p31flux_vars%sminp_to_plant_trans(p) = this%biogeo_flux(c)%p31flux_vars%sminp_to_plant_trans_patch(pi)
          !compute relative n return, note the following computation is different from ALM-ECA-CNP, because
          !betr includes transpiration incuded nitrogen uptake, which has not direct temperature sensitivity.
          pf_n14state_vars%pnup_pfrootc(p) = this%biogeo_flux(c)%pnup_pfrootc_patch(pi)
          pf_c12flux_vars%tempavg_agnpp(p) = this%biogeo_flux(c)%c12flux_vars%tempavg_agnpp_patch(pi)
          pf_c12flux_vars%annavg_agnpp(p)  = this%biogeo_flux(c)%c12flux_vars%annavg_agnpp_patch(pi)
          pf_c12flux_vars%tempavg_bgnpp(p) = this%biogeo_flux(c)%c12flux_vars%tempavg_bgnpp_patch(pi)
          pf_c12flux_vars%annavg_bgnpp(p)  = this%biogeo_flux(c)%c12flux_vars%annavg_bgnpp_patch(pi)
        else
          pf_n14flux_vars%smin_nh4_to_plant(p) = 0._r8
          pf_n14flux_vars%smin_no3_to_plant(p) = 0._r8
          pf_p31flux_vars%sminp_to_plant(p) = 0._r8
        endif
      enddo

      !recollect soil respirations, fire and hydraulic loss
      c12flux_vars%hr(c) = this%biogeo_flux(c)%c12flux_vars%hr_col(c_l)

      c12flux_vars%fire_decomp_closs(c) = this%biogeo_flux(c)%c12flux_vars%fire_decomp_closs_col(c_l)
      c12flux_vars%som_c_leached(c) = &
        this%biogeo_flux(c)%c12flux_vars%som_c_leached_col(c_l) + &
        this%biogeo_flux(c)%c12flux_vars%som_c_qdrain_col(c_l)
      c12flux_vars%som_c_runoff(c) = this%biogeo_flux(c)%c12flux_vars%som_c_runoff_col(c_l)
        !the following is for consistency with the ALM definitation, which computes
        !som_c_leached_col as a numerical roundoff
      c12flux_vars%som_c_leached(c)=-c12flux_vars%som_c_leached(c)
      if(use_c13_betr)then
        c13flux_vars%hr(c) = this%biogeo_flux(c)%c13flux_vars%hr_col(c_l)
        c13flux_vars%fire_decomp_closs(c) = this%biogeo_flux(c)%c13flux_vars%fire_decomp_closs_col(c_l)
      endif
      if(use_c14_betr)then
        c14flux_vars%hr(c) = this%biogeo_flux(c)%c14flux_vars%hr_col(c_l)
        c14flux_vars%fire_decomp_closs(c) = this%biogeo_flux(c)%c14flux_vars%fire_decomp_closs_col(c_l)
      endif

      !recollect  nitrifications, nitrifier-N2O loss, denitrifications
      n14flux_vars%f_nit(c) = this%biogeo_flux(c)%n14flux_vars%f_nit_col(c_l)
      n14flux_vars%f_denit(c)= this%biogeo_flux(c)%n14flux_vars%f_denit_col(c_l)
      n14flux_vars%denit(c)= n14flux_vars%f_denit(c)
      n14flux_vars%f_n2o_nit(c)=this%biogeo_flux(c)%n14flux_vars%f_n2o_nit_col(c_l)

      !hydraulic loss
      n14flux_vars%smin_nh4_leached(c)= &
        this%biogeo_flux(c)%n14flux_vars%smin_nh4_leached_col(c_l) + &
        this%biogeo_flux(c)%n14flux_vars%smin_nh4_qdrain_col(c_l)

      n14flux_vars%smin_nh4_runoff(c)=this%biogeo_flux(c)%n14flux_vars%smin_nh4_runoff_col(c_l)

      n14flux_vars%som_n_leached(c) = &
        this%biogeo_flux(c)%n14flux_vars%som_n_leached_col(c_l) + &
        this%biogeo_flux(c)%n14flux_vars%som_n_qdrain_col(c_l)
      n14flux_vars%som_n_runoff(c) = this%biogeo_flux(c)%n14flux_vars%som_n_runoff_col(c_l)
      n14flux_vars%nh3_soi_flx(c) = this%biogeo_flux(c)%n14flux_vars%nh3_soi_flx_col(c_l)

      !the following is for consistency with the ALM definitation, which computes
      !som_n_leached_col as a numerical roundoff
      n14flux_vars%som_n_leached(c) = - n14flux_vars%som_n_leached(c)

      !fire loss
      n14flux_vars%fire_decomp_nloss(c) = this%biogeo_flux(c)%n14flux_vars%fire_decomp_nloss_col(c_l)
      n14flux_vars%supplement_to_sminn(c) = this%biogeo_flux(c)%n14flux_vars%supplement_to_sminn_col(c_l)
      !no nh4 volatilization and runoff/leaching loss at this moment
      !the following is to ensure mass balance, with an attempt to overcome some issue in cpl bypass
      n14flux_vars%ndep_to_sminn(c) = n14flux_vars%ndep_to_sminn_nh3(c)+n14flux_vars%ndep_to_sminn_no3(c)


      p31flux_vars%supplement_to_sminp(c) = this%biogeo_flux(c)%p31flux_vars%supplement_to_sminp_col(c_l)
      p31flux_vars%secondp_to_occlp(c) = this%biogeo_flux(c)%p31flux_vars%secondp_to_occlp_col(c_l)
      p31flux_vars%fire_decomp_ploss(c) = this%biogeo_flux(c)%p31flux_vars%fire_decomp_ploss_col(c_l)
      p31flux_vars%som_p_leached(c) = &
        this%biogeo_flux(c)%p31flux_vars%som_p_leached_col(c_l) + &
        this%biogeo_flux(c)%p31flux_vars%som_p_qdrain_col(c_l)
      p31flux_vars%som_p_runoff(c) = this%biogeo_flux(c)%p31flux_vars%som_p_runoff_col(c_l)
      !the following is for consistency with the ALM definitation, which computes
      !som_p_leached_col as a numerical roundoff
      p31flux_vars%som_p_leached(c) = -p31flux_vars%som_p_leached(c)
      p31flux_vars%primp_to_labilep(c) = this%biogeo_flux(c)%p31flux_vars%pflx_minp_weathering_po4_col(c_l)

      !recollect soil organic carbon, soil organic nitrogen, and soil organic phosphorus
      c12state_vars%cwdc(c) = this%biogeo_state(c)%c12state_vars%cwdc_col(c_l)
      c12state_vars%totlitc(c) = this%biogeo_state(c)%c12state_vars%totlitc_col(c_l)
      c12state_vars%totsomc(c) = this%biogeo_state(c)%c12state_vars%totsomc_col(c_l)
      c12state_vars%totlitc_1m(c) = this%biogeo_state(c)%c12state_vars%totlitc_1m_col(c_l)
      c12state_vars%totsomc_1m(c) = this%biogeo_state(c)%c12state_vars%totsomc_1m_col(c_l)

      if(use_c13_betr)then
        c13state_vars%cwdc(c) = this%biogeo_state(c)%c13state_vars%cwdc_col(c_l)
        c13state_vars%totlitc(c) = this%biogeo_state(c)%c13state_vars%totlitc_col(c_l)
        c13state_vars%totsomc(c) = this%biogeo_state(c)%c13state_vars%totsomc_col(c_l)
        c13state_vars%totlitc_1m(c) = this%biogeo_state(c)%c13state_vars%totlitc_1m_col(c_l)
        c13state_vars%totsomc_1m(c) = this%biogeo_state(c)%c13state_vars%totsomc_1m_col(c_l)
      endif
      if(use_c14_betr)then
        c14state_vars%cwdc(c) = this%biogeo_state(c)%c14state_vars%cwdc_col(c_l)
        c14state_vars%totlitc(c) = this%biogeo_state(c)%c14state_vars%totlitc_col(c_l)
        c14state_vars%totsomc(c) = this%biogeo_state(c)%c14state_vars%totsomc_col(c_l)
        c14state_vars%totlitc_1m(c) = this%biogeo_state(c)%c14state_vars%totlitc_1m_col(c_l)
        c14state_vars%totsomc_1m(c) = this%biogeo_state(c)%c14state_vars%totsomc_1m_col(c_l)
      endif
      n14state_vars%cwdn(c) = this%biogeo_state(c)%n14state_vars%cwdn_col(c_l)
      n14state_vars%totlitn(c) = this%biogeo_state(c)%n14state_vars%totlitn_col(c_l)
      n14state_vars%totsomn(c) = this%biogeo_state(c)%n14state_vars%totsomn_col(c_l)
      n14state_vars%totlitn_1m(c) = this%biogeo_state(c)%n14state_vars%totlitn_1m_col(c_l)
      n14state_vars%totsomn_1m(c) = this%biogeo_state(c)%n14state_vars%totsomn_1m_col(c_l)

      p31state_vars%cwdp(c) = this%biogeo_state(c)%p31state_vars%cwdp_col(c_l)
      p31state_vars%totlitp(c) = this%biogeo_state(c)%p31state_vars%totlitp_col(c_l)
      p31state_vars%totsomp(c) = this%biogeo_state(c)%p31state_vars%totsomp_col(c_l)
      p31state_vars%totlitp_1m(c) = this%biogeo_state(c)%p31state_vars%totlitp_1m_col(c_l)
      p31state_vars%totsomp_1m(c) = this%biogeo_state(c)%p31state_vars%totsomp_1m_col(c_l)

      !recollect inorganic nitrogen (smin_nh4, smin_no3), and inorganic phosphorus (disolvable and protected)
      n14state_vars%sminn(c) = this%biogeo_state(c)%n14state_vars%sminn_col(c_l)
      n14state_vars%smin_nh4(c)=this%biogeo_state(c)%n14state_vars%sminn_nh4_col(c_l)
      n14state_vars%smin_no3(c)=this%biogeo_state(c)%n14state_vars%sminn_no3_col(c_l)

      p31state_vars%sminp(c) = this%biogeo_state(c)%p31state_vars%sminp_col(c_l)
      p31state_vars%occlp(c) = this%biogeo_state(c)%p31state_vars%occlp_col(c_l)

      if(index(reaction_method,'ecacnp')/=0 .or. index(reaction_method,'ch4soil')/=0)then
        c12state_vars%som1c(c) = this%biogeo_state(c)%c12state_vars%som1c_col(c_l)
        c12state_vars%som2c(c) = this%biogeo_state(c)%c12state_vars%som2c_col(c_l)
        c12state_vars%som3c(c) = this%biogeo_state(c)%c12state_vars%som3c_col(c_l)
        if(use_c13_betr)then
          c13state_vars%som1c(c) = this%biogeo_state(c)%c13state_vars%som1c_col(c_l)
          c13state_vars%som2c(c) = this%biogeo_state(c)%c13state_vars%som2c_col(c_l)
          c13state_vars%som3c(c) = this%biogeo_state(c)%c13state_vars%som3c_col(c_l)
        endif
        if(use_c14_betr)then
          c14state_vars%som1c(c) = this%biogeo_state(c)%c14state_vars%som1c_col(c_l)
          c14state_vars%som2c(c) = this%biogeo_state(c)%c14state_vars%som2c_col(c_l)
          c14state_vars%som3c(c) = this%biogeo_state(c)%c14state_vars%som3c_col(c_l)
        endif
        n14state_vars%som1n(c) = this%biogeo_state(c)%n14state_vars%som1n_col(c_l)
        n14state_vars%som2n(c) = this%biogeo_state(c)%n14state_vars%som2n_col(c_l)
        n14state_vars%som3n(c) = this%biogeo_state(c)%n14state_vars%som3n_col(c_l)

        p31state_vars%som1p(c) = this%biogeo_state(c)%p31state_vars%som1p_col(c_l)
        p31state_vars%som2p(c) = this%biogeo_state(c)%p31state_vars%som2p_col(c_l)
        p31state_vars%som3p(c) = this%biogeo_state(c)%p31state_vars%som3p_col(c_l)

      endif
    enddo
    do fc = 1, num_surfc
      c = filter_soilc(fc)

      p31flux_vars%sminp_runoff(c)=this%biogeo_flux(c)%p31flux_vars%sminp_runoff_col(c_l)
      p31flux_vars%sminp_leached(c) = &
        this%biogeo_flux(c)%p31flux_vars%sminp_leached_col(c_l) + &
        this%biogeo_flux(c)%p31flux_vars%sminp_qdrain_col(c_l)
    enddo
  endif
  do fc = 1, num_surfc
    c = filter_soilc(fc)
    n14flux_vars%smin_no3_leached(c)= &
        this%biogeo_flux(c)%n14flux_vars%smin_no3_leached_col(c_l) + &
        this%biogeo_flux(c)%n14flux_vars%smin_no3_qdrain_col(c_l)
    n14flux_vars%smin_no3_runoff(c)=this%biogeo_flux(c)%n14flux_vars%smin_no3_runoff_col(c_l)
  enddo
  do j = 1,betr_nlevtrc_soil
    do fc = 1, num_surfc
      c =filter_soilc(fc)
      n14state_vars%smin_no3_vr(c,j) = this%biogeo_state(c)%n14state_vars%sminn_no3_vr_col(c_l,j)
    enddo
  enddo
  return
  if (index(reaction_method,'v1eca')/=0 ) then
    associate(                                                             &
    c12_decomp_cpools_vr_col =>   c12state_vars%decomp_cpools_vr   , &
    decomp_npools_vr_col =>   n14state_vars%decomp_npools_vr    , &
    decomp_ppools_vr_col =>   p31state_vars%decomp_ppools_vr  ,  &
    solutionp_vr        => p31state_vars%solutionp_vr          &

    )
    do kk = 1, 7
      do j = 1,betr_nlevtrc_soil
        do fc = 1, num_surfc
          c =filter_soilc(fc)
          c12_decomp_cpools_vr_col(c,j,kk)=this%biogeo_state(c)%c12state_vars%decomp_cpools_vr(c_l,j,kk)
          decomp_npools_vr_col(c,j,kk) = this%biogeo_state(c)%n14state_vars%decomp_npools_vr(c_l,j,kk)
          decomp_ppools_vr_col(c,j,kk) = this%biogeo_state(c)%p31state_vars%decomp_ppools_vr(c_l,j,kk)
        enddo
      enddo
    enddo
    do j = 1,betr_nlevtrc_soil
      do fc = 1, num_surfc
        c =filter_soilc(fc)
        solutionp_vr(c,j) = this%biogeo_state(c)%p31state_vars%sminp_vr_col(c_l,j)
      enddo
    enddo
    end associate
  endif

  end subroutine ALMBetrPlantSoilBGCRecv
  !------------------------------------------------------------------------

  subroutine ALMCalcDewSubFlux(this,  &
       bounds, col, num_hydrologyc, filter_soilc_hydrologyc)
   !DESCRIPTION
    ! Calculate tracer flux from dew or/and sublimation
    !External interface called by ALM

    use WaterfluxType   , only : waterflux_type
    use elm_varcon      , only : denh2o,spval
    use landunit_varcon , only : istsoil, istcrop
    use betr_decompMod  , only : betr_bounds_type
    implicit none
    !ARGUMENTS
    class(betr_simulation_alm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds
    type(column_type)               , intent(in)    :: col ! column type
    integer                         , intent(in)    :: num_hydrologyc ! number of column soil points in column filter_soilc
    integer                         , intent(in)    :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points

    !temporary variables
    type(betr_bounds_type)     :: betr_bounds
    integer :: fc, c

    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col)
    do fc = 1, num_hydrologyc
      c = filter_soilc_hydrologyc(fc)
      if(.not. this%active_col(c))cycle
      call this%betr(c)%calc_dew_sub_flux(this%betr_time,           &
         betr_bounds, this%betr_col(c), this%num_surfc, this%filter_soilc, &
        this%biophys_forc(c), this%betr(c)%tracers, this%betr(c)%tracerfluxes, this%betr(c)%tracerstates)
    enddo
  end subroutine ALMCalcDewSubFlux

  !------------------------------------------------------------------------
  subroutine ALMCalcSmpL(this, bounds, lbj, ubj, numf, filter, t_soisno, &
     soilstate_vars, waterstate_vars, soil_water_retention_curve)
  !DESCRIPTION
  ! calculate soil suction potential
  !
  !USES
  use SoilStateType              , only : soilstate_type
  use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
  use elm_varcon                 , only : grav,hfus,tfrz
  implicit none
  !ARGUMENTS
  class(betr_simulation_alm_type), intent(inout)  :: this
  type(bounds_type)                      , intent(in)    :: bounds  ! bounds
  integer                                , intent(in)    :: lbj, ubj                                          ! lower and upper bounds, make sure they are > 0
  integer                                , intent(in)    :: numf                                              ! number of columns in column filter
  integer                                , intent(in)    :: filter(:)                                         ! column filter
  real(r8)                               , intent(in)    :: t_soisno(bounds%begc: , lbj: )                    ! soil temperature
  type(soilstate_type)                   , intent(in)    :: soilstate_vars
  type(waterstate_type)                  , intent(inout) :: waterstate_vars
  class(soil_water_retention_curve_type) , intent(in)    :: soil_water_retention_curve

  !local variables
  real(r8) :: s_node
  integer  :: fc, c, j
  real(r8) :: dsmpds_top
  SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(mod_filename,__LINE__))

  ! remove compiler warnings
  if (this%num_surfc > 0) continue

  associate(                                                     & !
    h2osoi_vol        =>    waterstate_vars%h2osoi_vol     , & ! Input:  [real(r8) (:,:) ]  volumetric soil moisture
    smp_l             =>    waterstate_vars%smp_l          , & ! Output: [real(r8) (:,:) ]  soil suction (mm)
    bsw               =>    soilstate_vars%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
    watsat            =>    soilstate_vars%watsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
    sucsat            =>    soilstate_vars%sucsat_col            & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
  )

  do j = lbj, ubj
    do fc = 1, numf
      c = filter(fc)
      if(.not. this%active_col(c))cycle
      if(j==1)then
        if(t_soisno(c,j)<tfrz)then
          smp_l(c,j)= -hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
        else
          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_l(c,j), dsmpds_top)
!  the following will be implemented later.
!          call soil_water_retention_curve%soil_hk(hksat, 1._r8, s_node, bsw(c,j), hk)
        endif

!        this%biophys_forc(c)%Dw_hk(c_l) = hk *
      else
        if(t_soisno(c,j)<tfrz)then
          smp_l(c,j)= -hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
        else
          s_node = max(h2osoi_vol(c,j)/watsat(c,j), 0.01_r8)
          call soil_water_retention_curve%soil_suction(sucsat(c,j), s_node, bsw(c,j), smp_l(c,j))
        endif
      endif
    enddo
  enddo
  end associate
  end subroutine ALMCalcSmpL

  !------------------------------------------------------------------------
  subroutine ALMSetBiophysForcing(this, bounds, col, pft, carbonflux_vars, pf_carbonflux_vars, &
    waterstate_vars, waterflux_vars, pf_waterflux_vars, temperature_vars, pf_temperature_vars, &
    soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
    chemstate_vars, soilstate_vars, cnstate_vars, carbonstate_vars, phosphorusstate_vars)
  !DESCRIPTION
  !pass in biogeophysical variables for running betr
  !USES
  use SoilStateType     , only : soilstate_type
  use ChemStateType     , only : chemstate_type
  use atm2lndType       , only : atm2lnd_type
  use SoilHydrologyType , only : soilhydrology_type
  use CNStateType       , only : cnstate_type
  use CanopyStateType   , only : canopystate_type
  use elm_varpar        , only : nlevsno, nlevsoi
  use tracer_varcon     , only : reaction_method
  use tracer_varcon     , only : betr_nlevsoi
  use tracer_varcon     , only : catomw
  implicit none
  !ARGUMENTS
  class(betr_simulation_alm_type) , intent(inout)        :: this
  type(bounds_type)               , intent(in)           :: bounds
  type(patch_type)            , intent(in) :: pft
  type(column_type)           , intent(in)    :: col ! column type
  type(cnstate_type)          , optional, intent(in) :: cnstate_vars
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
  type(carbonstate_type)      , optional, intent(in) :: carbonstate_vars
  type(phosphorusstate_type)  , optional, intent(in) :: phosphorusstate_vars

  integer :: p, pi, c, j, c_l, pp
  integer :: npft_loc

  c_l=1
  if(present(phosphorusstate_vars))then
    !the following is used for setting P upon exiting spinup
    do j = 1, betr_nlevsoi
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle
        this%biophys_forc(c)%solutionp_vr_col(c_l,j) = phosphorusstate_vars%solutionp_vr(c,j)
        this%biophys_forc(c)%labilep_vr_col(c_l,j) = phosphorusstate_vars%labilep_vr(c,j)
        this%biophys_forc(c)%secondp_vr_col(c_l,j) = phosphorusstate_vars%secondp_vr(c,j)
        this%biophys_forc(c)%occlp_vr_col(c_l,j) =  phosphorusstate_vars%occlp_vr(c,j)
      enddo
    enddo
  endif

  if(present(carbonflux_vars)) then
    call this%BeTRSetBiophysForcing(bounds, col, pft, 1, betr_nlevsoi, &
    carbonflux_vars, pf_carbonflux_vars, waterstate_vars,  waterflux_vars, pf_waterflux_vars, &
    temperature_vars, pf_temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
    chemstate_vars, soilstate_vars)
    do j = 1, betr_nlevsoi
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle
        this%biophys_forc(c)%c12flx%rt_vr_col(c_l,j) = carbonflux_vars%rr_vr(c,j)
      enddo
    enddo
  else
    return
  endif
  if(present(cnstate_vars)  .and. present(pf_carbonflux_vars)) then
    associate(                                                &
      cn_scalar            => cnstate_vars%cn_scalar        , &
      cp_scalar            => cnstate_vars%cp_scalar        , &
      rootfr               => soilstate_vars%rootfr_patch     &
    )
    !the following will be ALM specific
    !big leaf model
    !set profiles autotrohpic respiration
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      npft_loc = ubound(pf_carbonflux_vars%rr,1)-lbound(pf_carbonflux_vars%rr,1)+1
      if(npft_loc /= col%npfts(c) .and. col%pfti(c) /= lbound(pf_carbonflux_vars%rr,1)) then
        do pi = 1, betr_maxpatch_pft
          this%biophys_forc(c)%rr_patch(pi,1:betr_nlevsoi) = 0._r8
        enddo
      else
        if(use_cn)then
          pp = 0
          do pi = 1, betr_maxpatch_pft
            if (pi <= col%npfts(c)) then
              p = col%pfti(c) + pi - 1
              if (pft%active(p)) then
                pp = pp + 1
                this%biophys_forc(c)%cn_scalar_patch(pp) = cn_scalar(p)
                this%biophys_forc(c)%cp_scalar_patch(pp) = cp_scalar(p)
                do j = 1, betr_nlevsoi
                  this%biophys_forc(c)%rr_patch(pp,j) = pf_carbonflux_vars%rr(p) * rootfr(p,j)
                enddo
              endif
            endif
          enddo
        else
          do pi = 1, betr_maxpatch_pft
            this%biophys_forc(c)%rr_patch(pi,1:betr_nlevsoi) = 0._r8
          enddo
        endif
      endif
    enddo
  end associate
  endif
  !dvgm
  if(trim(reaction_method)=='doc_dic')then
     c_l=1
     do j = 1, betr_nlevsoi
        do c = bounds%begc, bounds%endc
          if(.not. this%active_col(c))cycle
             !for simplicity, atomic weight of carbon is set to 12._r8 g/mol
           this%biophys_forc(c)%dic_prod_vr_col(c_l,j) = (carbonflux_vars%hr_vr(c,j) + &
                cnstate_vars%nfixation_prof_col(c,j)*carbonflux_vars%rr(c))/catomw
           this%biophys_forc(c)%doc_prod_vr_col(c_l,j) = (carbonstate_vars%decomp_cpools_vr(c,j,6) - &
                carbonstate_vars%decomp_som2c_vr(c,j))/this%betr_time%delta_time/catomw
        enddo
      enddo
  endif
  end subroutine ALMSetBiophysForcing

  !------------------------------------------------------------------------

  subroutine set_vegpar_calibration(this, betr_bounds, col, pft, num_surfc, filter_soilc, calibration_vpars)
  !DESCRIPTION
  !set kinetic parameters for column c
  use PlantMicKineticsMod, only : PlantMicKinetics_type
  use tracer_varcon      , only : reaction_method,natomw,patomw
  implicit none
  class(betr_simulation_alm_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: betr_bounds
  type(column_type)     , intent(in)    :: col ! column type
  type(patch_type)      , intent(in) :: pft
  integer, intent(in) :: num_surfc
  integer, intent(in) :: filter_soilc(:)
  type(calibration_type), intent(in) :: calibration_vpars

  integer :: j, fc, c, p, pi, pp, g

  associate(                                                                  &
    plant_nh4_vmax_scalar  => calibration_vpars%plant_nh4_vmax_scalar       , &
    plant_no3_vmax_scalar  => calibration_vpars%plant_no3_vmax_scalar       , &
    plant_p_vmax_scalar    => calibration_vpars%plant_p_vmax_scalar         , &
    plant_nh4_km_scalar    => calibration_vpars%plant_nh4_km_scalar         , &
    plant_no3_km_scalar    => calibration_vpars%plant_no3_km_scalar         , &
    plant_p_km_scalar      => calibration_vpars%plant_p_km_scalar             &
  )

  do fc = 1, num_surfc
    c = filter_soilc(fc)
    pp = 0
    do pi = 1, betr_maxpatch_pft
      if (pi <= col%npfts(c)) then
        p = col%pfti(c) + pi - 1
        g=pft%gridcell(p)
        if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
          pp = pp + 1
          do j =1, betr_bounds%ubj
            this%betr(c)%plantNutkinetics%plant_nh4_vmax_vr_patch(pp,j) = this%betr(c)%plantNutkinetics%plant_nh4_vmax_vr_patch(pp,j) * plant_nh4_vmax_scalar(g,pft%itype(p))
            this%betr(c)%plantNutkinetics%plant_no3_vmax_vr_patch(pp,j) = this%betr(c)%plantNutkinetics%plant_no3_vmax_vr_patch(pp,j) * plant_no3_vmax_scalar(g,pft%itype(p))
            this%betr(c)%plantNutkinetics%plant_p_vmax_vr_patch(pp,j) = this%betr(c)%plantNutkinetics%plant_p_vmax_vr_patch(pp,j) * plant_p_vmax_scalar(g,pft%itype(p))
            this%betr(c)%plantNutkinetics%plant_nh4_km_vr_patch(pp,j) = this%betr(c)%plantNutkinetics%plant_nh4_km_vr_patch(pp,j) * plant_nh4_km_scalar(g,pft%itype(p))
            this%betr(c)%plantNutkinetics%plant_no3_km_vr_patch(pp,j) = this%betr(c)%plantNutkinetics%plant_no3_km_vr_patch(pp,j) * plant_no3_km_scalar(g,pft%itype(p))
            this%betr(c)%plantNutkinetics%plant_p_km_vr_patch(pp,j) = this%betr(c)%plantNutkinetics%plant_p_km_vr_patch(pp,j) * plant_p_km_scalar(g,pft%itype(p))
          enddo
        endif
      endif
    enddo
  enddo
  end associate
  end subroutine set_vegpar_calibration
  !------------------------------------------------------------------------
  subroutine set_transient_kinetics_par(this, betr_bounds, col, pft, num_surfc, filter_soilc, PlantMicKinetics_vars)
  !DESCRIPTION
  !set kinetic parameters for column c
  use PlantMicKineticsMod, only : PlantMicKinetics_type
  use tracer_varcon      , only : reaction_method,natomw,patomw
  use clm_time_manager   , only : get_nstep
  use tracer_varcon      , only : lbcalib
  implicit none
  class(betr_simulation_alm_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: betr_bounds
  type(column_type)     , intent(in)    :: col ! column type
  type(patch_type)      , intent(in) :: pft
  integer, intent(in) :: num_surfc
  integer, intent(in) :: filter_soilc(:)
  type(PlantMicKinetics_type), intent(in) :: PlantMicKinetics_vars

  integer :: j, fc, c, p, pi, pp, c_l, val

  associate(      &
    plant_nh4_vmax_vr_patch => PlantMicKinetics_vars%plant_nh4_vmax_vr_patch, &
    plant_no3_vmax_vr_patch => PlantMicKinetics_vars%plant_no3_vmax_vr_patch, &
    plant_p_vmax_vr_patch   => PlantMicKinetics_vars%plant_p_vmax_vr_patch, &
    plant_nh4_km_vr_patch   => PlantMicKinetics_vars%plant_nh4_km_vr_patch, &
    plant_no3_km_vr_patch   => PlantMicKinetics_vars%plant_no3_km_vr_patch, &
    plant_p_km_vr_patch     => PlantMicKinetics_vars%plant_p_km_vr_patch , &
    plant_eff_ncompet_b_vr_patch => PlantMicKinetics_vars%plant_eff_ncompet_b_vr_patch , &
    plant_eff_pcompet_b_vr_patch => PlantMicKinetics_vars%plant_eff_pcompet_b_vr_patch , &
    minsurf_nh4_compet_vr_col => PlantMicKinetics_vars%minsurf_nh4_compet_vr_col, &
    minsurf_p_compet_vr_col => PlantMicKinetics_vars%minsurf_p_compet_vr_col , &
    plant_eff_frootc_vr_patch => PlantMicKinetics_vars%plant_eff_frootc_vr_patch, &
    decomp_eff_ncompet_b_vr_col => PlantMicKinetics_vars%decomp_eff_ncompet_b_vr_col, &
    decomp_eff_pcompet_b_vr_col => PlantMicKinetics_vars%decomp_eff_pcompet_b_vr_col, &
    nit_eff_ncompet_b_vr_col => PlantMicKinetics_vars%nit_eff_ncompet_b_vr_col, &
    den_eff_ncompet_b_vr_col => PlantMicKinetics_vars%den_eff_ncompet_b_vr_col, &
    dsolutionp_dt_vr_col   => PlantMicKinetics_vars%dsolutionp_dt_vr_col, &
    vmax_minsurf_p_vr_col => PlantMicKinetics_vars%vmax_minsurf_p_vr_col, &
    km_minsurf_p_vr_col   => PlantMicKinetics_vars%km_minsurf_p_vr_col, &
    km_minsurf_nh4_vr_col => PlantMicKinetics_vars%km_minsurf_nh4_vr_col, &
    dlabp_dt_vr_col       => PlantMicKinetics_vars%dlabp_dt_vr_col &
  )
  c_l = 1
  do fc = 1, num_surfc
    c = filter_soilc(fc)
    pp = 0
    val=1._r8

    do pi = 1, betr_maxpatch_pft
      if (pi <= col%npfts(c)) then
        p = col%pfti(c) + pi - 1
        if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
          pp = pp + 1
          do j =1, betr_bounds%ubj
            this%betr(c)%plantNutkinetics%plant_nh4_vmax_vr_patch(pp,j) = plant_nh4_vmax_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_no3_vmax_vr_patch(pp,j) = plant_no3_vmax_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_p_vmax_vr_patch(pp,j) = plant_p_vmax_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_nh4_km_vr_patch(pp,j) = plant_nh4_km_vr_patch(p,j)/natomw
            this%betr(c)%plantNutkinetics%plant_no3_km_vr_patch(pp,j) = plant_no3_km_vr_patch(p,j)/natomw
            this%betr(c)%plantNutkinetics%plant_p_km_vr_patch(pp,j) = plant_p_km_vr_patch(p,j)/patomw
            this%betr(c)%plantNutkinetics%plant_eff_ncompet_b_vr_patch(pp,j)=plant_eff_ncompet_b_vr_patch(p,j)/natomw
            this%betr(c)%plantNutkinetics%plant_eff_pcompet_b_vr_patch(pp,j)=plant_eff_pcompet_b_vr_patch(p,j)/patomw
            this%betr(c)%plantNutkinetics%plant_eff_frootc_vr_patch(pp,j) = plant_eff_frootc_vr_patch(p,j)
          enddo
        endif
      endif
    enddo
    this%betr(c)%nactpft = pp
    do j = 1, betr_bounds%ubj
      this%betr(c)%plantNutkinetics%minsurf_p_compet_vr_col(c_l,j) = minsurf_p_compet_vr_col(c,j)/patomw
      this%betr(c)%plantNutkinetics%minsurf_nh4_compet_vr_col(c_l,j) = minsurf_nh4_compet_vr_col(c,j)/natomw
    enddo
  enddo

  !the following parameters are specific to ECACNP, and I assume they are
  !grid specific as they currently used in alm-cnp.
  if(index(reaction_method,'ecacnp')/=0 .or. index(reaction_method, 'ch4soil')/=0 &
     .or. index(reaction_method, 'v1eca')/=0)then
    do j =1, betr_bounds%ubj
      do fc = 1, num_surfc
        c = filter_soilc(fc)
        this%betr(c)%plantNutkinetics%km_minsurf_p_vr_col(c_l,j)  = km_minsurf_p_vr_col(c,j)/patomw
        this%betr(c)%plantNutkinetics%km_minsurf_nh4_vr_col(c_l,j)= km_minsurf_nh4_vr_col(c,j)/natomw
      enddo
    enddo
    if(lbcalib)then
      do j =1, betr_bounds%ubj
        do fc = 1, num_surfc
          c = filter_soilc(fc)
          this%betr(c)%plantNutkinetics%km_decomp_p_vr_col(c_l,j) = PlantMicKinetics_vars%km_decomp_p_vr_col(c,j)/patomw
          this%betr(c)%plantNutkinetics%km_decomp_nh4_vr_col(c_l,j)=PlantMicKinetics_vars%km_decomp_nh4_vr_col(c,j)/natomw
          this%betr(c)%plantNutkinetics%km_decomp_no3_vr_col(c_l,j)=PlantMicKinetics_vars%km_decomp_no3_vr_col(c,j)/natomw
          this%betr(c)%plantNutkinetics%km_nit_nh4_vr_col(c_l,j) = PlantMicKinetics_vars%km_nit_nh4_vr_col(c,j)/natomw
          this%betr(c)%plantNutkinetics%km_den_no3_vr_col(c_l,j) = PlantMicKinetics_vars%km_den_no3_vr_col(c,j)/natomw
        enddo
      enddo
    endif
  endif
  if(index(reaction_method,'v1eca')/=0)then
    do j =1, betr_bounds%ubj
      do fc = 1, num_surfc
        c = filter_soilc(fc)
        this%betr(c)%plantNutkinetics%decomp_eff_ncompet_b_vr_col(c_l,j)= decomp_eff_ncompet_b_vr_col(c,j)/natomw
        this%betr(c)%plantNutkinetics%decomp_eff_pcompet_b_vr_col(c_l,j)= decomp_eff_pcompet_b_vr_col(c,j)/patomw
        this%betr(c)%plantNutkinetics%nit_eff_ncompet_b_vr_col(c_l,j)   = nit_eff_ncompet_b_vr_col(c,j)/natomw
        this%betr(c)%plantNutkinetics%den_eff_ncompet_b_vr_col(c_l,j)   = den_eff_ncompet_b_vr_col(c,j)/natomw
        this%betr(c)%plantNutkinetics%km_nit_nh4_vr_col(c_l,j) = PlantMicKinetics_vars%km_nit_nh4_vr_col(c,j)/natomw
        this%betr(c)%plantNutkinetics%km_den_no3_vr_col(c_l,j) = PlantMicKinetics_vars%km_den_no3_vr_col(c,j)/natomw
        this%betr(c)%plantNutkinetics%dsolutionp_dt_vr_col(c_l,j)       = dsolutionp_dt_vr_col(c,j)/patomw   ! g/m2/s
        this%betr(c)%plantNutkinetics%vmax_minsurf_p_vr_col(c_l,j)      = vmax_minsurf_p_vr_col(c,j)/patomw   ! g/m3
        this%betr(c)%plantNutkinetics%dlabp_dt_vr_col(c_l,j)            = dlabp_dt_vr_col(c,j)/patomw
      enddo
    enddo
  endif
  end associate
  end subroutine set_transient_kinetics_par


!-------------------------------------------------------------------------------
  subroutine ALMOutLoopSoilBGC(this, bounds,  col, pft)

  implicit none
    ! !ARGUMENTS :
    class(betr_simulation_alm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds ! bounds
    type(column_type)               , intent(in)    :: col ! column type
    type(patch_type)                , intent(in)    :: pft

    !TEMPORARY VARIABLES
    type(betr_bounds_type)     :: betr_bounds
    integer :: c, c_l, begc_l, endc_l


    call this%BeTRSetBounds(betr_bounds)

    call this%BeTRSetcps(bounds, col, pft)

    c_l = 1; begc_l = betr_bounds%begc; endc_l=betr_bounds%endc;

    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      !this%betr(c)%tracers%debug=col%debug_flag(c)

      call this%betr(c)%OutLoopBGC(this%betr_time, betr_bounds, this%betr_col(c), &
         this%betr_pft(c), this%num_surfc, this%filter_soilc, this%num_soilp, this%filter_soilp, &
         this%biophys_forc(c), this%biogeo_flux(c), this%biogeo_state(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err(),c)
        exit
      endif

    enddo
  end subroutine ALMOutLoopSoilBGC

!-------------------------------------------------------------------------------
  function do_bgc_type(this, type_char)result(ans)
  use betr_ctrl, only : bgc_type
  implicit none
  class(betr_simulation_alm_type) , intent(inout) :: this
  character(len=*), intent(in) :: type_char


  logical :: ans

  ans = index(bgc_type,type_char)/=0

  end function do_bgc_type


!-------------------------------------------------------------------------------
  subroutine ALMEnterOutLoopBGC(this, bounds, col, pft, num_surfc, filter_soilc, &
   c12_cstate_vars, c12_cflx_vars, c13_cstate_vars, c14_cstate_vars,  &
   nitrogenstate_vars,  phosphorusstate_vars, PlantMicKinetics_vars)

  use tracer_varcon   , only : nlevtrc_soil  => betr_nlevtrc_soil
  use PlantMicKineticsMod, only : PlantMicKinetics_type
  implicit none
  class(betr_simulation_alm_type) , intent(inout) :: this
  type(bounds_type)               , intent(in)    :: bounds ! bounds
  type(column_type)               , intent(in)    :: col ! column type
  type(patch_type)                , intent(in)    :: pft
  integer                         , intent(in)    :: num_surfc        ! number of soil columns in filter
  integer                         , intent(in)    :: filter_soilc(:)  ! filter for soil columns
  type(PlantMicKinetics_type)     , intent(in)    :: PlantMicKinetics_vars
  type(carbonstate_type)          , intent(inout) :: c12_cstate_vars
  type(carbonflux_type)           , intent(inout) :: c12_cflx_vars
  type(carbonstate_type)          , intent(inout) :: c13_cstate_vars
  type(carbonstate_type)          , intent(inout) :: c14_cstate_vars
  type(nitrogenstate_type)        , intent(inout) :: nitrogenstate_vars
  type(phosphorusstate_type)      , intent(inout) :: phosphorusstate_vars

  !temporary variables
  type(betr_bounds_type) :: betr_bounds
  integer :: kk, c, j, c_l, fc

  associate(                                                             &
  decomp_k                 => c12_cflx_vars%decomp_k                   , &
  t_scalar                 => c12_cflx_vars%t_scalar                   , &
  w_scalar                 => c12_cflx_vars%w_scalar                   , &
  c12_decomp_cpools_vr_col =>   c12_cstate_vars%decomp_cpools_vr       , &
  decomp_npools_vr_col     =>   nitrogenstate_vars%decomp_npools_vr    , &
  decomp_ppools_vr_col     =>   phosphorusstate_vars%decomp_ppools_vr  , &
  smin_nh4_vr              =>   nitrogenstate_vars%smin_nh4_vr         , & ! Input:  [real(r8) (:,:)  ]  (gN/m3) soil mineral NH4 pool
  smin_no3_vr              =>   nitrogenstate_vars%smin_no3_vr         , & ! Input:  [real(r8) (:,:)  ]  (gN/m3) soil mineral NO3 pool
  solutionp_vr             => phosphorusstate_vars%solutionp_vr          & !
  )
  c_l=1
  call this%BeTRSetBounds(betr_bounds)
  do j = 1,nlevtrc_soil
    do fc = 1, num_surfc
      c = filter_soilc(fc)
      this%biophys_forc(c)%c12flx%in_t_scalar(c_l,j) = t_scalar(c,j)
      this%biophys_forc(c)%c12flx%in_w_scalar(c_l,j) = w_scalar(c,j)
      this%biophys_forc(c)%n14flx%in_sminn_no3_vr_col(c_l,j) = smin_no3_vr(c,j)
      this%biophys_forc(c)%n14flx%in_sminn_nh4_vr_col(c_l,j) = smin_nh4_vr(c,j)
      this%biophys_forc(c)%p31flx%in_sminp_vr_col(c_l,j) = solutionp_vr(c,j)
    enddo
  enddo

  do kk = 1, 7
    do j = 1,nlevtrc_soil
      do fc = 1, num_surfc
        c = filter_soilc(fc)
        this%biophys_forc(c)%c12flx%in_decomp_cpools_vr_col(c_l,j,kk)=c12_decomp_cpools_vr_col(c,j,kk)
        this%biophys_forc(c)%n14flx%in_decomp_npools_vr_col(c_l,j,kk)=decomp_npools_vr_col(c,j,kk)
        this%biophys_forc(c)%p31flx%in_decomp_ppools_vr_col(c_l,j,kk)=decomp_ppools_vr_col(c,j,kk)
        this%biogeo_flux(c)%c12flux_vars%decomp_k(c_l,j,kk)=decomp_k(c,j,kk)
      enddo
    enddo
  enddo

  !set autotrophic respiration

  !set kinetic parameters
  call this%set_transient_kinetics_par(betr_bounds, col, pft, num_surfc, filter_soilc, PlantMicKinetics_vars)

  end associate
  end subroutine ALMEnterOutLoopBGC

!-------------------------------------------------------------------------------
  subroutine ALMExitOutLoopBGC(this, bounds, col, pft, &
    c12_cstate_vars, c12_cflx_vars, &
    c13_cstate_vars, c13_cflx_vars, &
    c14_cstate_vars, c14_cflx_vars, &
    nitrogenstate_vars, nitrogenflux_vars, pf_nitrogenflux_vars, &
    phosphorusstate_vars, phosphorusflux_vars, pf_phosphorusflux_vars)

  use tracer_varcon   , only : nlevtrc_soil  => betr_nlevtrc_soil
  use clm_time_manager    , only : get_nstep
  implicit none
  class(betr_simulation_alm_type) , intent(inout) :: this
  type(bounds_type)               , intent(in)    :: bounds ! bounds
  type(column_type)               , intent(in)    :: col ! column type
  type(patch_type)                , intent(in)    :: pft
  type(carbonstate_type)          , intent(inout) :: c12_cstate_vars
  type(carbonflux_type)           , intent(inout) :: c12_cflx_vars
  type(carbonstate_type)          , intent(inout) :: c13_cstate_vars
  type(carbonflux_type)           , intent(inout) :: c13_cflx_vars
  type(carbonstate_type)          , intent(inout) :: c14_cstate_vars
  type(carbonflux_type)           , intent(inout) :: c14_cflx_vars
  type(nitrogenstate_type)        , intent(inout) :: nitrogenstate_vars
  type(nitrogenflux_type)         , intent(inout) :: nitrogenflux_vars
  type(phosphorusstate_type)      , intent(inout) :: phosphorusstate_vars
  type(phosphorusflux_type)       , intent(inout) :: phosphorusflux_vars
  type(pf_nitrogenflux_type)  , intent(inout) :: pf_nitrogenflux_vars
  type(pf_phosphorusflux_type), intent(inout) :: pf_phosphorusflux_vars

  integer :: kk, c, j, fc, c_l, p, pi
  associate(                                                         &
  decomp_k                 => c12_cflx_vars%decomp_k               , &
  c12_decomp_cpools_vr_col =>   c12_cstate_vars%decomp_cpools_vr   , &
  decomp_npools_vr_col =>   nitrogenstate_vars%decomp_npools_vr    , &
  decomp_ppools_vr_col =>   phosphorusstate_vars%decomp_ppools_vr    &

  )
  c_l=1
  do kk = 1, 7
    do j = 1,nlevtrc_soil
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle

        decomp_k(c,j,kk) = this%biogeo_flux(c)%c12flux_vars%decomp_k(c_l,j,kk)
        c12_decomp_cpools_vr_col(c,j,kk)=this%biogeo_state(c)%c12state_vars%decomp_cpools_vr(c_l,j,kk)
        decomp_npools_vr_col(c,j,kk) = this%biogeo_state(c)%n14state_vars%decomp_npools_vr(c_l,j,kk)
        decomp_ppools_vr_col(c,j,kk) = this%biogeo_state(c)%p31state_vars%decomp_ppools_vr(c_l,j,kk)
      enddo
    enddo
  enddo
  !extract plant nutrient uptake fluxes, soil respiration, denitrification, nitrification
  !
  do c = bounds%begc, bounds%endc
    if(.not. this%active_col(c))cycle
    pi = 0
    do p = col%pfti(c), col%pftf(c)
      if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
        pi = pi + 1
        pf_nitrogenflux_vars%smin_nh4_to_plant(p) = this%biogeo_flux(c)%n14flux_vars%smin_nh4_to_plant_patch(pi)
        pf_nitrogenflux_vars%smin_no3_to_plant(p) = this%biogeo_flux(c)%n14flux_vars%smin_no3_to_plant_patch(pi)
        pf_nitrogenflux_vars%smin_nh4_to_plant_vr(p,1:nlevtrc_soil) = this%biogeo_flux(c)%n14flux_vars%smin_nh4_to_plant_vr_patch(pi,1:nlevtrc_soil)
        pf_nitrogenflux_vars%smin_no3_to_plant_vr(p,1:nlevtrc_soil) = this%biogeo_flux(c)%n14flux_vars%smin_no3_to_plant_vr_patch(pi,1:nlevtrc_soil)

        pf_nitrogenflux_vars%sminn_to_plant(p) = pf_nitrogenflux_vars%smin_nh4_to_plant(p) + pf_nitrogenflux_vars%smin_no3_to_plant(p)
        pf_phosphorusflux_vars%sminp_to_plant(p)  = this%biogeo_flux(c)%p31flux_vars%sminp_to_plant_patch(pi)
      else
        pf_nitrogenflux_vars%smin_nh4_to_plant(p) = 0._r8
        pf_nitrogenflux_vars%smin_no3_to_plant(p) = 0._r8
        pf_phosphorusflux_vars%sminp_to_plant(p) = 0._r8
        pf_nitrogenflux_vars%sminn_to_plant(p) = 0._r8
        pf_nitrogenflux_vars%smin_nh4_to_plant_vr(p,1:nlevtrc_soil) = 0._r8
        pf_nitrogenflux_vars%smin_no3_to_plant_vr(p,1:nlevtrc_soil) = 0._r8
      endif
    enddo
    nitrogenflux_vars%actual_immob(c)=0._r8
!    print*,'immob',nitrogenflux_vars%actual_immob(c)
  enddo

  do j = 1,nlevtrc_soil
    do c = bounds%begc, bounds%endc
      if(.not. this%active_col(c))cycle
      phosphorusflux_vars%col_plant_pdemand_vr(c,j)  = this%biogeo_flux(c)%p31flux_vars%col_plant_pdemand_vr(c_l,j)
      nitrogenflux_vars%f_denit_vr(c,j)        = this%biogeo_flux(c)%n14flux_vars%f_denit_vr_col(c_l,j)
      nitrogenflux_vars%f_n2o_denit_vr(c,j)    = this%biogeo_flux(c)%n14flux_vars%f_n2o_denit_vr_col(c_l,j)
      nitrogenflux_vars%f_nit_vr(c,j)          = this%biogeo_flux(c)%n14flux_vars%f_nit_vr_col(c_l,j)
      nitrogenflux_vars%f_n2o_nit_vr(c,j)      = this%biogeo_flux(c)%n14flux_vars%f_n2o_nit_vr_col(c_l,j)
      phosphorusflux_vars%adsorb_to_labilep_vr(c,j)= this%biogeo_flux(c)%p31flux_vars%adsorb_to_labilep_vr_col(c_l,j)
      c12_cflx_vars%hr_vr(c,j)                 = this%biogeo_flux(c)%c12flux_vars%hr_vr_col(c_l,j)
      c12_cflx_vars%phr_vr(c,j)                = this%biogeo_flux(c)%c12flux_vars%phr_vr_col(c_l,j)
      c12_cflx_vars%o_scalar(c,j)              = this%biogeo_flux(c)%c12flux_vars%o_scalar_col(c_l,j)
      nitrogenstate_vars%smin_nh4_vr(c,j)      = this%biogeo_state(c)%n14state_vars%sminn_nh4_vr_col(c_l,j)
      nitrogenstate_vars%smin_no3_vr(c,j)      = this%biogeo_state(c)%n14state_vars%sminn_no3_vr_col(c_l,j)
      nitrogenflux_vars%supplement_to_sminn_vr(c,j) = this%biogeo_flux(c)%n14flux_vars%supplement_to_sminn_vr_col(c_l,j)
      nitrogenflux_vars%smin_nh4_to_plant_vr(c,j) = this%biogeo_flux(c)%n14flux_vars%smin_nh4_to_plant_vr_col(c_l,j)
      nitrogenflux_vars%smin_no3_to_plant_vr(c,j) = this%biogeo_flux(c)%n14flux_vars%smin_no3_to_plant_vr_col(c_l,j)
      phosphorusflux_vars%sminp_to_plant_vr(c,j) = this%biogeo_flux(c)%p31flux_vars%sminp_to_plant_vr_col(c_l,j)
      phosphorusflux_vars%supplement_to_sminp_vr(c,j) =this%biogeo_flux(c)%p31flux_vars%supplement_to_sminp_vr_col(c_l,j)
      phosphorusflux_vars%net_mineralization_p_vr(c,j) = this%biogeo_flux(c)%p31flux_vars%net_mineralization_p_vr_col(c_l,j)
      nitrogenflux_vars%actual_immob(c)= nitrogenflux_vars%actual_immob(c) + col%dz(c,j)*&
         (this%biogeo_flux(c)%n14flux_vars%smin_nh4_immob_vr_col(c_l,j) + this%biogeo_flux(c)%n14flux_vars%smin_no3_immob_vr_col(c_l,j))
      c12_cflx_vars%somhr(c) = c12_cflx_vars%somhr(c)  + col%dz(c,j)*&
         this%biogeo_flux(c)%c12flux_vars%somhr_vr_col(c_l,j)
      c12_cflx_vars%lithr(c) = c12_cflx_vars%lithr(c)  + col%dz(c,j)*&
         this%biogeo_flux(c)%c12flux_vars%lithr_vr_col(c_l,j)
      c12_cflx_vars%cwdc_hr(c) = c12_cflx_vars%cwdc_hr(c)  + col%dz(c,j)*&
         this%biogeo_flux(c)%c12flux_vars%cwdhr_vr_col(c_l,j)
    enddo
  enddo
  end associate

  end subroutine ALMExitOutLoopBGC
end module BeTRSimulationALM
