
module BetrType

#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  !  subroutines for betr application
  !
  !  !USES:
  use bshr_kind_mod            , only : r8 => shr_kind_r8
  use bshr_log_mod             , only : errMsg => shr_log_errMsg
  use BeTR_decompMod           , only : bounds_type  => betr_bounds_type
  use betr_ctrl                , only : iulog  => biulog
  use betr_constants           , only : betr_string_length
  use BetrStatusType           , only : betr_status_type
  use tracer_varcon            , only : nlevsoi  => betr_nlevsoi
  use tracer_varcon            , only : nlevsno => betr_nlevsno
  use betr_varcon              , only : spval => bspval
  use BGCReactionsMod          , only : bgc_reaction_type
  use PlantSoilBGCMod          , only : plant_soilbgc_type
  use BeTRTracerType           , only : betrtracer_type
  use TracerCoeffType          , only : TracerCoeff_type
  use TracerFluxType           , only : TracerFlux_type
  use TracerStateType          , only : TracerState_type
  use tracerboundarycondType   , only : tracerboundarycond_type
  use BeTR_aerocondType        , only : betr_aerecond_type
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use BeTR_biogeoStateType     , only : betr_biogeo_state_type
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  use PlantNutKineticsMod      , only : PlantNutKinetics_type
  implicit none

  private

  character(len=*), parameter :: filename = &
       __FILE__

  type, public :: betr_type
     ! namelist control variables
     character(len=betr_string_length) :: reaction_method
     logical                           :: diffusion_on = .true.
     logical                           :: advection_on = .true.
     logical                           :: reaction_on  = .true.
     logical                           :: ebullition_on= .true.
     integer                           :: nactpft
     logical                           :: active_soibgc
     ! internal types
     class(bgc_reaction_type)      , allocatable, public :: bgc_reaction
     class(plant_soilbgc_type)     , allocatable, public :: plant_soilbgc

     type(BeTRtracer_type)         , public              :: tracers
     type(TracerCoeff_type)        , public              :: tracercoeffs
     type(TracerFlux_type)         , public              :: tracerfluxes
     type(TracerState_type)        , public              :: tracerstates
     type(tracerboundarycond_type) , public              :: tracerboundaryconds
     type(PlantNutKinetics_type)   , public              :: plantNutkinetics
     type(betr_aerecond_type)      , private             :: aereconds

     real(r8)                      , private, pointer    :: h2osoi_liq_copy(:,:) => null()
     real(r8)                      , private, pointer    :: h2osoi_ice_copy(:,:) => null()

     ! FIXME(bja, 201603) replace LSM specific types!

   contains
     procedure, public  :: Init
     procedure, public  :: step_without_drainage
     procedure, public  :: step_with_drainage
     procedure, public  :: calc_dew_sub_flux
     procedure, public  :: betr_lsm_flux_state_sendback
     procedure, public  :: pre_diagnose_soilcol_water_flux
     procedure, public  :: diagnose_advect_water_flux
     procedure, public  :: diagnose_drainage_water_flux
     procedure, public  :: diagnose_dtracer_freeze_thaw
     procedure, public  :: check_mass_err
     procedure, public  :: Enter_tracer_LayerAdjustment
     procedure, public  :: Exit_tracer_LayerAdjustment
     procedure, public  :: tracer_snowcapping
     procedure, public  :: tracer_DivideSnowLayers
     procedure, public  :: tracer_CombineSnowLayers
     procedure, public  :: get_hist_size
     procedure, public  :: get_hist_info
     procedure, private :: ReadNamelist
     procedure, private :: create_betr_application
     procedure, public  :: HistRetrieveState
     procedure, public  :: HistRetrieveFlux
     procedure, public  :: set_restvar
     procedure, public  :: get_restartvar_size
     procedure, public  :: get_restartvar_info
     procedure, public  :: diagnoselnd2atm
     procedure, public  :: retrieve_biostates
     procedure, public  :: retrieve_biofluxes
     procedure, public  :: debug_info
  end type betr_type

  public :: create_betr_type
contains


  function create_betr_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_type), pointer :: create_betr_type
    class(betr_type), pointer :: betr

    allocate(betr)
    create_betr_type => betr

  end function create_betr_type

!-------------------------------------------------------------------------------
  subroutine Init(this, namelist_buffer, bounds, col, biophysforc, asoibgc, bstatus)

    ! FIXME(bja, 201604) need to remove waterstate, cnstate and
    ! ecophyscon from this routine.
    !USES
    use bshr_log_mod    , only : errMsg => shr_log_errMsg
    use BeTR_decompMod  , only : betr_bounds_type
    use betr_constants  , only : betr_namelist_buffer_size
    use TransportMod    , only : init_transportmod
    use TracerParamsMod , only : tracer_param_init
    use BetrStatusType  , only : betr_status_type
    use betr_columnType , only : betr_column_type
    implicit none
    !arguments
    class(betr_type)                         , intent(inout)        :: this
    character(len=betr_namelist_buffer_size) , intent(in)           :: namelist_buffer
    type(betr_bounds_type)                   , intent(in)           :: bounds
    type(betr_column_type)                   , intent(in)           :: col
    type(betr_biogeophys_input_type)         , intent(in)           :: biophysforc
    logical                                  , intent(out)          :: asoibgc
    type(betr_status_type)                   , intent(out)          :: bstatus

    !temporary variables
    character(len=1024) :: msg
    integer             :: lbj, ubj

    lbj = bounds%lbj;  ubj = bounds%ubj
    !read in top level control parameters
    call this%ReadNamelist(namelist_buffer, bstatus)
    if(bstatus%check_status())return

    !read in application specific parameters
    call this%create_betr_application(this%bgc_reaction, this%plant_soilbgc, this%reaction_method,&
        asoibgc, bstatus)
    if(bstatus%check_status())return

    this%active_soibgc = asoibgc
    if(asoibgc)then
      call this%plantNutkinetics%Init(bounds)
    endif

    call this%tracers%init_scalars()

    call this%bgc_reaction%Init_betrbgc(bounds, lbj, ubj, this%tracers, namelist_buffer, bstatus)
    if(bstatus%check_status())return

    call this%aereconds%Init(bounds)

    call init_transportmod()

    call this%tracerstates%Init(bounds, lbj, ubj, this%tracers)

    call this%tracerfluxes%Init(bounds,  lbj, ubj, this%tracers)

    call this%tracercoeffs%Init(bounds, lbj, ubj, this%tracers)

    call this%tracerboundaryconds%Init(bounds, this%tracers)

    !inside Init_plant_soilbgc, specific plant soil bgc coupler data type will be created
    call this%plant_soilbgc%Init_plant_soilbgc(bounds, lbj, ubj, namelist_buffer)

    !initialize state variable
    call this%bgc_reaction%initCold(bounds, col, this%tracers, biophysforc, this%tracerstates)

    !initialize boundary condition type
    call this%bgc_reaction%init_boundary_condition_type(bounds, this%tracers, this%tracerboundaryconds)

    !initialize the betr parameterization module
    call tracer_param_init(bounds)

   allocate(this%h2osoi_liq_copy(bounds%begc:bounds%endc, 1:nlevsoi));  this%h2osoi_liq_copy(:, :) = spval
   allocate(this%h2osoi_ice_copy(bounds%begc:bounds%endc, 1:nlevsoi));  this%h2osoi_ice_copy(:, :) = spval

  end subroutine Init

!-------------------------------------------------------------------------------
  subroutine ReadNamelist(this, namelist_buffer, bstatus)
    !
    ! !DESCRIPTION:
    ! read namelist for betr configuration
    ! !USES:
    use spmdMod        , only : masterproc, mpicom
    use betr_ctrl      , only : iulog => biulog
    use bshr_log_mod   , only : errMsg => shr_log_errMsg
    use betr_constants , only : stdout, betr_string_length_long, betr_namelist_buffer_size
    use BetrStatusType , only : betr_status_type
    use tracer_varcon  , only : advection_on, diffusion_on, reaction_on, ebullition_on, reaction_method
    implicit none
    ! !ARGUMENTS:
    class(betr_type)                         , intent(inout) :: this
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    type(betr_status_type)                   , intent(out) :: bstatus

    !
    ! !LOCAL VARIABLES:
    integer                                :: nml_error
    character(len=*), parameter            :: subname = 'ReadNamelist'
    character(len=betr_string_length_long) :: ioerror_msg

    !-----------------------------------------------------------------------

!    namelist / betr_parameters /                  &
!         reaction_method,                         &
!         advection_on, diffusion_on, reaction_on, &
!         ebullition_on, esm_on

    call bstatus%reset()
!    reaction_method = 'mock_run'
!    advection_on    = .true.
!    diffusion_on    = .true.
!    reaction_on     = .true.
!    ebullition_on   =.true.

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    !if ( index(trim(namelist_buffer),'betr_parameters')/=0 )then
    !   ioerror_msg=''
    !   read(namelist_buffer, nml=betr_parameters, iostat=nml_error, iomsg=ioerror_msg)
    !   if (nml_error /= 0) then
    !      call bstatus%set_msg(msg="ERROR reading betr_parameters namelist "//errmsg(filename, __LINE__),err=-1)
    !      return
    !   end if
    !end if

    !if (.not. esm_on) then
    !   write(stdout, *)
    !   write(stdout, *) '--------------------'
    !   write(stdout, *)
    !   write(stdout, *) ' betr bgc type :'
    !   write(stdout, *)
    !   write(stdout, *) ' betr_parameters namelist settings :'
    !   write(stdout, *)
    !   write(stdout, betr_parameters)
    !   write(stdout, *)
    !   write(stdout, *) '--------------------'
    !endif

    this%reaction_method = reaction_method
    this%advection_on    = advection_on
    this%diffusion_on    = diffusion_on
    this%reaction_on     = reaction_on
    this%ebullition_on   = ebullition_on
  end subroutine ReadNamelist

  !-------------------------------------------------------------------------------
  subroutine step_without_drainage(this, betr_time, bounds, col, pft, &
       num_soilc, filter_soilc, num_soilp, filter_soilp, &
       biophysforc, biogeo_flux, biogeo_state, betr_status)
    !
    ! !DESCRIPTION:
    ! run betr code one time step forward, without drainage calculation

    ! !USES:
    use tracerfluxType         , only : tracerflux_type
    use tracerstatetype        , only : tracerstate_type
    use tracercoeffType        , only : tracercoeff_type
    use TracerBoundaryCondType , only : TracerBoundaryCond_type
    use BetrTracerType         , only : betrtracer_type
    use BGCReactionsMod        , only : bgc_reaction_type
    use PlantSoilBGCMod        , only : plant_soilbgc_type
    use BeTR_aerocondType      , only : betr_aerecond_type
    use BetrBGCMod             , only : calc_ebullition
    use BetrBGCMod             , only : stage_tracer_transport
    use BetrBGCMod             , only : surface_tracer_hydropath_update
    use BetrBGCMod             , only : tracer_gws_transport
    use BeTR_TimeMod           , only : betr_time_type
    use BetrStatusType         , only : betr_status_type
    use betr_columnType        , only : betr_column_type
    use BeTR_PatchType         , only : betr_patch_type
    implicit none
    !
    ! !ARGUMENTS :
    class(betr_type)                 , intent(inout) :: this
    class(betr_time_type)            , intent(in)    :: betr_time
    type(bounds_type)                , intent(in)    :: bounds                     ! bounds
    type(betr_column_type)           , intent(in)    :: col
    type(betr_patch_type)            , intent(in)    :: pft
    integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
    integer                          , intent(in)    :: filter_soilc(:)            ! column filter_soilc
    integer                          , intent(in)    :: num_soilp
    integer                          , intent(in)    :: filter_soilp(:)            ! pft filter
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux
    type(betr_biogeo_state_type)     , intent(inout) :: biogeo_state
    type(betr_status_type)           , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    character(len=255) :: subname = 'run_betr_one_step_without_drainage'
    real(r8)           :: dtime
    real(r8)           :: Rfactor(bounds%begc:bounds%endc, bounds%lbj:bounds%ubj,1:this%tracers%ngwmobile_tracer_groups) !retardation factor
    integer            :: j
    integer            :: lbj, ubj


    call betr_status%reset()
    lbj = bounds%lbj; ubj = bounds%ubj

    dtime = betr_time%get_step_size()

    call stage_tracer_transport(betr_time, bounds, col, pft, num_soilc,&
         filter_soilc, num_soilp, filter_soilp, biophysforc,      &
         biogeo_state, biogeo_flux, this%aereconds, this%tracers, this%tracercoeffs, &
         this%tracerboundaryconds, this%tracerstates, this%tracerfluxes, this%bgc_reaction, &
         Rfactor, this%advection_on, betr_status)
    if(betr_status%check_status())return

    call surface_tracer_hydropath_update(betr_time, bounds, col, &
       num_soilc, filter_soilc,  biophysforc, this%tracers, this%tracerstates,    &
       this%tracercoeffs,  this%tracerfluxes, betr_status)
    if(betr_status%check_status())return

    if(this%active_soibgc)then
      !set up kinetic parameters that are passed in from the mother lsm. Mostly they
      !are plant-nutrient related parameters.
      call this%bgc_reaction%set_kinetics_par(1, ubj, this%nactpft, this%plantNutkinetics)
    endif
    if(this%reaction_on)                                       &
    call this%bgc_reaction%calc_bgc_reaction(bounds, col, lbj, ubj, &
         num_soilc,                                            &
         filter_soilc,                                         &
         num_soilp,                                            &
         filter_soilp,                                         &
         this%tracerboundaryconds%jtops_col,                   &
         dtime,                                                &
         this%tracers,                                         &
         this%tracercoeffs,                                    &
         biophysforc,                                          &
         this%tracerstates,                                    &
         this%tracerfluxes,                                    &
         this%tracerboundaryconds,                             &
         this%plant_soilbgc, biogeo_flux,  betr_status)
    if(betr_status%check_status())return

    call tracer_gws_transport(betr_time, bounds, col, pft, num_soilc, filter_soilc, &
      Rfactor, biophysforc, biogeo_flux, this%tracers, this%tracerboundaryconds  , &
      this%tracercoeffs,  this%tracerstates, this%tracerfluxes, this%bgc_reaction, &
      this%advection_on, this%diffusion_on, betr_status)
    if(betr_status%check_status())return

    call calc_ebullition(bounds, 1, ubj,                                                                  &
         this%tracerboundaryconds%jtops_col,                                                              &
         num_soilc,                                                                                       &
         filter_soilc,                                                                                    &
         biophysforc%forc_pbot_downscaled_col,                                                            &
         col%zi(bounds%begc:bounds%endc, 0:ubj),                                                          &
         col%dz(bounds%begc:bounds%endc, 1:ubj),                                                          &
         dtime,                                                                                           &
         biophysforc%fracice_col(bounds%begc:bounds%endc, 1:ubj),                                         &
         biogeo_state%zwts_col(bounds%begc:bounds%endc),                                                  &
         this%tracers,                                                                                    &
         this%tracercoeffs,                                                                               &
         this%tracerstates,                                                                               &
         this%tracerfluxes%tracer_flx_ebu_col(bounds%begc:bounds%endc, 1:this%tracers%nvolatile_tracers), &
         this%ebullition_on, betr_status)
    if(betr_status%check_status())return

    !update nutrient uptake fluxes
    call this%plant_soilbgc%plant_soilbgc_summary(bounds, lbj, ubj, pft, &
          num_soilc, filter_soilc,  dtime                              , &
          col%dz(bounds%begc:bounds%endc,1:ubj)                        , &
          this%tracers, this%tracerfluxes, biogeo_flux, betr_status)

  end subroutine step_without_drainage

  !--------------------------------------------------------------------------------
  subroutine debug_info(this, bounds, col, num_soilc, filter_soilc, header, betr_status)

  use betr_columnType        , only : betr_column_type

  implicit none
  ! !ARGUMENTS:
  class(betr_type)                     , intent(inout) :: this
  type(bounds_type)                    , intent(in)    :: bounds
  type(betr_column_type)               , intent(in)    :: col
  integer                              , intent(in)    :: num_soilc                   ! number of columns in column filter
  integer                              , intent(in)    :: filter_soilc(:)             ! column filter
  character(len=*)                     , intent(in)    :: header
  type(betr_status_type)               , intent(out)   :: betr_status

  call this%bgc_reaction%debug_info(bounds, num_soilc, filter_soilc, &
       col%dz(bounds%begc:bounds%endc, bounds%lbj:bounds%ubj), &
       this%tracers, this%tracerstates, header, betr_status)

  end subroutine debug_info

  !--------------------------------------------------------------------------------
  subroutine retrieve_biostates(this, bounds, lbj, ubj,  num_soilc, filter_soilc, jtops, &
    biogeo_state, betr_status)

  implicit none
  ! !ARGUMENTS:
  class(betr_type)                     , intent(inout) :: this
  type(bounds_type)                    , intent(in)  :: bounds                      ! bounds
  integer                              , intent(in) :: lbj, ubj
  integer                              , intent(in) :: jtops(bounds%begc:bounds%endc)
  integer                              , intent(in)    :: num_soilc                   ! number of columns in column filter
  integer                              , intent(in)    :: filter_soilc(:)             ! column filter
  type(betr_biogeo_state_type)         , intent(inout) :: biogeo_state
  type(betr_status_type)               , intent(out)   :: betr_status

  call this%bgc_reaction%retrieve_biostates(bounds, lbj, ubj, jtops, num_soilc, &
     filter_soilc, this%tracers, this%tracerstates, biogeo_state, betr_status)

  end subroutine retrieve_biostates

  !--------------------------------------------------------------------------------
  subroutine retrieve_biofluxes(this, num_soilc, filter_soilc, &
    biogeo_flux)

  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  implicit none
  ! !ARGUMENTS:
  class(betr_type)                     , intent(inout) :: this
  integer                              , intent(in)    :: num_soilc                   ! number of columns in column filter
  integer                              , intent(in)    :: filter_soilc(:)             ! column filter
  type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux

  call this%bgc_reaction%retrieve_biogeoflux(num_soilc, &
     filter_soilc, this%tracerfluxes, this%tracers, biogeo_flux)

  end subroutine retrieve_biofluxes

  !--------------------------------------------------------------------------------
  subroutine step_with_drainage(this, bounds, col, num_soilc, filter_soilc, jtops, &
    biogeo_flux, betr_status)
    !
    ! !DESCRIPTION:
    ! do tracer update due to drainage
    !
    ! !USES:
    use tracerfluxType  , only : tracerflux_type
    use tracerstatetype , only : tracerstate_type
    use tracercoeffType , only : tracercoeff_type
    use MathfuncMod     , only : safe_div
    use BetrBGCMod      , only : diagnose_gas_pressure
    use BetrStatusType  , only : betr_status_type
    use betr_columnType , only : betr_column_type
    implicit none
    ! !ARGUMENTS:
    class(betr_type)            , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    type(betr_column_type)      , intent(in)    :: col
    integer                     , intent(in)    :: num_soilc                          ! number of columns in column filter_soilc
    integer                     , intent(in)    :: filter_soilc(:)                    ! column filter_soilc
    integer                     , intent(in)    :: jtops(bounds%begc: )
    type(betr_biogeo_flux_type) , intent(in)    :: biogeo_flux
    type(betr_status_type)      , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    real(r8) :: aqucon
    integer  :: fc, c, j, k
    integer  :: lbj, ubj

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(jtops)  == (/bounds%endc/)), errMsg(filename,__LINE__), betr_status)
    if(betr_status%check_status())return
    associate(                                                                         & !
         ngwmobile_tracers     => this%tracers%ngwmobile_tracers         ,             & !
         groupid               => this%tracers%groupid                    ,            & !
         is_h2o                => this%tracers%is_h2o                       ,          & !
         is_advective          => this%tracers%is_advective                  ,         & !
         aqu2bulkcef_mobile    => this%tracercoeffs%aqu2bulkcef_mobile_col           , & !
         tracer_conc_mobile    => this%tracerstates%tracer_conc_mobile_col         ,   & !
         tracer_conc_grndwater => this%tracerstates%tracer_conc_grndwater_col   ,      & !
         dz                    => col%dz                                       ,       & !
         tracer_flx_drain      => this%tracerfluxes%tracer_flx_drain_col          ,    & !
         qflx_drain_vr         => biogeo_flux%qflx_drain_vr_col               ,        & ! Output  : [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/step) (to river +)
         qflx_totdrain         => biogeo_flux%qflx_totdrain_col                        & ! Output  : [real(r8) (: ,:) ]  (m H2o/step)
         )
      lbj = bounds%lbj; ubj = bounds%ubj

      do j = lbj, ubj
         do fc = 1, num_soilc
            c = filter_soilc(fc)
            if(j>=jtops(c))then
               do k = 1, ngwmobile_tracers
                  !obtain aqueous concentration
                  if(.not. is_advective(k))cycle
                  if(qflx_drain_vr(c,j) > 0._r8)then
                    aqucon = safe_div(tracer_conc_mobile(c,j,k),aqu2bulkcef_mobile(c,j,groupid(k)))
                  else
                    !when drainage is negative, tracer comes from groundwater
                    if(is_h2o(k))then
                      aqucon = tracer_conc_grndwater(c,k)
                    else
                      aqucon = 0._r8
                    endif
                  endif
                  !when drainage is negative, assume the flux is magically coming from other groundwater sources
                  tracer_flx_drain(c,k)     = tracer_flx_drain(c,k)  + aqucon * qflx_drain_vr(c,j)
                  tracer_conc_mobile(c,j,k) =  tracer_conc_mobile(c,j,k) - aqucon * qflx_drain_vr(c,j)/dz(c,j)
                  if(tracer_conc_mobile(c,j,k)<0._r8)then
                     tracer_flx_drain(c,k) = tracer_flx_drain(c,k)+tracer_conc_mobile(c,j,k)*dz(c,j)
                     tracer_conc_mobile(c,j,k)=0._r8
                  endif
               enddo
            endif
         enddo
      enddo
      if(betr_status%check_status())return
      !diagnose gas pressure
      call diagnose_gas_pressure(bounds, lbj, ubj, num_soilc, filter_soilc, &
           this%tracers, this%tracercoeffs, this%tracerstates, betr_status)

    end associate
  end subroutine step_with_drainage


  !--------------------------------------------------------------------------------
  subroutine diagnoselnd2atm(this, bounds, num_soilc, filter_soilc,  &
    biogeo_flux)
    !
    ! !DESCRIPTION:
    ! do tracer update due to drainage
    !
    ! !USES:
    use BetrStatusType  , only : betr_status_type
    use betr_columnType , only : betr_column_type
    implicit none
    ! !ARGUMENTS:
    class(betr_type)            , intent(inout) :: this
    type(bounds_type)           , intent(in)    :: bounds
    integer                     , intent(in)    :: num_soilc                          ! number of columns in column filter_soilc
    integer                     , intent(in)    :: filter_soilc(:)                    ! column filter_soilc
    type(betr_biogeo_flux_type) , intent(inout) :: biogeo_flux


    call this%bgc_reaction%retrieve_lnd2atm(bounds, num_soilc, filter_soilc, &
       this%tracerfluxes, this%tracers, biogeo_flux)

  end subroutine diagnoselnd2atm
  !--------------------------------------------------------------------------------
  subroutine calc_dew_sub_flux(this, betr_time, bounds, col, num_hydrologyc, &
      filter_soilc_hydrologyc,  biophysforc, betrtracer_vars, tracerflux_vars, tracerstate_vars)
    !
    ! DESCRIPTION:
    ! calculate water flux from dew formation, and sublimation
    ! !USES:
    use tracerfluxType      , only : tracerflux_type
    use tracerstatetype     , only : tracerstate_type
    use betr_varcon         , only : spval => bspval
    use tracer_varcon       , only : nlevtrc_soil  => betr_nlevtrc_soil
    use BeTR_landvarconType , only : landvarcon  => betr_landvarcon
    use BeTR_TimeMod        , only : betr_time_type
    use betr_columnType     , only : betr_column_type
    ! !ARGUMENTS:
    class(betr_type)                 , intent(inout) :: this
    class(betr_time_type)            , intent(in)    :: betr_time
    type(bounds_type)                , intent(in)    :: bounds
    type(betr_column_type)           , intent(in)    :: col
    integer                          , intent(in)    :: num_hydrologyc             ! number of column soil points in column filter_soilc
    integer                          , intent(in)    :: filter_soilc_hydrologyc(:) ! column filter_soilc for soil points
    type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
    class(betrtracer_type)           , intent(in)    :: betrtracer_vars            ! betr configuration information
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars            ! tracer flux
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars           ! tracer state variables data structure

    ! !LOCAL VARIABLES:
    real(r8) :: dtime
    real(r8) :: tot1, totz, tot0, tot2, tot3, tot4
    integer  :: fc, c, j, l, ll, ll2

    ! remove compiler warnings about unused dummy args
    if (this%advection_on) continue
    if (bounds%begc > 0)   continue

   associate(                                                                & !
         snl                =>    col%snl                                 ,  & ! Input:  [integer  (:)   ]  number of snow layers
         dz                 =>    col%dz                                  ,  & ! Input:  [real(r8) (:,:) ]  layer depth (m)
         h2osoi_ice         =>    biophysforc%h2osoi_ice_col          ,      & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)
         h2osoi_liq         =>    biophysforc%h2osoi_liq_col          ,      & ! Output: [real(r8) (:,:) ]  ice lens (kg/m2)
         frac_h2osfc        =>    biophysforc%frac_h2osfc_col         ,      & ! Input:  [real(r8) (:)   ]
         qflx_dew_grnd      =>    biophysforc%qflx_dew_grnd_col        ,     & ! Input:  [real(r8) (:)   ]  ground surface dew formation (mm H2O /s) [+]
         qflx_dew_snow      =>    biophysforc%qflx_dew_snow_col        ,     & ! Input:  [real(r8) (:)   ]  surface dew added to snow pack (mm H2O /s
         qflx_sub_snow_vol  =>    biophysforc%qflx_sub_snow_vol_col    ,     & ! Output: [real(r8) (:)   ]  sublimation rate from snow pack (mm H2O /s)
         qflx_snow2topsoi   =>    biophysforc%qflx_snow2topsoi_col     ,     & ! Output: [real(r8) (:)   ]
         qflx_h2osfc2topsoi =>    biophysforc%qflx_h2osfc2topsoi_col   ,     & !
         tracer_flx_dew_grnd=>    tracerflux_vars%tracer_flx_dew_grnd_col ,  & !
         tracer_flx_dew_snow=>    tracerflux_vars%tracer_flx_dew_snow_col ,  & !
         tracer_flx_sub_snow=>    tracerflux_vars%tracer_flx_sub_snow_col ,  & !
         tracer_conc_mobile =>    tracerstate_vars%tracer_conc_mobile_col ,  & !
         tracer_conc_frozen =>    tracerstate_vars%tracer_conc_frozen_col ,  & !
         is_h2o             =>    betrtracer_vars%is_h2o                  ,  & !
         tracernames        =>    betrtracer_vars%tracernames             ,  &
         frozenid           =>    betrtracer_vars%frozenid                ,  &
         ngwmobile_tracers  =>    betrtracer_vars%ngwmobile_tracers          &
         )


         !initialize the following fluxes to zero
      do fc = 1, num_hydrologyc
        c = filter_soilc_hydrologyc(fc)
        qflx_snow2topsoi(c) = 0._r8
        qflx_h2osfc2topsoi(c) = 0._r8
      enddo

      !do tracer update from dew and sublimation
      dtime = betr_time%get_step_size()

      do j = 1, ngwmobile_tracers
         !now only do water isotope tracer
         if(.not. is_h2o(j))cycle
         do fc = 1, num_hydrologyc
            c = filter_soilc_hydrologyc(fc)
            if(snl(c)+1>=1)then
               tracer_flx_dew_grnd(c, j) = (1._r8 - frac_h2osfc(c))*qflx_dew_grnd(c) * dtime
               tracer_flx_dew_snow(c, j) = (1._r8 - frac_h2osfc(c))*qflx_dew_snow(c) * dtime
               tracer_flx_sub_snow(c,j) = tracer_conc_frozen(c,1,frozenid(j)) * dz(c,1) * qflx_sub_snow_vol(c)
            else
               tracer_flx_dew_grnd(c, j) = 0._r8
               tracer_flx_dew_snow(c, j) = 0._r8
               tracer_flx_sub_snow(c, j) = 0._r8
            endif
         enddo
      enddo

      !apply those fluxes
      do j = 1, ngwmobile_tracers
         if(.not. is_h2o(j))cycle
         do fc = 1, num_hydrologyc
            c = filter_soilc_hydrologyc(fc)
            tracer_conc_mobile(c,1,j) = tracer_conc_mobile(c,1,j) + tracer_flx_dew_grnd(c, j)/dz(c,1)
            tracer_conc_frozen(c,1,frozenid(j)) = tracer_conc_frozen(c,1,frozenid(j)) + &
                 (tracer_flx_dew_snow(c, j)-tracer_flx_sub_snow(c,j))/dz(c,1)
         enddo
      enddo
    end associate

  end subroutine calc_dew_sub_flux


  !--------------------------------------------------------------------------------

  subroutine check_mass_err(this, col, betr_time, c, trcid, ubj, dz, betrtracer_vars,&
    tracerstate_vars, tracerflux_vars, bstatus)

  !DESCRIPTION
  !temporary mass balance check
  !USES
  use tracerfluxType  , only : tracerflux_type
  use tracerstatetype , only : tracerstate_type
  use BeTR_TimeMod    , only : betr_time_type
  use BetrStatusType  , only : betr_status_type
  use betr_columnType , only : betr_column_type
  implicit none
  !ARGUMENTS
  class(betr_type)        , intent(inout) :: this
  type(betr_column_type)  , intent(in)    :: col
  class(betr_time_type)   , intent(in)    :: betr_time
  integer                 , intent(in)    :: ubj
  integer                 , intent(in)    :: c, trcid
  real(r8)                , intent(in)    :: dz(1:ubj)
  class(betrtracer_type)  , intent(inout) :: betrtracer_vars  ! tracer info data structure
  type(tracerflux_type)   , intent(inout) :: tracerflux_vars  ! tracer flux
  class(tracerstate_type) , intent(inout) :: tracerstate_vars ! tracer state variables data structure
  type(betr_status_type)  , intent(out)   :: bstatus
  real(r8)                                :: totmass , err

  call bstatus%reset()
  ! remove compiler warnings about unused dummy args
  if (this%advection_on) continue

  associate(                                                                            &
       beg_tracer_molarmass      => tracerstate_vars%beg_tracer_molarmass_col         , &
       tracer_flx_netpro         => tracerflux_vars%tracer_flx_netpro_col             , &
       tracer_flx_netphyloss     => tracerflux_vars%tracer_flx_netphyloss_col         , &
       ntracers                  => betrtracer_vars%ntracers                          , &
       is_adsorb                 => betrtracer_vars%is_adsorb                         , &
       adsorbid                  => betrtracer_vars%adsorbid                          , &
       is_frozen                 => betrtracer_vars%is_frozen                         , &
       frozenid                  => betrtracer_vars%frozenid                            &
       )

  totmass=tracerstate_vars%int_mass_mobile_col(1,ubj,c,trcid,dz(1:ubj), bstatus)
  if(bstatus%check_status())return
  if(is_frozen(trcid))then
     totmass = totmass + &
          tracerstate_vars%int_mass_frozen_col(1,ubj,c,frozenid(trcid),dz(1:ubj),bstatus)
     if(bstatus%check_status())return
  endif
  call tracerflux_vars%flux_summary(col, betr_time, c, betrtracer_vars,bstatus)
  if(bstatus%check_status())return
  err=beg_tracer_molarmass(c,trcid)-totmass  &
       + tracer_flx_netpro(c,trcid)-tracer_flx_netphyloss(c,trcid)

  end associate
  end subroutine check_mass_err


   !------------------------------------------------------------------------
   subroutine diagnose_dtracer_freeze_thaw(this, bounds, num_nolakec, filter_nolakec, &
     biophysforc)
   !
   ! DESCRIPTION
   ! aqueous tracer partition based on freeze-thaw
   !
   ! USES
   !
   use tracerstatetype       , only : tracerstate_type
   use BeTRTracerType        , only : betrtracer_type
   use BeTR_landvarconType   , only : landvarcon => betr_landvarcon
   !
   ! Arguments
   implicit none
   class(betr_type)                 , intent(inout) :: this
   type(bounds_type)                , intent(in)    :: bounds
   integer                          , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
   integer                          , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points
   type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
   !temporary variables
   integer :: j, fc, c, l, k
   real(r8) :: thaw_frac, freeze_frac
   real(r8) :: dtracer

   if (bounds%begc > 0) continue

   associate(                                                        &
     frozenid         => this%tracers%frozenid                     , &
     is_frozen        => this%tracers%is_frozen                    , &
     ngwmobile_tracers=> this%tracers%ngwmobile_tracers            , &
     h2osoi_liq       => biophysforc%h2osoi_liq_col                , &
     h2osoi_ice       => biophysforc%h2osoi_ice_col                , &
     tracer_conc_mobile=> this%tracerstates%tracer_conc_mobile_col , &
     tracer_conc_frozen=> this%tracerstates%tracer_conc_frozen_col   &

   )
   ! the tracer concentration change between frozen and liquid
   !  phases due to freeze and thaw only occurs to non-volatile tracers
   ! and water tracers

   do j = 1, nlevsoi
     do fc = 1, num_nolakec
       c =  filter_nolakec(fc)
       do k = 1, ngwmobile_tracers
        !if it is a frozenable tracer, do it
        if(is_frozen(k))then
          if(h2osoi_liq(c,j) > this%h2osoi_liq_copy(c,j))then
            !thaw, solid to aqueous
            thaw_frac = min(1._r8-h2osoi_ice(c,j)/this%h2osoi_ice_copy(c,j),1._r8)
            dtracer = tracer_conc_frozen(c,j,frozenid(k)) * thaw_frac
            tracer_conc_frozen(c,j,frozenid(k)) = tracer_conc_frozen(c,j,frozenid(k)) - dtracer
            tracer_conc_mobile(c,j, k) = tracer_conc_mobile(c,j, k) + dtracer
          else
            !freeze, aqueous to solid
            freeze_frac = min(1._r8 - h2osoi_liq(c,j)/this%h2osoi_liq_copy(c,j),1._r8)
            dtracer = tracer_conc_mobile(c,j,k) * freeze_frac          !some modifier are needed to account for change in solubility
            tracer_conc_frozen(c,j,frozenid(k)) = tracer_conc_frozen(c,j,frozenid(k)) + dtracer
            tracer_conc_mobile(c,j, k) = tracer_conc_mobile(c,j, k) - dtracer
          endif
        endif
       enddo
     enddo
   enddo

   end associate
   end subroutine diagnose_dtracer_freeze_thaw

  !------------------------------------------------------------------------
  subroutine betr_lsm_flux_state_sendback(this, bounds,  num_soilc, filter_soilc )

   implicit none
   class(betr_type)  , intent(inout) :: this
   type(bounds_type) , intent(in)    :: bounds
   integer           , intent(in)    :: num_soilc                        ! number of column non-lake points in column filter
   integer           , intent(in)    :: filter_soilc(:)                  ! column filter for non-lake points

   ! remove compiler warnings about unused dummy args
   if (this%advection_on)      continue
   if (bounds%begc > 0)        continue
   if (num_soilc > 0)          continue
   if (size(filter_soilc) > 0) continue

!x   call this%bgc_reaction%betr_lsm_flux_state_sendback(bounds, &
!x       num_soilc, filter_soilc,                    &
!x       this%tracerstate_vars, this%tracerflux_vars,  this%betrtracer_vars)
  end subroutine betr_lsm_flux_state_sendback


   !------------------------------------------------------------------------
   subroutine pre_diagnose_soilcol_water_flux(this, bounds, num_nolakec, filter_nolakec, biophysforc)
   !
   ! DESCRIPTION
   ! pre diagnose advective water fluxes at different soil interfaces

   implicit none
   !ARGUMENTS
   class(betr_type)                 , intent(inout) :: this
   type(bounds_type)                , intent(in)    :: bounds
   integer                          , intent(in)    :: num_nolakec                        ! number of column non-lake points in column filter
   integer                          , intent(in)    :: filter_nolakec(:)                  ! column filter for non-lake points
   type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
   integer :: j, fc, c, l

   ! remove compiler warnings about unused dummy args
   if (bounds%begc > 0) continue

   do j = 1, nlevsoi
     do fc = 1, num_nolakec
       c =  filter_nolakec(fc)
       this%h2osoi_liq_copy(c,j) = biophysforc%h2osoi_liq_col(c,j)
       this%h2osoi_ice_copy(c,j) = biophysforc%h2osoi_ice_col(c,j)
     enddo
   enddo
   end subroutine pre_diagnose_soilcol_water_flux

   !------------------------------------------------------------------------
   subroutine diagnose_advect_water_flux(this, betr_time, &
        bounds, num_hydrologyc, filter_hydrologyc, &
        biophysforc, biogeo_flux)
   !
   ! DESCRIPTION
   ! diagnose advective water fluxes between different soil layers
   !
   use betr_varcon  , only : denh2o => bdenh2o
   use BeTR_TimeMod , only : betr_time_type
   implicit none
   !ARGUMENTS
   class(betr_type)                 , intent(inout) :: this
   class(betr_time_type)            , intent(in) :: betr_time
   type(bounds_type)                , intent(in)    :: bounds               ! bounds
   integer                          , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
   integer                          , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
   type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
   type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux

   !local variables
   integer :: j, fc, c
   real(r8):: dtime
   real(r8):: diff
   real(r8):: infl_tmp
   real(r8):: scal
   logical :: tf

   ! remove compiler warnings about unused dummy args
   if (bounds%begc > 0) continue
   associate(                                                          & !
     h2osoi_liq          =>    biophysforc%h2osoi_liq_col            , &
     h2osoi_ice          =>    biophysforc%h2osoi_ice_col            , &
     qflx_rootsoi        =>    biophysforc%qflx_rootsoi_col          , & ! Input  : [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/s) (+ = to atm)
     qflx_bot            =>    biophysforc%qflx_bot_col              , & ! Input : [real(r8)]
     qflx_adv            =>    biogeo_flux%qflx_adv_col              , & ! Output: [real(r8) (:,:) ]  water flux at interfaces       (m H2O/s) (- = to atm)
     qflx_gross_infl_soil=>    biogeo_flux%qflx_gross_infl_soil_col  , & ! Output: [real(r8) (:)] gross infiltration (mm H2O/s)
     qflx_infl           =>    biogeo_flux%qflx_infl_col             , & ! Output: [real(r8) (:)] infiltration, mm H2O/s
     qflx_gross_evap_soil=>    biogeo_flux%qflx_gross_evap_soil_col    & ! Output: [real(r8) (:)] gross evaporation (mm H2O/s)
   )

   ! get time step
   dtime = betr_time%get_step_size()
   !start from the bottom layer, because the water exchange between vadose zone soil and aquifer and plant root is known
   !the water flux at uppper surface can be inferred using the mass balance approach
   do j = nlevsoi, 1, -1
     do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       if(j==nlevsoi)then
         qflx_adv(c,j) = qflx_bot(c) * 1.e-3_r8                                 ! m/s
       else
         qflx_adv(c,j) = 1.e-3_r8 * (h2osoi_liq(c,j+1)-this%h2osoi_liq_copy(c,j+1))/dtime &
           + qflx_adv(c,j+1) + qflx_rootsoi(c,j+1)
       endif
     enddo
   enddo

   ! correct gross infiltration and gross evaporation
   ! (h2osoi_liq(c,1)-h2osoi_liq_copy(c,1))/dtime=qflx_infl-q_out-qflx_rootsoi
   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)

     !obtain the corrected infiltration
     qflx_infl(c) = (h2osoi_liq(c,1)-this%h2osoi_liq_copy(c,1))/dtime + (qflx_rootsoi(c,1)+qflx_adv(c,1))*1.e3_r8
     !the predicted net infiltration
     infl_tmp=qflx_gross_infl_soil(c)-qflx_gross_evap_soil(c)
     diff=qflx_infl(c)-infl_tmp
     if(abs(diff)>0._r8)then
       if(abs(infl_tmp)>0._r8)then
         scal = (1._r8+diff/infl_tmp)
         qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c) * scal
         qflx_gross_evap_soil(c) = qflx_gross_evap_soil(c) * scal
         if(qflx_gross_evap_soil(c)<0._r8)then
           !no negative evaporation allowed
           qflx_gross_infl_soil(c) = qflx_gross_infl_soil(c)-qflx_gross_evap_soil(c)
           qflx_gross_evap_soil(c) = 0._r8
         endif
         if(qflx_gross_infl_soil(c)<0._r8)then
           qflx_gross_evap_soil(c) = qflx_gross_evap_soil(c)-qflx_gross_infl_soil(c)
           qflx_gross_infl_soil(c) = 0._r8
         endif
       else
         if(diff>0._r8)then
           !the corrected infiltration > net infiltration
           qflx_gross_infl_soil(c)=qflx_gross_infl_soil(c) + diff
         else
           qflx_gross_evap_soil(c)=qflx_gross_evap_soil(c)-diff
         endif
       endif
     endif
     qflx_adv(c,0) = qflx_gross_infl_soil(c) *1.e-3_r8  !surface infiltration
   enddo
   end associate
   end subroutine diagnose_advect_water_flux

   !------------------------------------------------------------------------
   subroutine diagnose_drainage_water_flux(this, betr_time, &
        bounds, num_hydrologyc, filter_hydrologyc, &
        biophysforc, biogeo_flux)
   !
   ! DESCRIPTION
   ! diagnose advective water fluxes between different soil layers
   !
   !USES
   use betr_varcon  , only : denh2o => bdenh2o
   use BeTR_TimeMod , only : betr_time_type
   implicit none
   !ARGUMENTS
   class(betr_type)                 , intent(inout) :: this
   class(betr_time_type)            , intent(in)    :: betr_time
   type(bounds_type)                , intent(in)    :: bounds               ! bounds
   integer                          , intent(in)    :: num_hydrologyc       ! number of column soil points in column filter
   integer                          , intent(in)    :: filter_hydrologyc(:) ! column filter for soil points
   type(betr_biogeophys_input_type) , intent(in)    :: biophysforc
   type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux

   !local variables
   integer :: j, fc, c
   real(r8):: dtime

   ! remove compiler warnings about unused dummy args
   if (bounds%begc > 0) continue

   associate(                                                       & !
     h2osoi_liq         =>    biophysforc%h2osoi_liq_col          , & ! Output: [real(r8) (:,:) ] liquid water (kg/m2)
     h2osoi_ice         =>    biophysforc%h2osoi_ice_col          , & ! Output: [real(r8) (:,:) ] ice lens (kg/m2)
     qflx_drain_vr      =>    biogeo_flux%qflx_drain_vr_col       , & ! Output  : [real(r8) (:,:) ]  vegetation/soil water exchange (m H2O/step) (to river +)
     qflx_totdrain      =>    biogeo_flux%qflx_totdrain_col         & ! Output  : [real(r8) (:,:) ]  (m H2o/step)
   )

   ! get time step
   dtime = betr_time%get_step_size()
   !start from the bottom layer, because the water exchange between vadose zone soil and aquifer and plant root is known
   !the water flux at uppper surface can be inferred using the mass balance approach
   !also, the flux through transpiration has been already calcualted in the no drainage code, therefore no
   !the change in water content is solely due to subsurface draiange.
   do fc = 1, num_hydrologyc
     c = filter_hydrologyc(fc)
     qflx_totdrain(c) = 0._r8
   enddo
   do j = nlevsoi, 1, -1
     do fc = 1, num_hydrologyc
       c = filter_hydrologyc(fc)
       qflx_drain_vr(c,j) = this%h2osoi_liq_copy(c,j)-h2osoi_liq(c,j)   !kg/m2/step
       if(j==1)then
         !the following line is due to the ice check that is added in drainage calculation
         qflx_drain_vr(c,j) = qflx_drain_vr(c,j) + (this%h2osoi_ice_copy(c,j)-h2osoi_ice(c,j))
       endif
       !the following line will allow negative drainage
       qflx_drain_vr(c,j) =qflx_drain_vr(c,j) * 1.e-3_r8               !kg/m2/(kg/m3)/step = m/step
       qflx_totdrain(c) = qflx_totdrain(c) + qflx_drain_vr(c,j)
     enddo
   enddo

   end associate
   end subroutine diagnose_drainage_water_flux

   !------------------------------------------------------------------------

   subroutine tracer_DivideSnowLayers(this, bounds, col, num_snowc, filter_snowc, divide_matrix, bstatus)
   !
   ! DESCRIPTIONS
   ! divide tracers in snow
   !
   ! USES
   use BetrBGCMod, only : tracer_col_mapping_div
   use BetrBGCMod, only : tracer_copy_a2b_div
   use BetrStatusType      , only : betr_status_type
   use betr_columnType          , only : betr_column_type
   implicit none
   !ARGUMENTS
   class(betr_type)  , intent(inout) :: this
   type(bounds_type) , intent(in)    :: bounds               ! bounds
   type(betr_column_type), intent(in):: col
   integer           , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer           , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)          , intent(in)    :: divide_matrix(bounds%begc: , 1: , 1: )
   type(betr_status_type),intent(out):: bstatus

   !temporary variable
   real(r8) :: tracer_copy(bounds%begc:bounds%endc, 1:nlevsno)
   integer  :: jj, jj1

   call bstatus%reset()
   SHR_ASSERT_ALL((ubound(divide_matrix ,1)  == bounds%endc)  , errMsg(filename,__LINE__), bstatus)
   if(bstatus%check_status())return
   SHR_ASSERT_ALL((ubound(divide_matrix ,2)  == nlevsno)      , errMsg(filename,__LINE__), bstatus)
   if(bstatus%check_status())return
   SHR_ASSERT_ALL((ubound(divide_matrix ,3)  == nlevsno)      , errMsg(filename,__LINE__), bstatus)
   if(bstatus%check_status())return
   associate(                                                                           &
      tracer_conc_frozen_col        => this%tracerstates%tracer_conc_frozen_col ,       &
      tracer_conc_mobile_col        => this%tracerstates%tracer_conc_mobile_col ,       &
      ntracers                      => this%tracers%ntracers,                           &
      ngwmobile_tracers             => this%tracers%ngwmobile_tracers ,                 &
      is_frozen                     => this%tracers%is_frozen                           &
   )

   do jj = 1, ntracers
     if(jj<= ngwmobile_tracers)then
        !make a copy
        call tracer_copy_a2b_div(bounds, num_snowc, filter_snowc, col%snl,                     &
           tracer_conc_mobile_col(bounds%begc:bounds%endc, -nlevsno+1:0,jj),tracer_copy)
        !remapping
        call tracer_col_mapping_div(bounds, num_snowc, filter_snowc, col%snl, divide_matrix,   &
           tracer_copy, tracer_conc_mobile_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj),bstatus)
        if(bstatus%check_status())return
        if(is_frozen(jj))then
          !make a copy
          call tracer_copy_a2b_div(bounds, num_snowc, filter_snowc, col%snl,                   &
            tracer_conc_frozen_col(bounds%begc:bounds%endc, -nlevsno+1:0,jj),tracer_copy)
          !remapping
          call tracer_col_mapping_div(bounds, num_snowc, filter_snowc, col%snl, divide_matrix, &
             tracer_copy, tracer_conc_frozen_col(bounds%begc:bounds%endc, -nlevsno+1:0,jj),bstatus)
          if(bstatus%check_status())return
        endif
     else
       !x!currently, no solid tracer is allowed to mix with snow
       !x!copy
       !x!remapping
     endif
   enddo
   end associate
   end subroutine tracer_DivideSnowLayers
   !------------------------------------------------------------------------

   subroutine tracer_CombineSnowLayers(this, bounds, col, num_snowc, filter_snowc, combine_matrix, bstatus)
   !
   ! DESCRIPTIONS
   ! combine tracers in snow
   !
   !! USES
   use tracer_varcon , only : nlevsno => betr_nlevsno
   use BetrBGCMod    , only : tracer_col_mapping_div
   use BetrBGCMod    , only : tracer_col_mapping_comb
   use BetrBGCMod    , only : tracer_copy_a2b_comb
   use BetrStatusType      , only : betr_status_type
   use betr_columnType    , only : betr_column_type
   !!
   ! ARGUMENTS
   implicit none
   class(betr_type)  , intent(inout) :: this
   type(bounds_type) , intent(in)    :: bounds               ! bounds
   type(betr_column_type),intent(in) :: col
   integer           , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer           , intent(in)    :: filter_snowc(:) ! column filter for soil points
   real(r8)          , intent(in)    :: combine_matrix(bounds%begc: ,-nlevsno+1: ,-nlevsno+1: )
   type(betr_status_type),intent(out):: bstatus

   !temporary variables
   real(r8) :: tracer_copy(bounds%begc:bounds%endc,-nlevsno+1:1)
   integer  :: jj, jj1

   call bstatus%reset()
   SHR_ASSERT_ALL((ubound(combine_matrix,1) == bounds%endc)  , errMsg(filename,__LINE__),bstatus)
   if(bstatus%check_status())return
   SHR_ASSERT_ALL((ubound(combine_matrix,2) == 1)            , errMsg(filename,__LINE__),bstatus)
   if(bstatus%check_status())return
   SHR_ASSERT_ALL((ubound(combine_matrix,3) == 1)            , errMsg(filename,__LINE__),bstatus)
   if(bstatus%check_status())return

   associate(                                                                           &
      tracer_conc_frozen_col        => this%tracerstates%tracer_conc_frozen_col ,       &
      tracer_conc_mobile_col        => this%tracerstates%tracer_conc_mobile_col ,       &
      ntracers                      => this%tracers%ntracers,                           &
      ngwmobile_tracers             => this%tracers%ngwmobile_tracers ,                 &
      is_frozen                     => this%tracers%is_frozen                           &
   )


   do jj = 1, ntracers
     if(jj<= ngwmobile_tracers)then
        !make a copy
        call tracer_copy_a2b_comb(bounds, num_snowc, filter_snowc, col%snl,                      &
           tracer_conc_mobile_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj),tracer_copy)
        !remapping
        call tracer_col_mapping_comb(bounds, num_snowc, filter_snowc, col%snl, combine_matrix,   &
           tracer_copy, tracer_conc_mobile_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj),bstatus)
        if(bstatus%check_status())return
        if(is_frozen(jj))then
          !make a copy
          call tracer_copy_a2b_comb(bounds, num_snowc, filter_snowc, col%snl,                    &
            tracer_conc_frozen_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj),tracer_copy)
          !remapping
          call tracer_col_mapping_comb(bounds, num_snowc, filter_snowc, col%snl, combine_matrix, &
             tracer_copy, tracer_conc_frozen_col(bounds%begc:bounds%endc, -nlevsno+1:1,jj),bstatus)
          if(bstatus%check_status())return
        endif
     else
       !x!currently, no solid tracer is allowed to mix with snow
       !x!copy
       !x!remapping
     endif
   enddo

   end associate
   end subroutine tracer_CombineSnowLayers

   !------------------------------------------------------------------------
   subroutine Enter_tracer_LayerAdjustment(this, bounds, col, num_snowc, filter_snowc )
   !
   !! DESCRIPTION
   ! prepare tracer for entering mass adjustment in response to snow dynamics
   use betr_columnType          , only : betr_column_type
   implicit none
   !arguments
   class(betr_type)  , intent(inout) :: this
   type(bounds_type) , intent(in)    :: bounds               ! bounds
   type(betr_column_type),intent(in) :: col
   integer           , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer           , intent(in)    :: filter_snowc(:) ! column filter for soil points

   !temporary variables
   integer :: fc, c, j, jj

   ! remove compiler warnings about unused dummy args
   if (bounds%begc > 0) continue

   associate(                                                                           &
      tracer_conc_frozen_col        => this%tracerstates%tracer_conc_frozen_col ,       &
      tracer_conc_mobile_col        => this%tracerstates%tracer_conc_mobile_col ,       &
      ntracers                      => this%tracers%ntracers,                           &
      ngwmobile_tracers             => this%tracers%ngwmobile_tracers ,                 &
      is_frozen                     => this%tracers%is_frozen ,                         &
      snl                           => col%snl                                          &
   )
   do jj = 1, ntracers
     if(jj<= ngwmobile_tracers)then
       do j = -nlevsno + 1, 1
         do fc = 1, num_snowc
           c = filter_snowc(fc)
           if(j >= snl(c) + 1)then
             tracer_conc_mobile_col(c,j,jj) = tracer_conc_mobile_col(c,j,jj) /col%dz(c,j)
             if(is_frozen(jj))then
               tracer_conc_frozen_col(c,j,jj) = tracer_conc_frozen_col(c,j,jj)/col%dz(c,j)
             endif
           endif
         enddo
       enddo
     else
       !no passive tracer adjustment at this moment, can be used for aerosols or mineral particles
     endif
   enddo
   end associate
   end subroutine Enter_tracer_LayerAdjustment

   !------------------------------------------------------------------------
   subroutine Exit_tracer_LayerAdjustment(this, bounds, col, num_snowc, filter_snowc )
   !! DESCRIPTION
   ! prepare tracer for exit mass adjustment in response to snow dynamics
   use betr_columnType          , only : betr_column_type
   implicit none
   !ARGUMENTS
   class(betr_type)  , intent(inout) :: this
   type(bounds_type) , intent(in)    :: bounds               ! bounds
   type(betr_column_type),intent(in) :: col
   integer           , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer           , intent(in)    :: filter_snowc(:) ! column filter for soil points

   !temporary variables
   integer :: fc, c, j, jj

   ! remove compiler warnings about unused dummy args
   if (bounds%begc > 0) continue

   associate(                                                                           &
      tracer_conc_frozen_col        => this%tracerstates%tracer_conc_frozen_col ,       &
      tracer_conc_mobile_col        => this%tracerstates%tracer_conc_mobile_col ,       &
      ntracers                      => this%tracers%ntracers,                           &
      ngwmobile_tracers             => this%tracers%ngwmobile_tracers ,                 &
      is_frozen                     => this%tracers%is_frozen ,                         &
      snl                           => col%snl                                          &
   )
   do jj = 1, ntracers
     if(jj<= ngwmobile_tracers)then
       do j = -nlevsno + 1, 1
         do fc = 1, num_snowc
           c = filter_snowc(fc)
           if(j >= snl(c) + 1)then
             tracer_conc_mobile_col(c,j,jj) = tracer_conc_mobile_col(c,j,jj) *col%dz(c,j)
             if(is_frozen(jj))then
               tracer_conc_frozen_col(c,j,jj) = tracer_conc_frozen_col(c,j,jj)*col%dz(c,j)
             endif
           endif
         enddo
       enddo
     else
       !no passive tracer adjustment at this moment, can be used for aerosols or mineral particles
     endif
   enddo
   end associate
   end subroutine Exit_tracer_LayerAdjustment
   !------------------------------------------------------------------------

   subroutine tracer_snowcapping(this, bounds, num_snowc, filter_snowc, biophysforc)

   !do tracer update due to snow capping

   implicit none
   class(betr_type)                 , intent(inout) :: this
   type(bounds_type)                , intent(in)    :: bounds               ! bounds
   integer                          , intent(in)    :: num_snowc      ! number of column soil points in column filter
   integer                          , intent(in)    :: filter_snowc(:) ! column filter for soil points
   type(betr_biogeophys_input_type) , intent(in)    :: biophysforc


   ! remove compiler warnings about unused dummy args
   if (this%advection_on)                    continue
   if (bounds%begc > 0)                      continue
   if (num_snowc > 0)                        continue
   if (size(filter_snowc) > 0)               continue
   if (size(biophysforc%h2osoi_liq_col) > 0) continue

  end subroutine tracer_snowcapping
   !------------------------------------------------------------------------
  subroutine create_betr_application(this, bgc_reaction, plant_soilbgc, method, asoibgc, bstatus)
  !DESCRIPTION
  !create betr applications based on method
  !USES
  use ReactionsFactory    , only : create_betr_def_application
  use ApplicationsFactory , only : create_betr_usr_application
  use BetrStatusType      , only : betr_status_type
  implicit none
  !ARGUMENTS
  class(betr_type)          ,               intent(inout) :: this
  class(bgc_reaction_type)  ,  allocatable, intent(out)  :: bgc_reaction
  class(plant_soilbgc_type) ,  allocatable, intent(out)  :: plant_soilbgc
  character(len=*)          ,               intent(in)   :: method
  logical                   , intent(out) :: asoibgc
  type(betr_status_type)    , intent(out) :: bstatus

  !temporary variable
  logical :: yesno

  call bstatus%reset()

  asoibgc = .false.
  !if it is a default case, create it
  call create_betr_def_application(bgc_reaction, plant_soilbgc, method, yesno)

  if(.not. yesno)then
    call create_betr_usr_application(bgc_reaction, plant_soilbgc, method, bstatus)
  endif
  if(trim(method) =='eca_cnp')then
    asoibgc = .true.
  endif
  end subroutine create_betr_application

   !------------------------------------------------------------------------

  subroutine get_hist_size(this, num_state1d, num_state2d, num_flux1d, num_flux2d)
  !
  ! DESCRIPTION
  ! return number of variables for history output
  !
  ! USES
  use betr_ctrl, only : max_betr_hist_type
  use BetrStatusType      , only : betr_status_type
  use betr_constants , only : stdout
  implicit none
  !ARGUMENTS
  class(betr_type)  ,     intent(inout) :: this
  integer           ,     intent(out)   :: num_state1d
  integer           ,     intent(out)   :: num_state2d
  integer           ,     intent(out)   :: num_flux1d
  integer           ,     intent(out)   :: num_flux2d
  integer :: j

  num_state1d = this%tracerstates%num_hist1d

  num_state2d = this%tracerstates%num_hist2d

  num_flux1d  = this%tracerfluxes%num_hist1d

  num_flux2d  = this%tracerfluxes%num_hist2d

  end subroutine get_hist_size
   !------------------------------------------------------------------------

  subroutine get_hist_info(this, num_state1d, num_state2d, num_flux1d, num_flux2d, &
         state_hist1d_var, state_hist2d_var, &
         flux_hist1d_var, flux_hist2d_var)
  !
  ! DESCRIPTION
  ! return number of variables for history output
  !
  ! USES
  use betr_ctrl           , only : max_betr_hist_type
  use BetrStatusType      , only : betr_status_type
  use betr_constants      , only : stdout
  use BeTRHistVarType     , only : betr_hist_var_type, hist_var_copy
  implicit none
  !ARGUMENTS
  class(betr_type)  ,     intent(inout) :: this
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d
  type(betr_hist_var_type), intent(out) :: state_hist1d_var(1: )
  type(betr_hist_var_type), intent(out) :: state_hist2d_var(1: )
  type(betr_hist_var_type), intent(out) :: flux_hist1d_var(1: )
  type(betr_hist_var_type), intent(out) :: flux_hist2d_var(1: )
  integer :: j

  do j = 1, num_state1d
    call hist_var_copy(state_hist1d_var(j), this%tracerstates%hist1d_var(j))
  enddo

  do j = 1, num_state2d
    call hist_var_copy(state_hist2d_var(j),this%tracerstates%hist2d_var(j))
  enddo

  do j = 1, num_flux1d
    call hist_var_copy(flux_hist1d_var(j), this%tracerfluxes%hist1d_var(j))
  enddo

  do j = 1, num_flux2d
    call hist_var_copy(flux_hist2d_var(j), this%tracerfluxes%hist2d_var(j))
  enddo
  end subroutine get_hist_info

  !------------------------------------------------------------------------
  subroutine HistRetrieveState(this, bounds, lbj, ubj, &
     num_state1d, num_state2d, states_1d, states_2d)
  implicit none
  class(betr_type)  ,     intent(inout) :: this
  type(bounds_type) , intent(in)    :: bounds               ! bounds
  integer, intent(in) :: lbj, ubj
  integer           ,     intent(in)   :: num_state1d
  integer           ,     intent(in)   :: num_state2d
  real(r8)          , intent(out) :: states_1d(bounds%begc:bounds%endc,1:num_state1d)
  real(r8)          , intent(out) :: states_2d(bounds%begc:bounds%endc,lbj:ubj, 1:num_state2d)

  call this%tracerstates%retrieve_hist(bounds, lbj, ubj, states_2d, states_1d, this%tracers)

  end subroutine HistRetrieveState

  !------------------------------------------------------------------------
  subroutine HistRetrieveFlux(this, bounds, lbj, ubj, &
     num_flux1d, num_flux2d,  fluxes_1d, fluxes_2d)
  implicit none
  class(betr_type)  ,     intent(inout) :: this
  type(bounds_type) , intent(in)    :: bounds               ! bounds
  integer, intent(in) :: lbj, ubj
  integer           ,     intent(in)   :: num_flux1d
  integer           ,     intent(in)   :: num_flux2d
  real(r8)          , intent(out) :: fluxes_1d(bounds%begc:bounds%endc,1:num_flux1d)
  real(r8)          , intent(out) :: fluxes_2d(bounds%begc:bounds%endc,lbj:ubj, 1:num_flux2d)

  call this%tracerfluxes%retrieve_hist(bounds, lbj, ubj, fluxes_2d, fluxes_1d, this%tracers)

  end subroutine HistRetrieveFlux

  !------------------------------------------------------------------------
  subroutine get_restartvar_size(this, nrest_1d, nrest_2d)
  !
  !DESCRIPTION
  use betr_ctrl, only : max_betr_rest_type
  implicit none
  class(betr_type)  ,     intent(inout) :: this
  integer, intent(out) :: nrest_1d, nrest_2d

  call this%tracerstates%get_restartvars_size(nrest_1d, nrest_2d, this%tracers)

  end subroutine get_restartvar_size

  !------------------------------------------------------------------------
  subroutine get_restartvar_info(this, nrest_1d, nrest_2d,rest_varname_1d, &
   rest_varname_2d)
  !
  !DESCRIPTION
  use betr_ctrl, only : max_betr_rest_type
  implicit none
  class(betr_type)  ,     intent(inout) :: this
  integer, intent(in) :: nrest_1d, nrest_2d
  character(len=255), intent(out) :: rest_varname_1d(max_betr_rest_type)
  character(len=255), intent(out) :: rest_varname_2d(max_betr_rest_type)

  call this%tracerstates%get_restartvars_info(nrest_1d, nrest_2d,rest_varname_1d, &
   rest_varname_2d, this%tracers)

  end subroutine get_restartvar_info


  !------------------------------------------------------------------------
  subroutine set_restvar(this, bounds, lbj, ubj, nrest_1d, nrest_2d, states_1d, states_2d, flag)
  !
  !DESCRIPTIONS
  !set initial conditions based on restart file
  implicit none
  class(betr_type)  ,     intent(inout) :: this
  type(bounds_type)                , intent(in)    :: bounds               ! bounds
  integer, intent(in) :: nrest_1d, nrest_2d
  integer, intent(in) :: lbj, ubj
  real(r8), intent(inout) :: states_1d(bounds%begc:bounds%endc, 1:nrest_1d)
  real(r8), intent(inout) :: states_2d(bounds%begc:bounds%endc, lbj:ubj, 1:nrest_2d)
  character(len=*), intent(in) :: flag

  call this%tracerstates%restart(bounds, lbj, ubj, nrest_1d, nrest_2d, states_1d, states_2d, this%tracers, flag)
  end subroutine set_restvar
end module BetrType
