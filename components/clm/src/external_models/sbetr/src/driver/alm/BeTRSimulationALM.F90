module BeTRSimulationALM
  !
  ! !DESCRIPTION:
  !  API for using BeTR in ALM
  !
  ! !USES:
  !
#include "shr_assert.h"
  use abortutils          , only : endrun
  use clm_varctl          , only : iulog,use_cn
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use BeTRSimulation      , only : betr_simulation_type
  use decompMod           , only : bounds_type
  use BeTRSimulation      , only : betr_simulation_type
  use BeTR_TimeMod        , only : betr_time_type
  use tracer_varcon       , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
  use betr_decompMod      , only : betr_bounds_type
  use betr_varcon         , only : betr_maxpatch_pft
#if (defined SBETR)
  use PatchType      , only : patch_type
  use ColumnType     , only : column_type
  use LandunitType   , only : landunit_type
#else
  use VegetationType      , only : patch_type => vegetation_physical_properties_type
  use ColumnType          , only : column_type => column_physical_properties_type
  use LandunitType        , only : landunit_type => landunit_physical_properties_type
#endif
  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public, extends(betr_simulation_type) :: betr_simulation_alm_type
   contains
     procedure :: InitOnline                        => ALMInit
     procedure :: Init                              => ALMInitOffline
     procedure, public :: StepWithoutDrainage       => ALMStepWithoutDrainage
     procedure, public :: StepWithDrainage          => ALMStepWithDrainage
     procedure, public :: SetBiophysForcing         => ALMSetBiophysForcing
     !unique subroutines
     procedure, public :: CalcDewSubFlux            => ALMCalcDewSubFlux
     procedure, public :: SoilFluxStateRecv         => ALMBetrSoilFluxStateRecv
     procedure, public :: CalcSmpL                  => ALMCalcSmpL
     procedure, public :: PlantSoilBGCSend          => ALMBetrPlantSoilBGCSend
     procedure, public :: PlantSoilBGCRecv          => ALMBetrPlantSoilBGCRecv
     procedure, public :: DiagnoseLnd2atm           => ALMDiagnoseLnd2atm
     procedure, public :: set_active                => ALMset_active
     procedure, private:: set_transient_kinetics_par
  end type betr_simulation_alm_type

  public :: create_betr_simulation_alm

contains

!-------------------------------------------------------------------------------

  function create_betr_simulation_alm()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_simulation_alm_type), pointer :: create_betr_simulation_alm
    class(betr_simulation_alm_type), pointer :: simulation

    allocate(simulation)
    create_betr_simulation_alm => simulation

  end function create_betr_simulation_alm

!-------------------------------------------------------------------------------

  subroutine ALMInit(this, bounds, lun, col, pft, waterstate, namelist_buffer, masterproc)
    !DESCRIPTION
    !Initialize BeTR for ALM
    !
    !USES
    !data types from alm
    use pftvarcon       , only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
    use WaterStateType  , only : waterstate_type
    use landunit_varcon , only : istcrop, istice, istsoil
    use clm_varpar      , only : nlevsno, nlevsoi, nlevtrc_soil
    !betr types
    use betr_constants      , only : betr_filename_length
    use betr_constants      , only : betr_namelist_buffer_size
    use BeTR_pftvarconType  , only : betr_pftvarcon
    use BeTR_landvarconType , only : betr_landvarcon
    use BeTR_decompMod      , only : betr_bounds_type

    implicit none
    class(betr_simulation_alm_type)          , intent(inout) :: this
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(in) :: col
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

    ! now call the base simulation init to continue initialization
    if(present(masterproc))then
      call this%BeTRInit(bounds, lun, col, pft, waterstate, namelist_buffer, masterproc=masterproc)
    else
      call this%BeTRInit(bounds, lun, col, pft, waterstate, namelist_buffer)
    endif
  end subroutine ALMInit
!-------------------------------------------------------------------------------

  subroutine ALMInitOffline(this, bounds, lun, col, pft, waterstate, namelist_buffer, base_filename)
    !DESCRIPTION
    !Initialize BeTR for ALM
    !
    !USES
    !data types from alm
    use pftvarcon       , only : noveg, nc4_grass, nc3_arctic_grass, nc3_nonarctic_grass
    use WaterStateType  , only : waterstate_type
    use landunit_varcon , only : istcrop, istice, istsoil
    use clm_varpar      , only : nlevsno, nlevsoi, nlevtrc_soil
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
    type(bounds_type)                        , intent(in)    :: bounds
    type(landunit_type)                      , intent(in) :: lun
    type(column_type)                        , intent(in) :: col
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
    use clm_varpar        , only : nlevsno, nlevsoi, nlevtrc_soil
    use tracer_varcon     , only : betr_nlevsoi, betr_nlevsno, betr_nlevtrc_soil
    implicit none
    ! !ARGUMENTS :
    class(betr_simulation_alm_type) , intent(inout) :: this
    type(bounds_type)               , intent(in)    :: bounds ! bounds
    type(column_type)               , intent(in)    :: col ! column type
    type(patch_type)                , intent(in)    :: pft
    !TEMPORARY VARIABLES
    type(betr_bounds_type)     :: betr_bounds
    integer :: c, c_l, begc_l, endc_l

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
      call this%biogeo_flux(c)%reset(value_column=0._r8, active_soibgc=this%active_soibgc)

!!
!--------------
!  debug
      call this%biogeo_state(c)%reset(value_column=0._r8, active_soibgc=this%active_soibgc)

      call this%betr(c)%retrieve_biostates(betr_bounds,      &
         1, betr_nlevsoi, this%num_soilc, this%filter_soilc, this%jtops, this%biogeo_state(c),this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif

        call this%biogeo_state(c)%summary(betr_bounds, 1, betr_nlevtrc_soil,this%betr_col(c)%dz(begc_l:endc_l,1:betr_nlevtrc_soil),&
          this%betr_col(c)%zi(begc_l:endc_l,1:betr_nlevtrc_soil), this%active_soibgc)

      if(.false.)then
        call this%betr(c)%debug_info(betr_bounds, this%betr_col(c), this%num_soilc, this%filter_soilc, 'bef w/o drain',&
             this%bstatus(c))
        this%betr(c)%tracers%debug=.true.
      endif
!--------
      call this%betr(c)%step_without_drainage(this%betr_time, betr_bounds, this%betr_col(c), &
         this%betr_pft(c), this%num_soilc, this%filter_soilc, this%num_soilp, this%filter_soilp, &
         this%biophys_forc(c), this%biogeo_flux(c), this%biogeo_state(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif
!--------------
!  debug
      call this%biogeo_state(c)%reset(value_column=0._r8, active_soibgc=this%active_soibgc)

      call this%betr(c)%retrieve_biostates(betr_bounds,      &
         1, betr_nlevsoi, this%num_soilc, this%filter_soilc, this%jtops, this%biogeo_state(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif

      call this%biogeo_state(c)%summary(betr_bounds, 1, betr_nlevtrc_soil,this%betr_col(c)%dz(begc_l:endc_l,1:betr_nlevtrc_soil), &
         this%betr_col(c)%zi(begc_l:endc_l,1:betr_nlevtrc_soil), this%active_soibgc)

      if(.false.)call this%betr(c)%debug_info(betr_bounds, this%betr_col(c), this%num_soilc, this%filter_soilc, 'aft w/o drain',&
           this%bstatus(c))
!--------
    enddo
    if(this%bsimstatus%check_status()) &
      call endrun(msg=this%bsimstatus%print_msg())
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
  subroutine ALMStepWithDrainage(this, bounds,  col)
   !
   !DESCRIPTION
   !interface for using diagnose land fluxes to atm and river copmonents
   !
   !USES
    use MathfuncMod   , only : safe_div
    use lnd2atmType    , only : lnd2atm_type
    use clm_varpar     , only : nlevsno, nlevsoi, nlevtrc_soil
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
      if(.false.)call this%betr(c)%debug_info(betr_bounds, this%betr_col(c), this%num_soilc, this%filter_soilc, 'bfdrain',&
           this%bstatus(c))
      call this%betr(c)%step_with_drainage(betr_bounds,      &
         this%betr_col(c),this%num_soilc, this%filter_soilc, this%jtops, &
         this%biogeo_flux(c), this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif

      call this%biogeo_state(c)%reset(value_column=0._r8, active_soibgc=this%active_soibgc)

      call this%betr(c)%retrieve_biostates(betr_bounds,      &
         1, betr_nlevsoi, this%num_soilc, this%filter_soilc, this%jtops, this%biogeo_state(c),this%bstatus(c))

      if(this%bstatus(c)%check_status())then
        call this%bsimstatus%setcol(c)
        call this%bsimstatus%set_msg(this%bstatus(c)%print_msg(),this%bstatus(c)%print_err())
        exit
      endif

      call this%biogeo_state(c)%summary(betr_bounds, 1, betr_nlevtrc_soil,this%betr_col(c)%dz(begc_l:endc_l,1:betr_nlevtrc_soil), &
          this%betr_col(c)%zi(begc_l:endc_l,1:betr_nlevtrc_soil),this%active_soibgc)

! debug
      if(.false.)call this%betr(c)%debug_info(betr_bounds, this%betr_col(c), this%num_soilc, this%filter_soilc, 'afdrain', &
           this%bstatus(c))
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
    use clm_varpar     , only : nlevsno, nlevsoi, nlevtrc_soil
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
           this%num_soilc, this%filter_soilc, this%biogeo_flux(c))
    enddo


    if(trim(reaction_method)=='doc_dic')then

      c_l = 1
      do c = bounds%begc, bounds%endc
        if(.not. this%active_col(c))cycle
        qflx_rofliq_qsur_doc_col(c)=this%biogeo_flux(c)%qflx_rofliq_qsur_doc_col(c_l)
        qflx_rofliq_qsur_dic_col(c)=this%biogeo_flux(c)%qflx_rofliq_qsur_dic_col(c_l)
        qflx_rofliq_qsub_doc_col(c)=this%biogeo_flux(c)%qflx_rofliq_qsub_doc_col(c_l)
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
  subroutine ALMBetrPlantSoilBGCSend(this, bounds, col, pft, num_soilc,  filter_soilc, cnstate_vars, &
    carbonflux_vars,  c13_cflx_vars, c14_cflx_vars, nitrogenflux_vars, phosphorusflux_vars, &
    PlantMicKinetics_vars)

  !read in biogeochemical fluxes from alm for soil bgc modeling
  !these are C, N and P fluxes from root input, surface litter input
  !atmospheric deposition, fire (negative), and fertilization
  !Because of possible harvest activity that is
  !related to dynamic land use, input profiles are computed in alm.
  !
  use CNCarbonFluxType, only : carbonflux_type
  use CNNitrogenFluxType, only : nitrogenflux_type
  use PhosphorusFluxType, only : phosphorusflux_type
  use CNStateType, only : cnstate_type
  use clm_varpar, only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use PlantMicKineticsMod, only : PlantMicKinetics_type
  use mathfuncMod, only : apvb
  use tracer_varcon, only : use_c13_betr, use_c14_betr
  implicit none
  class(betr_simulation_alm_type), intent(inout)  :: this
  type(bounds_type) , intent(in)  :: bounds
  type(column_type) , intent(in)  :: col ! column type
  type(patch_type)  , intent(in)  :: pft ! pft type
  integer           , intent(in)  :: num_soilc
  integer           , intent(in)  :: filter_soilc(:)
  type(cnstate_type), intent(in)  :: cnstate_vars
  type(carbonflux_type), intent(in):: carbonflux_vars
  type(carbonflux_type), intent(in):: c13_cflx_vars
  type(carbonflux_type), intent(in):: c14_cflx_vars
  type(nitrogenflux_type), intent(in):: nitrogenflux_vars
  type(phosphorusflux_type), intent(in):: phosphorusflux_vars
  type(PlantMicKinetics_type), intent(in) :: PlantMicKinetics_vars

  !temporary variables
  type(betr_bounds_type) :: betr_bounds
  integer :: c, fc, j, c_l
  ! remove compiler warnings
  if (this%num_soilc > 0)     continue
  if (bounds%begc > 0)        continue
  if (num_soilc > 0)          continue
  if (size(filter_soilc) > 0) continue

  associate(                                           &
    ndep_prof     => cnstate_vars%ndep_prof_col     ,  &
    pdep_prof     => cnstate_vars%pdep_prof_col     ,  &
    nfixation_prof=> cnstate_vars%nfixation_prof_col,  &
    frac_loss_lit_to_fire_col=> cnstate_vars%frac_loss_lit_to_fire_col, &
    frac_loss_cwd_to_fire_col=> cnstate_vars%frac_loss_cwd_to_fire_col  &
  )
  call this%BeTRSetBounds(betr_bounds)

  !set kinetic parameters
  call this%set_transient_kinetics_par(betr_bounds, col, pft, num_soilc, filter_soilc, PlantMicKinetics_vars)

  !set biophysical forcing
  c_l = 1
  do fc = 1, num_soilc
    c = filter_soilc(fc)
    call this%biophys_forc(c)%reset(value_column=0._r8)
    this%biophys_forc(c)%isoilorder(c_l) = 1                 !this needs update
    this%biophys_forc(c)%frac_loss_lit_to_fire_col(c_l) =frac_loss_lit_to_fire_col(c)
    this%biophys_forc(c)%frac_loss_cwd_to_fire_col(c_l) =frac_loss_cwd_to_fire_col(c)

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
  do j = betr_bounds%lbj, betr_bounds%ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      !!------------------------------------------------------------------------
      !carbon input
      !metabolic carbon
      call apvb(this%biophys_forc(c)%c12flx%cflx_input_litr_met_vr_col(1,j), & !
         (/carbonflux_vars%phenology_c_to_litr_met_c_col(c,j)       , & !phenology
         carbonflux_vars%dwt_frootc_to_litr_met_c_col(c,j)          , & !dynamic land cover
         carbonflux_vars%gap_mortality_c_to_litr_met_c_col(c,j)     , & !gap mortality
         carbonflux_vars%harvest_c_to_litr_met_c_col(c,j)           , & !harvest
         carbonflux_vars%m_c_to_litr_met_fire_col(c,j)/))              ! fire mortality

      !cellulose carbon
      call apvb(this%biophys_forc(c)%c12flx%cflx_input_litr_cel_vr_col(1,j)   , &
         (/carbonflux_vars%phenology_c_to_litr_cel_c_col(c,j)          , &  !phenology
         carbonflux_vars%dwt_frootc_to_litr_cel_c_col(c,j)             , &  !dynamic land cover
         carbonflux_vars%gap_mortality_c_to_litr_cel_c_col(c,j)        , & !gap mortality
         carbonflux_vars%harvest_c_to_litr_cel_c_col(c,j)              , & !harvest
         carbonflux_vars%m_c_to_litr_cel_fire_col(c,j)/))              ! fire mortality

      !lignin carbon
      call apvb(this%biophys_forc(c)%c12flx%cflx_input_litr_lig_vr_col(1,j) , &
         (/carbonflux_vars%phenology_c_to_litr_lig_c_col(c,j)        , & !phenology
         carbonflux_vars%dwt_frootc_to_litr_lig_c_col(c,j)           , & !dynamic land cover
         carbonflux_vars%gap_mortality_c_to_litr_lig_c_col(c,j)      , & !gap mortality
         carbonflux_vars%harvest_c_to_litr_lig_c_col(c,j)            , & !harvest
         carbonflux_vars%m_c_to_litr_lig_fire_col(c,j)/))                ! fire mortality

      !cwd carbon
      call apvb(this%biophys_forc(c)%c12flx%cflx_input_litr_cwd_vr_col(1,j) , &
        (/carbonflux_vars%dwt_livecrootc_to_cwdc_col(c,j)            , &
        carbonflux_vars%dwt_deadcrootc_to_cwdc_col(c,j)              , &
        carbonflux_vars%gap_mortality_c_to_cwdc_col(c,j)             , &
        carbonflux_vars%harvest_c_to_cwdc_col(c,j)                   , &
        carbonflux_vars%fire_mortality_c_to_cwdc_col(c,j)/))

      !fire carbon loss
      call apvb(this%biophys_forc(c)%c12flx%cflx_output_litr_met_vr_col(1,j), &
         carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_met_lit))

      call apvb(this%biophys_forc(c)%c12flx%cflx_output_litr_cel_vr_col(1,j), &
         carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_cel_lit))

      call apvb(this%biophys_forc(c)%c12flx%cflx_output_litr_lig_vr_col(1,j), &
         carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_lig_lit))

      call apvb(this%biophys_forc(c)%c12flx%cflx_output_litr_cwd_vr_col(1,j), &
         carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_cwd))
      !!------------------------------------------------------------------------
      if(use_c13_betr)then
        !metabolic carbon
        call apvb(this%biophys_forc(c)%c13flx%cflx_input_litr_met_vr_col(1,j), & !
           (/carbonflux_vars%phenology_c_to_litr_met_c_col(c,j)       , & !phenology
           carbonflux_vars%dwt_frootc_to_litr_met_c_col(c,j)          , & !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_met_c_col(c,j)     , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_met_c_col(c,j)           , & !harvest
           carbonflux_vars%m_c_to_litr_met_fire_col(c,j)/))              ! fire mortality

        !cellulose carbon
        call apvb(this%biophys_forc(c)%c13flx%cflx_input_litr_cel_vr_col(1,j)   , &
           (/carbonflux_vars%phenology_c_to_litr_cel_c_col(c,j)          , &  !phenology
           carbonflux_vars%dwt_frootc_to_litr_cel_c_col(c,j)             , &  !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_cel_c_col(c,j)        , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_cel_c_col(c,j)              , & !harvest
           carbonflux_vars%m_c_to_litr_cel_fire_col(c,j)/))              ! fire mortality

        !lignin carbon
        call apvb(this%biophys_forc(c)%c13flx%cflx_input_litr_lig_vr_col(1,j) , &
           (/carbonflux_vars%phenology_c_to_litr_lig_c_col(c,j)        , & !phenology
           carbonflux_vars%dwt_frootc_to_litr_lig_c_col(c,j)           , & !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_lig_c_col(c,j)      , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_lig_c_col(c,j)            , & !harvest
           carbonflux_vars%m_c_to_litr_lig_fire_col(c,j)/))                ! fire mortality

        !cwd carbon
        call apvb(this%biophys_forc(c)%c13flx%cflx_input_litr_cwd_vr_col(1,j) , &
          (/carbonflux_vars%dwt_livecrootc_to_cwdc_col(c,j)            , &
          carbonflux_vars%dwt_deadcrootc_to_cwdc_col(c,j)              , &
          carbonflux_vars%gap_mortality_c_to_cwdc_col(c,j)             , &
          carbonflux_vars%harvest_c_to_cwdc_col(c,j)                   , &
          carbonflux_vars%fire_mortality_c_to_cwdc_col(c,j)/))

        !fire carbon loss
        call apvb(this%biophys_forc(c)%c13flx%cflx_output_litr_met_vr_col(1,j), &
           carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_met_lit))

        call apvb(this%biophys_forc(c)%c13flx%cflx_output_litr_cel_vr_col(1,j), &
           carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_cel_lit))

        call apvb(this%biophys_forc(c)%c13flx%cflx_output_litr_lig_vr_col(1,j), &
           carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_lig_lit))

        call apvb(this%biophys_forc(c)%c13flx%cflx_output_litr_cwd_vr_col(1,j), &
           carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_cwd))
      endif
      if(use_c14_betr)then
        !metabolic carbon
        call apvb(this%biophys_forc(c)%c14flx%cflx_input_litr_met_vr_col(1,j), & !
           (/carbonflux_vars%phenology_c_to_litr_met_c_col(c,j)       , & !phenology
           carbonflux_vars%dwt_frootc_to_litr_met_c_col(c,j)          , & !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_met_c_col(c,j)     , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_met_c_col(c,j)           , & !harvest
           carbonflux_vars%m_c_to_litr_met_fire_col(c,j)/))              ! fire mortality

        !cellulose carbon
        call apvb(this%biophys_forc(c)%c14flx%cflx_input_litr_cel_vr_col(1,j)   , &
           (/carbonflux_vars%phenology_c_to_litr_cel_c_col(c,j)          , &  !phenology
           carbonflux_vars%dwt_frootc_to_litr_cel_c_col(c,j)             , &  !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_cel_c_col(c,j)        , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_cel_c_col(c,j)              , & !harvest
           carbonflux_vars%m_c_to_litr_cel_fire_col(c,j)/))              ! fire mortality

        !lignin carbon
        call apvb(this%biophys_forc(c)%c14flx%cflx_input_litr_lig_vr_col(1,j) , &
           (/carbonflux_vars%phenology_c_to_litr_lig_c_col(c,j)        , & !phenology
           carbonflux_vars%dwt_frootc_to_litr_lig_c_col(c,j)           , & !dynamic land cover
           carbonflux_vars%gap_mortality_c_to_litr_lig_c_col(c,j)      , & !gap mortality
           carbonflux_vars%harvest_c_to_litr_lig_c_col(c,j)            , & !harvest
           carbonflux_vars%m_c_to_litr_lig_fire_col(c,j)/))                ! fire mortality

        !cwd carbon
        call apvb(this%biophys_forc(c)%c14flx%cflx_input_litr_cwd_vr_col(1,j) , &
          (/carbonflux_vars%dwt_livecrootc_to_cwdc_col(c,j)            , &
          carbonflux_vars%dwt_deadcrootc_to_cwdc_col(c,j)              , &
          carbonflux_vars%gap_mortality_c_to_cwdc_col(c,j)             , &
          carbonflux_vars%harvest_c_to_cwdc_col(c,j)                   , &
          carbonflux_vars%fire_mortality_c_to_cwdc_col(c,j)/))

        !fire carbon loss
        call apvb(this%biophys_forc(c)%c14flx%cflx_output_litr_met_vr_col(1,j), &
           carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_met_lit))

        call apvb(this%biophys_forc(c)%c14flx%cflx_output_litr_cel_vr_col(1,j), &
           carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_cel_lit))

        call apvb(this%biophys_forc(c)%c14flx%cflx_output_litr_lig_vr_col(1,j), &
           carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_lig_lit))

        call apvb(this%biophys_forc(c)%c14flx%cflx_output_litr_cwd_vr_col(1,j), &
           carbonflux_vars%m_decomp_cpools_to_fire_vr_col(c,j,i_cwd))
      endif

      !nitrogen input
      !metabolic nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_input_litr_met_vr_col(c_l,j) , &
         (/nitrogenflux_vars%phenology_n_to_litr_met_n_col(c,j)      , & !phenology
         nitrogenflux_vars%dwt_frootn_to_litr_met_n_col(c,j)         , & !dynamic land cover
         nitrogenflux_vars%gap_mortality_n_to_litr_met_n_col(c,j)    , & !gap mortality
         nitrogenflux_vars%harvest_n_to_litr_met_n_col(c,j)          , & !harvest
         nitrogenflux_vars%m_n_to_litr_met_fire_col(c,j)/))              ! fire mortality

      !cellulose nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_input_litr_cel_vr_col(c_l,j), &
         (/nitrogenflux_vars%phenology_n_to_litr_cel_n_col(c,j)     , & !phenology
         nitrogenflux_vars%dwt_frootn_to_litr_cel_n_col(c,j)        , & !dynamic land cover
         nitrogenflux_vars%gap_mortality_n_to_litr_cel_n_col(c,j)   , & !gap mortality
         nitrogenflux_vars%harvest_n_to_litr_cel_n_col(c,j)         , & !harvest
         nitrogenflux_vars%m_n_to_litr_cel_fire_col(c,j)/))             ! fire mortality

      !lignin nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_input_litr_lig_vr_col(c_l,j) , &
         (/nitrogenflux_vars%phenology_n_to_litr_lig_n_col(c,j)      , &  !phenology
         nitrogenflux_vars%dwt_frootn_to_litr_lig_n_col(c,j)         , &   !dynamic land cover
         nitrogenflux_vars%gap_mortality_n_to_litr_lig_n_col(c,j)    , & !gap mortality
         nitrogenflux_vars%harvest_n_to_litr_lig_n_col(c,j)          , & !harvest
         nitrogenflux_vars%m_n_to_litr_lig_fire_col(c,j)/))              ! fire mortality

      !cwd nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_input_litr_cwd_vr_col(c_l,j) , &
        (/nitrogenflux_vars%dwt_livecrootn_to_cwdn_col(c,j)          , &
        nitrogenflux_vars%dwt_deadcrootn_to_cwdn_col(c,j)            , &
        nitrogenflux_vars%gap_mortality_n_to_cwdn_col(c,j)           , &
        nitrogenflux_vars%harvest_n_to_cwdn_col(c,j)                 , &
        nitrogenflux_vars%fire_mortality_n_to_cwdn_col(c,j)/))

      !fire nitrogen loss
      call apvb(this%biophys_forc(c)%n14flx%nflx_output_litr_met_vr_col(c_l,j) , &
         nitrogenflux_vars%m_decomp_npools_to_fire_vr_col(c,j,i_met_lit))

      call apvb(this%biophys_forc(c)%n14flx%nflx_output_litr_cel_vr_col(c_l,j) , &
         nitrogenflux_vars%m_decomp_npools_to_fire_vr_col(c,j,i_cel_lit))

      call apvb(this%biophys_forc(c)%n14flx%nflx_output_litr_lig_vr_col(c_l,j) , &
         nitrogenflux_vars%m_decomp_npools_to_fire_vr_col(c,j,i_lig_lit))

      call apvb(this%biophys_forc(c)%n14flx%nflx_output_litr_cwd_vr_col(c_l,j) , &
         nitrogenflux_vars%m_decomp_npools_to_fire_vr_col(c,j,i_cwd))

      !phosphorus input
      !metabolic phosphorus
      call apvb(this%biophys_forc(c)%p31flx%pflx_input_litr_met_vr_col(c_l,j) , &
         (/phosphorusflux_vars%phenology_p_to_litr_met_p_col(c,j)    , & !phenology
         phosphorusflux_vars%dwt_frootp_to_litr_met_p_col(c,j)       , & !dynamic land cover
         phosphorusflux_vars%gap_mortality_p_to_litr_met_p_col(c,j)  , & !gap mortality
         phosphorusflux_vars%harvest_p_to_litr_met_p_col(c,j)        , & !harvest
         phosphorusflux_vars%m_p_to_litr_met_fire_col(c,j)/))            ! fire mortality

      !cellulose phosphorus
      call apvb(this%biophys_forc(c)%p31flx%pflx_input_litr_cel_vr_col(c_l,j) , &
         (/phosphorusflux_vars%phenology_p_to_litr_cel_p_col(c,j)    , & !phenology
         phosphorusflux_vars%dwt_frootp_to_litr_cel_p_col(c,j)       , & !dynamic land cover
         phosphorusflux_vars%gap_mortality_p_to_litr_cel_p_col(c,j)  , & !gap mortality
         phosphorusflux_vars%harvest_p_to_litr_cel_p_col(c,j)        , & !harvest
         phosphorusflux_vars%m_p_to_litr_cel_fire_col(c,j)/))            ! fire mortality

      !lignin phosphorus
      call apvb(this%biophys_forc(c)%p31flx%pflx_input_litr_lig_vr_col(c_l,j) , &
         (/phosphorusflux_vars%phenology_p_to_litr_lig_p_col(c,j)    , & !phenology
         phosphorusflux_vars%dwt_frootp_to_litr_lig_p_col(c,j)       , & !dynamic land cover
         phosphorusflux_vars%gap_mortality_p_to_litr_lig_p_col(c,j)  , & !gap mortality
         phosphorusflux_vars%harvest_p_to_litr_lig_p_col(c,j)        , & !harvest
         phosphorusflux_vars%m_p_to_litr_lig_fire_col(c,j)/))            ! fire mortality

      !cwd phosphorus
      call apvb(this%biophys_forc(c)%p31flx%pflx_input_litr_cwd_vr_col(c_l,j) , &
        (/phosphorusflux_vars%dwt_livecrootp_to_cwdp_col(c,j) , &
        phosphorusflux_vars%dwt_deadcrootp_to_cwdp_col(c,j)   , &
        phosphorusflux_vars%gap_mortality_p_to_cwdp_col(c,j)  , &
        phosphorusflux_vars%harvest_p_to_cwdp_col(c,j)        , &
        phosphorusflux_vars%fire_mortality_p_to_cwdp_col(c,j)/))

      !fire phosphorus loss
      call apvb(this%biophys_forc(c)%p31flx%pflx_output_litr_met_vr_col(c_l,j) , &
         phosphorusflux_vars%m_decomp_ppools_to_fire_vr_col(c,j,i_met_lit))

      call apvb(this%biophys_forc(c)%p31flx%pflx_output_litr_cel_vr_col(c_l,j) , &
         phosphorusflux_vars%m_decomp_ppools_to_fire_vr_col(c,j,i_cel_lit))

      call apvb(this%biophys_forc(c)%p31flx%pflx_output_litr_lig_vr_col(c_l,j) , &
         phosphorusflux_vars%m_decomp_ppools_to_fire_vr_col(c,j,i_lig_lit))

      call apvb(this%biophys_forc(c)%p31flx%pflx_output_litr_cwd_vr_col(c_l,j) , &
         phosphorusflux_vars%m_decomp_ppools_to_fire_vr_col(c,j,i_cwd))

      !mineral nitrogen
      call apvb(this%biophys_forc(c)%n14flx%nflx_minn_input_nh4_vr_col(c_l,j) , &
         (/nitrogenflux_vars%ndep_to_sminn_col(c)                    , &
         nitrogenflux_vars%fert_to_sminn_col(c)/),  ndep_prof(c,j))

      !the following could be commented out if a fixation model is done in betr
      call apvb(this%biophys_forc(c)%n14flx%nflx_minn_nh4_fix_nomic_vr_col(c_l,j) , &
         (/nitrogenflux_vars%nfix_to_sminn_col(c)                        , &
         nitrogenflux_vars%soyfixn_to_sminn_col(c)/), nfixation_prof(c,j))

      !mineral phosphorus, the deposition is assumed to be of primary form
      call apvb(this%biophys_forc(c)%p31flx%pflx_minp_input_po4_vr_col(c_l,j) , &
         (/phosphorusflux_vars%fert_p_to_sminp_col(c)/),   pdep_prof(c,j))

      call apvb(this%biophys_forc(c)%p31flx%pflx_minp_weathering_po4_vr_col(c_l,j), &
         phosphorusflux_vars%primp_to_labilep_vr_col(c,j))

    enddo
  enddo

  end associate
  !pull in all state variables and update tracers
  end subroutine ALMBetrPlantSoilBGCSend

  !------------------------------------------------------------------------
  subroutine ALMBetrPlantSoilBGCRecv(this, bounds, col, pft, num_soilc,  filter_soilc,&
   c12state_vars, c12flux_vars, c13state_vars, c13flux_vars, c14state_vars, c14flux_vars, &
   n14state_vars, n14flux_vars, p31state_vars, p31flux_vars)
  !this returns the flux back to ALM after doing soil BGC
  !this specifically returns plant nutrient yield
  use CNCarbonFluxType    , only : carbonflux_type
  use CNCarbonStateType   , only : carbonstate_type
  use CNNitrogenFluxType  , only : nitrogenflux_type
  use CNNitrogenStateType , only : nitrogenstate_type
  !!! add phosphorus
  use PhosphorusFluxType  , only : phosphorusflux_type
  use PhosphorusStateType , only : phosphorusstate_type
  use tracer_varcon       , only : use_c13_betr, use_c14_betr
  use pftvarcon           , only : noveg
  use MathfuncMod         , only : safe_div
  implicit none
  class(betr_simulation_alm_type), intent(inout)  :: this
  type(bounds_type) , intent(in)  :: bounds
  type(patch_type)            , intent(in) :: pft
  type(column_type)           , intent(in)    :: col ! column type
  integer           , intent(in)  :: num_soilc
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
  integer :: c, fc, p, pi, c_l

    !TEMPORARY VARIABLES
  type(betr_bounds_type)     :: betr_bounds
  integer :: begc_l, endc_l

  !summarize the fluxes and state variables
  c_l = 1
  call this%BeTRSetBounds(betr_bounds)
  begc_l = betr_bounds%begc; endc_l=betr_bounds%endc;

  if(this%active_soibgc)then
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      call this%betr(c)%retrieve_biofluxes(this%num_soilc, this%filter_soilc, this%biogeo_flux(c))
      call this%biogeo_flux(c)%summary(betr_bounds, 1, betr_nlevtrc_soil, this%betr_col(c)%dz(begc_l:endc_l,1:betr_nlevtrc_soil))
    enddo


    !retrieve plant nutrient uptake from biogeo_flux
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      pi = 0
      do p = col%pfti(c), col%pftf(c)
        if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
          pi = pi + 1
          n14flux_vars%smin_nh4_to_plant_patch(p) = this%biogeo_flux(c)%n14flux_vars%smin_nh4_to_plant_patch(pi)
          n14flux_vars%smin_no3_to_plant_patch(p) = this%biogeo_flux(c)%n14flux_vars%smin_no3_to_plant_patch(pi)
          p31flux_vars%sminp_to_plant_patch(p)  = this%biogeo_flux(c)%p31flux_vars%sminp_to_plant_patch(pi)
          !compute relative n return.
          n14state_vars%pnup_pfrootc_patch(p) = safe_div(n14flux_vars%smin_nh4_to_plant_patch(p)&
               + n14flux_vars%smin_no3_to_plant_patch(p), &
            c12state_vars%frootc_patch(p))
        else
          n14flux_vars%smin_nh4_to_plant_patch(p) = 0._r8
          n14flux_vars%smin_no3_to_plant_patch(p) = 0._r8
          p31flux_vars%sminp_to_plant_patch(p) = 0._r8
        endif
      enddo

      !recollect soil respirations, fire and hydraulic loss
      c12flux_vars%hr_col(c) = this%biogeo_flux(c)%c12flux_vars%hr_col(c_l)
      c12flux_vars%fire_decomp_closs_col(c) = this%biogeo_flux(c)%c12flux_vars%fire_decomp_closs_col(c_l)
      c12flux_vars%som_c_leached_col(c) = &
        this%biogeo_flux(c)%c12flux_vars%som_c_leached_col(c_l) + &
        this%biogeo_flux(c)%c12flux_vars%som_c_runoff_col(c_l) + &
        this%biogeo_flux(c)%c12flux_vars%som_c_qdrain_col(c_l)
      !the following is for consistency with the ALM definitation, which computes

      !som_c_leached_col as a numerical roundoff
      c12flux_vars%som_c_leached_col(c)=-c12flux_vars%som_c_leached_col(c)
      if(use_c13_betr)then
        c13flux_vars%hr_col(c) = this%biogeo_flux(c)%c13flux_vars%hr_col(c_l)
        c13flux_vars%fire_decomp_closs_col(c) = this%biogeo_flux(c)%c13flux_vars%fire_decomp_closs_col(c_l)
      endif
      if(use_c14_betr)then
        c14flux_vars%hr_col(c) = this%biogeo_flux(c)%c14flux_vars%hr_col(c_l)
        c14flux_vars%fire_decomp_closs_col(c) = this%biogeo_flux(c)%c14flux_vars%fire_decomp_closs_col(c_l)
      endif

      !recollect  nitrifications, nitrifier-N2O loss, denitrifications
      n14flux_vars%f_nit_col(c) = this%biogeo_flux(c)%n14flux_vars%f_nit_col(c_l)
      n14flux_vars%f_denit_col(c)= this%biogeo_flux(c)%n14flux_vars%f_denit_col(c_l)
      n14flux_vars%denit_col(c)= n14flux_vars%f_denit_col(c)
      n14flux_vars%f_n2o_nit_col(c)=this%biogeo_flux(c)%n14flux_vars%f_n2o_nit_col(c_l)

      !hydraulic loss
      n14flux_vars%smin_no3_leached_col(c)= &
          this%biogeo_flux(c)%n14flux_vars%smin_no3_leached_col(c_l) + &
          this%biogeo_flux(c)%n14flux_vars%smin_no3_qdrain_col(c_l)
      n14flux_vars%som_n_leached_col(c) = &
          this%biogeo_flux(c)%n14flux_vars%som_n_leached_col(c_l) + &
          this%biogeo_flux(c)%n14flux_vars%som_n_runoff_col(c_l) + &
          this%biogeo_flux(c)%n14flux_vars%som_n_qdrain_col(c_l)

      if(.false.)then
        write(*,*)'betr smin_no3_leach=',this%biogeo_flux(c)%n14flux_vars%smin_no3_leached_col(c_l)*1800._r8
        write(*,*)'betr smin_no3_qdrai=',this%biogeo_flux(c)%n14flux_vars%smin_no3_qdrain_col(c_l)*1800._r8
        write(*,*)'betr som_n_leach   =',this%biogeo_flux(c)%n14flux_vars%som_n_leached_col(c_l)*1800._r8
        write(*,*)'betr som_n_runoff  =',this%biogeo_flux(c)%n14flux_vars%som_n_runoff_col(c_l)*1800._r8
        write(*,*)'betr som_n_drain   =',this%biogeo_flux(c)%n14flux_vars%som_n_qdrain_col(c_l)*1800._r8
      endif

      !the following is for consistency with the ALM definitation, which computes
      !som_n_leached_col as a numerical roundoff
      n14flux_vars%som_n_leached_col(c) = - n14flux_vars%som_n_leached_col(c)
      !fire loss
      n14flux_vars%smin_no3_runoff_col(c)=this%biogeo_flux(c)%n14flux_vars%smin_no3_runoff_col(c_l)
      n14flux_vars%fire_decomp_nloss_col(c) = this%biogeo_flux(c)%n14flux_vars%fire_decomp_nloss_col(c_l)

      !no nh4 volatilization and runoff/leaching loss at this moment

      !recollect mineral phosphorus loss
      !Remark: now hydraulic mineral p loss lumps all three fluxes, Jinyun Tang
      p31flux_vars%sminp_leached_col(c) = &
         this%biogeo_flux(c)%p31flux_vars%sminp_leached_col(c_l) + &
         this%biogeo_flux(c)%p31flux_vars%sminp_runoff_col(c_l) + &
         this%biogeo_flux(c)%p31flux_vars%sminp_qdrain_col(c_l)

      p31flux_vars%supplement_to_sminp_col(c) = this%biogeo_flux(c)%p31flux_vars%supplement_to_sminp_col(c_l)
      p31flux_vars%secondp_to_occlp_col(c) = this%biogeo_flux(c)%p31flux_vars%secondp_to_occlp_col(c_l)
      p31flux_vars%fire_decomp_ploss_col(c) = this%biogeo_flux(c)%p31flux_vars%fire_decomp_ploss_col(c_l)

      p31flux_vars%som_p_leached_col(c) = &
          this%biogeo_flux(c)%p31flux_vars%som_p_leached_col(c_l) + &
          this%biogeo_flux(c)%p31flux_vars%som_p_runoff_col(c_l) + &
          this%biogeo_flux(c)%p31flux_vars%som_p_qdrain_col(c_l)

      !the following is for consistency with the ALM definitation, which computes
      !som_p_leached_col as a numerical roundoff
      p31flux_vars%som_p_leached_col(c) = -p31flux_vars%som_p_leached_col(c)

      !recollect soil organic carbon, soil organic nitrogen, and soil organic phosphorus
      c12state_vars%cwdc_col(c) = this%biogeo_state(c)%c12state_vars%cwdc_col(c_l)
      c12state_vars%totlitc_col(c) = this%biogeo_state(c)%c12state_vars%totlitc_col(c_l)
      c12state_vars%totsomc_col(c) = this%biogeo_state(c)%c12state_vars%totsomc_col(c_l)
      c12state_vars%totlitc_1m_col(c) = this%biogeo_state(c)%c12state_vars%totlitc_1m_col(c_l)
      c12state_vars%totsomc_1m_col(c) = this%biogeo_state(c)%c12state_vars%totsomc_1m_col(c_l)

      if(use_c13_betr)then
        c13state_vars%cwdc_col(c) = this%biogeo_state(c)%c13state_vars%cwdc_col(c_l)
        c13state_vars%totlitc_col(c) = this%biogeo_state(c)%c13state_vars%totlitc_col(c_l)
        c13state_vars%totsomc_col(c) = this%biogeo_state(c)%c13state_vars%totsomc_col(c_l)
        c13state_vars%totlitc_1m_col(c) = this%biogeo_state(c)%c13state_vars%totlitc_1m_col(c_l)
        c13state_vars%totsomc_1m_col(c) = this%biogeo_state(c)%c13state_vars%totsomc_1m_col(c_l)
      endif
      if(use_c14_betr)then
        c14state_vars%cwdc_col(c) = this%biogeo_state(c)%c14state_vars%cwdc_col(c_l)
        c14state_vars%totlitc_col(c) = this%biogeo_state(c)%c14state_vars%totlitc_col(c_l)
        c14state_vars%totsomc_col(c) = this%biogeo_state(c)%c14state_vars%totsomc_col(c_l)
        c13state_vars%totlitc_1m_col(c) = this%biogeo_state(c)%c13state_vars%totlitc_1m_col(c_l)
        c13state_vars%totsomc_1m_col(c) = this%biogeo_state(c)%c13state_vars%totsomc_1m_col(c_l)
      endif
      n14state_vars%cwdn_col(c) = this%biogeo_state(c)%n14state_vars%cwdn_col(c_l)
      n14state_vars%totlitn_col(c) = this%biogeo_state(c)%n14state_vars%totlitn_col(c_l)
      n14state_vars%totsomn_col(c) = this%biogeo_state(c)%n14state_vars%totsomn_col(c_l)
      n14state_vars%totlitn_1m_col(c) = this%biogeo_state(c)%n14state_vars%totlitn_1m_col(c_l)
      n14state_vars%totsomn_1m_col(c) = this%biogeo_state(c)%n14state_vars%totsomn_1m_col(c_l)

      p31state_vars%cwdp_col(c) = this%biogeo_state(c)%p31state_vars%cwdp_col(c_l)
      p31state_vars%totlitp_col(c) = this%biogeo_state(c)%p31state_vars%totlitp_col(c_l)
      p31state_vars%totsomp_col(c) = this%biogeo_state(c)%p31state_vars%totsomp_col(c_l)
      p31state_vars%totlitp_1m_col(c) = this%biogeo_state(c)%p31state_vars%totlitp_1m_col(c_l)
      p31state_vars%totsomp_1m_col(c) = this%biogeo_state(c)%p31state_vars%totsomp_1m_col(c_l)

      !recollect inorganic nitrogen (smin_nh4, smin_no3), and inorganic phosphorus (disolvable and protected)
      n14state_vars%sminn_col(c) = this%biogeo_state(c)%n14state_vars%sminn_col(c_l)
      p31state_vars%sminp_col(c) = this%biogeo_state(c)%p31state_vars%sminp_col(c_l)
      p31state_vars%occlp_col(c) = this%biogeo_state(c)%p31state_vars%occlp_col(c_l)
    enddo
  endif
  end subroutine ALMBetrPlantSoilBGCRecv
  !------------------------------------------------------------------------

  subroutine ALMCalcDewSubFlux(this,  &
       bounds, col, num_hydrologyc, filter_soilc_hydrologyc)
   !DESCRIPTION
    ! Calculate tracer flux from dew or/and sublimation
    !External interface called by ALM

    use WaterfluxType   , only : waterflux_type
    use WaterstateType  , only : waterstate_type
    use clm_varcon      , only : denh2o,spval
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
         betr_bounds, this%betr_col(c), this%num_soilc, this%filter_soilc, &
        this%biophys_forc(c), this%betr(c)%tracers, this%betr(c)%tracerfluxes, this%betr(c)%tracerstates)
    enddo
  end subroutine ALMCalcDewSubFlux

  !------------------------------------------------------------------------
  subroutine ALMBetrSoilFluxStateRecv(this, num_soilc, filter_soilc)
  !this should be expanded and called after tracer update with drainage
  implicit none
  ! !ARGUMENTS:
  class(betr_simulation_alm_type), intent(inout) :: this
  integer           , intent(in)  :: num_soilc
  integer           , intent(in)  :: filter_soilc(:)

  integer :: fc, c
  type(betr_bounds_type)     :: betr_bounds

  call this%BeTRSetBounds(betr_bounds)


!  do fc = 1, num_soilc
!    c = filter_soilc(fc)
!    if(.not. this%active_col(c))cycle
!    call this%betr(c)%bgc_reaction%lsm_betr_flux_state_receive(betr_bounds, &
!       this%num_soilc, this%filter_soilc,                                   &
!       this%betr(c)%tracerstates, this%betr(c)%tracerfluxes,  this%betr(c)%tracers)
!  enddo

  end subroutine ALMBetrSoilFluxStateRecv

  !------------------------------------------------------------------------
  subroutine ALMCalcSmpL(this, bounds, lbj, ubj, numf, filter, t_soisno, &
     soilstate_vars, waterstate_vars, soil_water_retention_curve)
  !DESCRIPTION
  ! calculate soil suction potential
  !
  !USES
  use SoilStateType              , only : soilstate_type
  use WaterStateType             , only : waterstate_type
  use SoilWaterRetentionCurveMod , only : soil_water_retention_curve_type
  use clm_varcon                 , only : grav,hfus,tfrz
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

  SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, ubj/)),errMsg(mod_filename,__LINE__))

  ! remove compiler warnings
  if (this%num_soilc > 0) continue

  associate(                                                     & !
    h2osoi_vol        =>    waterstate_vars%h2osoi_vol_col     , & ! Input:  [real(r8) (:,:) ]  volumetric soil moisture
    smp_l             =>    waterstate_vars%smp_l_col          , & ! Output: [real(r8) (:,:) ]  soil suction (mm)
    bsw               =>    soilstate_vars%bsw_col             , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"
    watsat            =>    soilstate_vars%watsat_col          , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
    sucsat            =>    soilstate_vars%sucsat_col            & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm)
  )

  do j = lbj, ubj
    do fc = 1, numf
      c = filter(fc)
      if(.not. this%active_col(c))cycle
      if(j>=1)then
        if(t_soisno(c,j)<tfrz)then
          smp_l(c,j)= hfus*(tfrz-t_soisno(c,j))/(grav*t_soisno(c,j)) * 1000._r8  !(mm)
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
  subroutine ALMSetBiophysForcing(this, bounds, col, pft, carbonflux_vars, waterstate_vars, &
    waterflux_vars, temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
    chemstate_vars, soilstate_vars, cnstate_vars, carbonstate_vars)
  !DESCRIPTION
  !pass in biogeophysical variables for running betr
  !USES
  use SoilStateType     , only : soilstate_type
  use WaterStateType    , only : Waterstate_Type
  use TemperatureType   , only : temperature_type
  use ChemStateType     , only : chemstate_type
  use WaterfluxType     , only : waterflux_type
  use atm2lndType       , only : atm2lnd_type
  use SoilHydrologyType , only : soilhydrology_type
  use CNStateType       , only : cnstate_type
  use CNCarbonFluxType  , only : carbonflux_type
  use CanopyStateType   , only : canopystate_type
  use clm_varpar        , only : nlevsno, nlevsoi
  use tracer_varcon     , only : reaction_method
  use CNCarbonStateType , only : carbonstate_type
  use tracer_varcon     , only : catomw
  implicit none
  !ARGUMENTS
  class(betr_simulation_alm_type) , intent(inout)        :: this
  type(bounds_type)               , intent(in)           :: bounds
  type(patch_type)            , intent(in) :: pft
  type(column_type)           , intent(in)    :: col ! column type
  type(cnstate_type)          , optional, intent(in) :: cnstate_vars
  type(carbonflux_type)       , optional, intent(in) :: carbonflux_vars
  type(Waterstate_Type)       , optional, intent(in) :: Waterstate_vars
  type(waterflux_type)        , optional, intent(in) :: waterflux_vars
  type(temperature_type)      , optional, intent(in) :: temperature_vars
  type(soilhydrology_type)    , optional, intent(in) :: soilhydrology_vars
  type(atm2lnd_type)          , optional, intent(in) :: atm2lnd_vars
  type(canopystate_type)      , optional, intent(in) :: canopystate_vars
  type(chemstate_type)        , optional, intent(in) :: chemstate_vars
  type(soilstate_type)        , optional, intent(in) :: soilstate_vars
  type(carbonstate_type)      , optional, intent(in) :: carbonstate_vars


  integer :: p, pi, c, j, c_l
  integer :: npft_loc

  call this%BeTRSetBiophysForcing(bounds, col, pft, 1, nlevsoi, carbonflux_vars, waterstate_vars, &
      waterflux_vars, temperature_vars, soilhydrology_vars, atm2lnd_vars, canopystate_vars, &
      chemstate_vars, soilstate_vars)


  !the following will be ALM specific
  !big leaf model
  !set profiles autotrohpic respiration
  do c = bounds%begc, bounds%endc
    npft_loc = ubound(carbonflux_vars%rr_patch,1)-lbound(carbonflux_vars%rr_patch,1)+1
    if(npft_loc /= col%npfts(c) .and. col%pfti(c) /= lbound(carbonflux_vars%rr_patch,1)) then
      do pi = 1, betr_maxpatch_pft
        this%biophys_forc(c)%rr_patch(pi,1:nlevsoi) = 0._r8
      enddo
    else
      if(use_cn)then
        do pi = 1, betr_maxpatch_pft
          if (pi <= col%npfts(c)) then
            p = col%pfti(c) + pi - 1
            if (pft%active(p)) then
              this%biophys_forc(c)%rr_patch(pi,1:nlevsoi) = carbonflux_vars%rr_patch(p) !* root_prof(p,1:nlevsoi)
            endif
          endif
        enddo
      else
        do pi = 1, betr_maxpatch_pft
          this%biophys_forc(c)%rr_patch(pi,1:nlevsoi) = 0._r8
        enddo
      endif
    endif
  enddo
  !dvgm
  if(trim(reaction_method)=='doc_dic')then
     c_l=1
     do j = 1, nlevsoi
        do c = bounds%begc, bounds%endc
          if(col%active(c))then
             !for simplicity, atomic weight of carbon is set to 12._r8 g/mol
             this%biophys_forc(c)%dic_prod_vr_col(c_l,j) = (carbonflux_vars%hr_vr_col(c,j) + &
                cnstate_vars%nfixation_prof_col(c,j)*carbonflux_vars%rr_col(c))/catomw
             this%biophys_forc(c)%doc_prod_vr_col(c_l,j) = (carbonstate_vars%decomp_cpools_vr_col(c,j,6) - &
                carbonstate_vars%decomp_som2c_vr_col(c,j))/this%betr_time%delta_time/catomw
          endif
        enddo
      enddo
  endif
  end subroutine ALMSetBiophysForcing
  !------------------------------------------------------------------------
  subroutine set_transient_kinetics_par(this, betr_bounds, col, pft, num_soilc, filter_soilc, PlantMicKinetics_vars)
  !DESCRIPTION
  !set kinetic parameters for column c
  use PlantMicKineticsMod, only : PlantMicKinetics_type
  use tracer_varcon      , only : reaction_method
  use pftvarcon             , only : noveg
  implicit none
  class(betr_simulation_alm_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: betr_bounds
  type(column_type)     , intent(in)    :: col ! column type
  type(patch_type)      , intent(in) :: pft
  integer, intent(in) :: num_soilc
  integer, intent(in) :: filter_soilc(:)
  type(PlantMicKinetics_type), intent(in) :: PlantMicKinetics_vars

  integer :: j, fc, c, p, pi, pp, c_l

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
    minsurf_p_compet_vr_col => PlantMicKinetics_vars%minsurf_p_compet_vr_col &
  )
  c_l = 1
  do j =1, betr_bounds%ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      pp = 0
      do pi = 1, betr_maxpatch_pft
        if (pi <= col%npfts(c)) then
          p = col%pfti(c) + pi - 1
          if (pft%active(p) .and. (pft%itype(p) .ne. noveg)) then
            pp = pp + 1
            this%betr(c)%plantNutkinetics%plant_nh4_vmax_vr_patch(pp,j) = plant_nh4_vmax_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_no3_vmax_vr_patch(pp,j) = plant_no3_vmax_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_p_vmax_vr_patch(pp,j) = plant_p_vmax_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_nh4_km_vr_patch(pp,j) = plant_nh4_km_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_no3_km_vr_patch(pp,j) = plant_no3_km_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_p_km_vr_patch(pp,j) = plant_p_km_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_eff_ncompet_b_vr_patch(pp,j)=plant_eff_ncompet_b_vr_patch(p,j)
            this%betr(c)%plantNutkinetics%plant_eff_pcompet_b_vr_patch(pp,j)=plant_eff_pcompet_b_vr_patch(p,j)
          endif
        endif
      enddo
      this%betr(c)%nactpft = pp
      this%betr(c)%plantNutkinetics%minsurf_p_compet_vr_col(c_l,j) = minsurf_p_compet_vr_col(c,j)
      this%betr(c)%plantNutkinetics%minsurf_nh4_compet_vr_col(c_l,j) = minsurf_nh4_compet_vr_col(c,j)
    enddo
  enddo

  !the following parameters are specific to ECACNP, and I assume they are
  !grid specific as they currently used in alm-cnp.
  if(trim(reaction_method)=='eca_cnp')then
    do j =1, betr_bounds%ubj
      do fc = 1, num_soilc
        c = filter_soilc(fc)
        this%betr(c)%plantNutkinetics%km_minsurf_p_vr_col(c_l,j) = PlantMicKinetics_vars%km_minsurf_p_vr_col(c,j)
        this%betr(c)%plantNutkinetics%km_decomp_nh4_vr_col(c_l,j) = PlantMicKinetics_vars%km_decomp_nh4_vr_col(c,j)
        this%betr(c)%plantNutkinetics%km_decomp_no3_vr_col(c_l,j) = PlantMicKinetics_vars%km_decomp_no3_vr_col(c,j)
        this%betr(c)%plantNutkinetics%km_decomp_p_vr_col(c_l,j)=PlantMicKinetics_vars%km_decomp_p_vr_col(c,j)
        this%betr(c)%plantNutkinetics%km_nit_nh4_vr_col(c_l,j)=PlantMicKinetics_vars%km_nit_nh4_vr_col(c,j)
        this%betr(c)%plantNutkinetics%km_den_no3_vr_col(c_l,j)=PlantMicKinetics_vars%km_den_no3_vr_col(c,j)


        !effective p competing decomposers
        this%betr(c)%plantNutkinetics%decomp_eff_ncompet_b_vr_col(c_l,j) = PlantMicKinetics_vars%decomp_eff_ncompet_b_vr_col(c,j)
        this%betr(c)%plantNutkinetics%decomp_eff_pcompet_b_vr_col(c_l,j) = PlantMicKinetics_vars%decomp_eff_pcompet_b_vr_col(c,j)
        this%betr(c)%plantNutkinetics%den_eff_ncompet_b_vr_col(c_l,j) = PlantMicKinetics_vars%den_eff_ncompet_b_vr_col(c,j)
        this%betr(c)%plantNutkinetics%nit_eff_ncompet_b_vr_col(c_l,j) = PlantMicKinetics_vars%nit_eff_ncompet_b_vr_col(c,j)

      enddo
    enddo
  endif
  end associate
  end subroutine set_transient_kinetics_par
end module BeTRSimulationALM
