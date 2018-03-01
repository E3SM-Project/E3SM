module DIOCBGCReactionsType

#include "bshr_assert.h"
  !
  ! !DESCRIPTION:
  ! This is an example on how to use polymorphism to create your own bgc modules that will be run with BeTR
  !
  ! HISTORY:
  ! Created by Jinyun Tang, Oct 2nd, 2014
  ! !USES:
  use bshr_log_mod             , only : errMsg => shr_log_errMsg
  use bshr_kind_mod            , only : r8 => shr_kind_r8
  use BeTR_decompMod           , only : bounds_type  => betr_bounds_type
  use BGCReactionsMod          , only : bgc_reaction_type
  use tracer_varcon            , only : bndcond_as_conc, bndcond_as_flux
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use BetrStatusType           , only : betr_status_type
  implicit none

  private

  character(len=*), parameter :: mod_filename = &
       __FILE__

  public :: bgc_reaction_dioc_run_type

  type, extends(bgc_reaction_type) :: &
       bgc_reaction_dioc_run_type
  private
  contains
     procedure :: Init_betrbgc                          ! initialize betr bgc
     procedure :: set_boundary_conditions               ! set top/bottom boundary conditions for various tracers
     procedure :: calc_bgc_reaction                     ! doing bgc calculation
     procedure :: init_boundary_condition_type          ! initialize type of top boundary conditions
     procedure :: do_tracer_equilibration               ! do equilibrium tracer chemistry
     procedure :: InitCold                              ! do cold initialization
     procedure :: retrieve_biogeoflux           !
     procedure :: set_kinetics_par
     procedure, private :: readParams                   ! read in parameters
     procedure :: retrieve_lnd2atm
     procedure :: retrieve_biostates
     procedure :: debug_info
  end type bgc_reaction_dioc_run_type

  interface bgc_reaction_dioc_run_type
     module procedure constructor
  end interface bgc_reaction_dioc_run_type

contains
  !-------------------------------------------------------------------------------
  type(bgc_reaction_dioc_run_type) function constructor()
  !
  ! !DESCRIPTION:
  ! create an object of type bgc_reaction_dioc_run_type.
  ! Right now it is purposely empty
    type(bgc_reaction_dioc_run_type), allocatable :: bgc
    allocate(bgc)
    constructor = bgc
  end function constructor
  !----------------------------------------------------------------------
  subroutine set_kinetics_par(this, lbj, ubj, nactpft, plantNutkinetics)
  use PlantNutKineticsMod, only : PlantNutKinetics_type

  ! !ARGUMENTS:
  class(bgc_reaction_dioc_run_type)         , intent(inout)    :: this                       !
  class(PlantNutKinetics_type), intent(in) :: plantNutkinetics
  integer, intent(in) :: lbj, ubj
  integer, intent(in) :: nactpft

  end subroutine set_kinetics_par
  !-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
    !
    ! !DESCRIPTION:
    ! initialize boundary condition types
    !
    ! !USES:
    use BeTRTracerType        , only : betrtracer_type
    use TracerBoundaryCondType, only : tracerboundarycond_type
    use tracer_varcon         , only : bndcond_as_conc, bndcond_as_flux
    use BeTRTracerType        , only : betrtracer_type

    ! !ARGUMENTS:
    class(bgc_reaction_dioc_run_type), intent(inout) :: this
    type(BeTRtracer_type),             intent(in) :: betrtracer_vars
    type(bounds_type),                 intent(in) :: bounds
    type(tracerboundarycond_type),     intent(in) :: tracerboundarycond_vars

    ! !LOCAL VARIABLES:

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning) continue
    if (bounds%begc > 0) continue
    if (len(betrtracer_vars%betr_simname) > 0) continue
    tracerboundarycond_vars%topbc_type(:) = bndcond_as_conc

    !when bottom BC is not given, it is specified as constant flux
    ! FIXME(bja, 201604) Don't we need a bottom BC?
    !X!tracerboundarycond_vars%botbc_type(:) = bndcond_as_flux

  end subroutine init_boundary_condition_type

  !-------------------------------------------------------------------------------
  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars, namelist_buffer, bstatus)
    !
    ! DESCRIPTION:
    ! initialize the betrbgc
    !
    ! !USES:
    use BeTRTracerType , only : betrtracer_type
    use MathfuncMod    , only : addone
    use BetrStatusType , only : betr_status_type
    use gbetrType      , only : gbetr_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_dioc_run_type), intent(inout)    :: this
    type(bounds_type)                , intent(in)    :: bounds
    integer                          , intent(in)    :: lbj, ubj
    type(BeTRtracer_type )           , intent(inout) :: betrtracer_vars
    character(len=*)                 , intent(in)    :: namelist_buffer
    type(betr_status_type)           , intent(out)   :: bstatus
    character(len=*), parameter                      :: subname ='Init_betrbgc'

    integer :: itemp_gwm
    integer :: itemp_g
    integer :: itemp_s
    integer :: itemp_gwm_grp
    integer :: dum, itemp
    integer :: itemp_grp, itemp_v, itemp_vgrp, itemp_trc

    call bstatus%reset()
    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)           continue
    if (bounds%begc > 0)                       continue
    if (ubj > lbj)                             continue
    if (len(betrtracer_vars%betr_simname) > 0) continue

    itemp_gwm     = 0;
    itemp_g       = 0 ;
    itemp_s       = 0;
    itemp_gwm_grp = 0

    !volatile tracers
    itemp = 0; itemp_trc=0
    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_n2, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_n2,&
      trc_grp_end=betrtracer_vars%id_trc_end_n2, is_trc_gw=.true., is_trc_volatile = .true.)

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_o2, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_o2, &
      trc_grp_end=betrtracer_vars%id_trc_end_o2, &
      is_trc_gw=.true., is_trc_volatile = .true.)

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_ar,&
      trc_grp_beg= betrtracer_vars%id_trc_beg_ar, &
      trc_grp_end= betrtracer_vars%id_trc_end_ar, &
      is_trc_gw=.true., is_trc_volatile = .true.)

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_co2x, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_co2x, &
      trc_grp_end=betrtracer_vars%id_trc_end_co2x, &
      is_trc_gw=.true., is_trc_volatile = .true.)

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp= betrtracer_vars%id_trc_ch4, &
      trc_grp_beg= betrtracer_vars%id_trc_beg_ch4, &
      trc_grp_end= betrtracer_vars%id_trc_end_ch4, &
      is_trc_gw=.true., is_trc_volatile = .true.)

    call betrtracer_vars%add_tracer_group(trc_grp_cnt=addone(itemp), mem = 1, &
      trc_cnt=itemp_trc, trc_grp=betrtracer_vars%id_trc_doc, &
      trc_grp_beg=betrtracer_vars%id_trc_beg_doc, &
      trc_grp_end=betrtracer_vars%id_trc_end_doc, &
      is_trc_gw=.true., is_trc_volatile = .false.)


    betrtracer_vars%nmem_max               = 1

    call betrtracer_vars%Init()

    itemp_grp = 0    !group id
    itemp_v = 0      !volatile id
    itemp_vgrp = 0   !volatile group

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_n2, trc_name='N2'  ,      &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem= 1,  is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_o2, trc_name='O2'  ,      &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ar, trc_name='AR'  ,      &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_co2x, trc_name='CO2x',    &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp)  , &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)       , &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_ch4, trc_name='CH4',      &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)     ,   &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(bstatus=bstatus,trc_id = betrtracer_vars%id_trc_doc, trc_name='DOC',      &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem = 1)

  end subroutine Init_betrbgc

  !-------------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, dz_top, betrtracer_vars, &
       biophysforc, biogeo_flux, tracerboundarycond_vars, betr_status)
    !
    ! !DESCRIPTION:
    ! set up boundary conditions for tracer movement
    !
    ! !USES:
    use betr_ctrl              , only : iulog  => biulog
    use TracerBoundaryCondType , only : tracerboundarycond_type
    use bshr_log_mod           , only : errMsg => shr_log_errMsg
    use BeTRTracerType         , only : betrtracer_type
    use betr_varcon            , only : rgas => brgas
    use BeTR_biogeoFluxType    , only : betr_biogeo_flux_type
    use BetrStatusType         , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_dioc_run_type) , intent(inout)    :: this                       !
    type(bounds_type)                 , intent(in)    :: bounds                     !
    integer                           , intent(in)    :: num_soilc                  ! number of columns in column filter_soilc
    integer                           , intent(in)    :: filter_soilc(:)            ! column filter_soilc
    type(betrtracer_type)             , intent(in)    :: betrtracer_vars            !
    real(r8)                          , intent(in)    :: dz_top(bounds%begc: )      !
    type(betr_biogeophys_input_type)  , intent(in)    :: biophysforc
    type(betr_biogeo_flux_type)       , intent(in)    :: biogeo_flux
    type(tracerboundarycond_type)     , intent(inout) :: tracerboundarycond_vars !
    type(betr_status_type)            , intent(out)   :: betr_status

    ! !LOCAL VARIABLES:
    integer            :: fc, c
    character(len=255) :: subname = 'set_boundary_conditions'
    real(r8) :: irt   !the inverse of R*T

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(dz_top)  == (/bounds%endc/)),   errMsg(mod_filename,__LINE__), betr_status)
    if(betr_status%check_status())return
    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)      continue
    if (size(biogeo_flux%qflx_adv_col)>0) continue
    associate(                                                             &
         forc_pbot            => biophysforc%forc_pbot_downscaled_col    , &
         forc_tbot            => biophysforc%forc_t_downscaled_col       , &
         groupid              => betrtracer_vars%groupid                   &
         )

      do fc = 1, num_soilc
         c = filter_soilc(fc)
         irt = 1.e3_r8/(forc_tbot(c)*rgas)
         !eventually, the following code will be implemented using polymorphism
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_n2)   = forc_pbot(c)*0.78084_r8*irt  !mol m-3, contant boundary condition, as concentration
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_o2)   = forc_pbot(c)*0.20946_r8*irt  !mol m-3, contant boundary condition, as concentration
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ar)   = forc_pbot(c)*0.009340_r8*irt !mol m-3, contant boundary condition, as concentration
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_co2x) = forc_pbot(c)*367e-6_r8*irt   !mol m-3, contant boundary condition, as concentration
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ch4)  = forc_pbot(c)*1.79e-6_r8*irt  !mol m-3, contant boundary condition, as concentration
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_doc)  = 0._r8                        !mol m-3, contant boundary condition, as concentration

         tracerboundarycond_vars%bot_concflux_col(c,1,:)                                          = 0._r8                        !zero flux boundary condition for diffusion
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_n2))           = 2._r8*1.837e-5_r8/dz_top(c)  !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_o2))           = 2._r8*1.713e-5_r8/dz_top(c)  !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_ar))           = 2._r8*1.532e-5_r8/dz_top(c)  !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_co2x))         = 2._r8*1.399e-5_r8/dz_top(c)  !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_ch4))          = 2._r8*1.808e-5_r8/dz_top(c)  !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_doc))          = 0._r8                        !m/s surface conductance
      enddo

    end associate
  end subroutine set_boundary_conditions

  !-------------------------------------------------------------------------------
  subroutine calc_bgc_reaction(this, bounds, col, lbj, ubj, num_soilc, filter_soilc,               &
       num_soilp,filter_soilp, jtops, dtime, betrtracer_vars, tracercoeff_vars,  biophysforc, &
       tracerstate_vars, tracerflux_vars, tracerboundarycond_vars, plant_soilbgc, &
       biogeo_flux,  betr_status)
    !
    ! !DESCRIPTION:
    ! do bgc reaction
    !
    ! !USES:
    use TracerBoundaryCondType , only : tracerboundarycond_type
    use tracerfluxType         , only : tracerflux_type
    use tracerstatetype        , only : tracerstate_type
    use tracercoeffType        , only : tracercoeff_type
    use BetrTracerType         , only : betrtracer_type
    use PlantSoilBGCMod        , only : plant_soilbgc_type
    use BetrStatusType         , only : betr_status_type
    use betr_columnType        , only : betr_column_type
    use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
    use BeTR_biogeoStateType     , only : betr_biogeo_state_type
    implicit none
    !ARGUMENTS
    class(bgc_reaction_dioc_run_type) , intent(inout) :: this                       !
    type(bounds_type)                 , intent(in)    :: bounds                      ! bounds
    type(betr_column_type)            , intent(in)    :: col
    integer                           , intent(in)    :: num_soilc                   ! number of columns in column filter
    integer                           , intent(in)    :: filter_soilc(:)             ! column filter
    integer                           , intent(in)    :: num_soilp
    integer                           , intent(in)    :: filter_soilp(:)
    integer                           , intent(in)    :: jtops( : )                  ! top index of each column
    integer                           , intent(in)    :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
    real(r8)                          , intent(in)    :: dtime                       ! model time step
    type(betrtracer_type)             , intent(in)    :: betrtracer_vars             ! betr configuration information
    type(betr_biogeophys_input_type)  , intent(in)    :: biophysforc
    type(tracercoeff_type)            , intent(in)    :: tracercoeff_vars
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars
    type(tracerflux_type)             , intent(inout) :: tracerflux_vars
    type(tracerboundarycond_type)     , intent(inout) :: tracerboundarycond_vars !
    class(plant_soilbgc_type)         , intent(inout) ::  plant_soilbgc
    type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux
    type(betr_status_type)            , intent(out)   :: betr_status
    character(len=*)                 , parameter     :: subname ='calc_bgc_reaction'

    integer :: c, fc, ll

    call betr_status%reset()
    associate(                                                                    &
    tracer_mobile_phase            => tracerstate_vars%tracer_conc_mobile_col  ,  &
    tracer_flx_netpro_vr           => tracerflux_vars%tracer_flx_netpro_vr_col ,  &
    dic_prod_vr                    => biophysforc%dic_prod_vr_col              ,  &
    doc_prod_vr                    => biophysforc%doc_prod_vr_col              ,  &
    id_trc_doc                     => betrtracer_vars%id_trc_doc               ,  &
    id_trc_co2x                    => betrtracer_vars%id_trc_co2x                 &
    )
    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)                          continue
    if (bounds%begc > 0)                                      continue
    if (ubj > lbj)                                            continue
    if (size(jtops) > 0)                                      continue
    if (num_soilp > 0)                                        continue
    if (size(filter_soilp) > 0)                               continue
    if (dtime > 0.0)                                          continue
    if (len(betrtracer_vars%betr_simname) > 0)                continue
    if (size(biophysforc%isoilorder) > 0)                     continue
    if (size(tracercoeff_vars%annsum_counter_col) > 0)        continue
    if (size(tracerstate_vars%tracer_conc_surfwater_col) > 0) continue
    if (size(tracerflux_vars%tracer_flx_top_soil_col) > 0)    continue
    if (size(tracerboundarycond_vars%jtops_col) > 0)          continue
    if (plant_soilbgc%dummy_compiler_warning)                 continue

    !now assume doc decays with a turnover rate 1.e-6_r8
    do ll = 1, ubj
      do fc = 1, num_soilc
        c = filter_soilc(fc)
        tracer_flx_netpro_vr(c,ll,id_trc_doc)= doc_prod_vr(c,ll)*dtime
        tracer_mobile_phase(c,ll,id_trc_doc) = tracer_mobile_phase(c,ll,id_trc_doc)+tracer_flx_netpro_vr(c,ll,id_trc_doc)
        if(tracer_mobile_phase(c,ll,id_trc_doc)<0._r8)then
           tracer_flx_netpro_vr(c,ll,id_trc_doc)= tracer_flx_netpro_vr(c,ll,id_trc_doc)-tracer_mobile_phase(c,ll,id_trc_doc)
           tracer_mobile_phase(c,ll,id_trc_doc)=0._r8
        endif
        tracer_flx_netpro_vr(c,ll,id_trc_co2x)= dic_prod_vr(c,ll)*dtime
        tracer_mobile_phase(c,ll,id_trc_co2x) = tracer_mobile_phase(c,ll,id_trc_co2x)+tracer_flx_netpro_vr(c,ll,id_trc_co2x)
      enddo
    enddo
   end associate
  end subroutine calc_bgc_reaction

  !-------------------------------------------------------------------------------
  subroutine do_tracer_equilibration(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       betrtracer_vars, tracercoeff_vars, tracerstate_vars, betr_status)
    !
    ! DESCRIPTION:
    ! requilibrate tracers that has solid and mobile phases
    ! using the theory of mass action.
    !
    ! !USES:
    !
    use tracerstatetype , only : tracerstate_type
    use tracercoeffType , only : tracercoeff_type
    use BeTRTracerType  , only : betrtracer_type
    use BetrStatusType  , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_dioc_run_type) , intent(inout)    :: this
    type(bounds_type)                 , intent(in)    :: bounds
    integer                           , intent(in)    :: lbj, ubj
    integer                           , intent(in)    :: jtops(bounds%begc: )        ! top label of each column
    integer                           , intent(in)    :: num_soilc
    integer                           , intent(in)    :: filter_soilc(:)
    type(betrtracer_type)             , intent(in)    :: betrtracer_vars
    type(tracercoeff_type)            , intent(in)    :: tracercoeff_vars
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars
    type(betr_status_type)            , intent(out)   :: betr_status
    !local variables
    character(len=255) :: subname = 'do_tracer_equilibration'

    call betr_status%reset()
    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(mod_filename,__LINE__), betr_status)
    if(betr_status%check_status())return

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)                          continue
    if (bounds%begc > 0)                                      continue
    if (ubj > lbj)                                            continue
    if (size(jtops) > 0)                                      continue
    if (num_soilc > 0)                                        continue
    if (size(filter_soilc) > 0)                               continue
    if (len(betrtracer_vars%betr_simname) > 0)                continue
    if (size(tracerstate_vars%tracer_conc_surfwater_col) > 0) continue
    if (size(tracercoeff_vars%annsum_counter_col) > 0)        continue

    !continue on the simulation type, an implementation of aqueous chemistry will be
    !employed to separate out the adsorbed phase
    !It should be noted that this formulation excludes the use of linear isotherm, which
    !can be integrated through the retardation factor

  end subroutine do_tracer_equilibration

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, col, betrtracer_vars, biophysforc, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! do cold initialization
    !
    ! !USES:
    use BeTRTracerType      , only : BeTRTracer_Type
    use tracerstatetype     , only : tracerstate_type
    use betr_varcon         , only : spval => bspval, ispval => bispval
    use BeTR_landvarconType , only : landvarcon  => betr_landvarcon
    use betr_columnType     , only : betr_column_type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_dioc_run_type) , intent(inout)    :: this
    type(bounds_type)                 , intent(in)    :: bounds
    type(betr_column_type)            , intent(in)    :: col
    type(BeTRTracer_Type)             , intent(in)    :: betrtracer_vars
    type(betr_biogeophys_input_type)  , intent(in)    :: biophysforc
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars

    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l, k, j
    integer :: fc                                        ! filter_soilc index
    integer               :: begc, endc
    integer               :: begg, endg
    !-----------------------------------------------------------------------

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)          continue
    if (size(biophysforc%h2osoi_liq_col) > 0) continue

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg
    !-----------------------------------------------------------------------


    do c = bounds%begc, bounds%endc

      !dual phase tracers

      tracerstate_vars%tracer_conc_mobile_col(c,:, :)          = 0._r8
      tracerstate_vars%tracer_conc_surfwater_col(c,:)          = 0._r8
      tracerstate_vars%tracer_conc_aquifer_col(c,:)            = 0._r8
      tracerstate_vars%tracer_conc_grndwater_col(c,:)          = 0._r8
      tracerstate_vars%tracer_conc_mobile_col(c,7, betrtracer_vars%id_trc_doc) = 0._r8  !point source

      if(betrtracer_vars%nsolid_equil_tracers>0)then
        tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = 0._r8
      endif
      tracerstate_vars%tracer_soi_molarmass_col(c,:)          = 0._r8
    enddo

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine readParams(this, name_list_buffer, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! read in module specific parameters
    !
    ! !USES:
    use BeTRTracerType , only : BeTRTracer_Type
    implicit none
    ! !ARGUMENTS:
    class(bgc_reaction_dioc_run_type) , intent(inout)    :: this
    type(BeTRTracer_Type)             , intent(inout) :: betrtracer_vars
    character(len=*)                  , intent(in)  :: name_list_buffer

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)           continue
    if (len(betrtracer_vars%betr_simname) > 0) continue

    !do nothing here
  end subroutine readParams

  !-------------------------------------------------------------------------------
  subroutine retrieve_biogeoflux(this, num_soilc, filter_soilc, tracerflux_vars, &
  betrtracer_vars, biogeo_flux)

  use tracerfluxType           , only : tracerflux_type
  use BeTR_decompMod           , only : betr_bounds_type
  use BeTRTracerType           , only : BeTRTracer_Type
  use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
  implicit none
   class(bgc_reaction_dioc_run_type) , intent(inout) :: this
  integer                          , intent(in)    :: num_soilc                   ! number of columns in column filter
  integer                          , intent(in)    :: filter_soilc(:)             ! column filter
  type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
  type(tracerflux_type)            , intent(in)    :: tracerflux_vars
  type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux

    ! remove compiler warnings for unused dummy args
    if (this%dummy_compiler_warning)                          continue
    if (num_soilc > 0)                                        continue
    if (size(filter_soilc) > 0)                               continue
    if (len(betrtracer_vars%betr_simname) > 0)                continue
    if (size(tracerflux_vars%tracer_flx_top_soil_col) > 0)    continue

  end subroutine retrieve_biogeoflux


   !----------------------------------------------------------------------
   subroutine retrieve_lnd2atm(this, bounds, num_soilc, filter_soilc, tracerflux_vars, &
   betrtracer_vars, biogeo_flux)

   use tracerfluxType           , only : tracerflux_type
   use BeTR_decompMod           , only : betr_bounds_type
   use BeTRTracerType           , only : BeTRTracer_Type
   use BeTR_biogeoFluxType      , only : betr_biogeo_flux_type
   implicit none
   class(bgc_reaction_dioc_run_type) , intent(inout) :: this               !
   type(betr_bounds_type)           , intent(in)    :: bounds                      ! bounds
   integer                          , intent(in)    :: num_soilc                   ! number of columns in column filter
   integer                          , intent(in)    :: filter_soilc(:)             ! column filter
   type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
   type(tracerflux_type)            , intent(in)    :: tracerflux_vars
   type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux

   integer :: c, fc

   associate( &
    id_trc_doc                     => betrtracer_vars%id_trc_doc               ,  &
    id_trc_co2x                    => betrtracer_vars%id_trc_co2x                 &
   )
   do fc = 1, num_soilc
     c = filter_soilc(fc)
     biogeo_flux%qflx_rofliq_qsur_doc_col(c) = max(tracerflux_vars%tracer_flx_surfrun_col(c,id_trc_doc),0._r8)*12._r8
     biogeo_flux%qflx_rofliq_qsur_dic_col(c) = max(tracerflux_vars%tracer_flx_surfrun_col(c,id_trc_co2x),0._r8)*12._r8
     biogeo_flux%qflx_rofliq_qsub_doc_col(c) = max(tracerflux_vars%tracer_flx_leaching_col(c,id_trc_doc) + &
        tracerflux_vars%tracer_flx_drain_col(c,id_trc_doc),0._r8)*12._r8
     biogeo_flux%qflx_rofliq_qsub_dic_col(c) = max(tracerflux_vars%tracer_flx_leaching_col(c,id_trc_co2x) + &
        tracerflux_vars%tracer_flx_drain_col(c,id_trc_co2x),0._r8)*12._r8
   enddo
   end associate
   end subroutine retrieve_lnd2atm

   !-------------------------------------------------------------------------------
   subroutine debug_info(this, bounds, num_soilc, filter_soilc, dzsoi, betrtracer_vars, tracerstate_vars,  header, betr_status)

   use BeTRTracerType           , only : BeTRTracer_Type
   use tracerstatetype          , only : tracerstate_type
   use BeTR_decompMod           , only : betr_bounds_type
     ! !ARGUMENTS:
   implicit none
   class(bgc_reaction_dioc_run_type) , intent(inout) :: this     !
   type(betr_bounds_type)               , intent(in) :: bounds                      ! bounds
   integer                              , intent(in) :: num_soilc                   ! number of columns in column filter
   integer                              , intent(in) :: filter_soilc(:)             ! column filter
   real(r8)                             , intent(in) :: dzsoi(bounds%begc: ,bounds%lbj: )
   type(betrtracer_type)                , intent(in) :: betrtracer_vars             ! betr configuration information
   type(tracerstate_type)               , intent(in) :: tracerstate_vars
   character(len=*)                     , intent(in) :: header
   type(betr_status_type)               , intent(out):: betr_status
   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(dzsoi)==(/bounds%endc,bounds%ubj/)), errMsg(mod_filename,__LINE__), betr_status)
   if(betr_status%check_status())return

   if (this%dummy_compiler_warning) continue
     end subroutine debug_info
   !----------------------------------------------------------------------
   subroutine retrieve_biostates(this, bounds, lbj, ubj, jtops,num_soilc, filter_soilc, &
      betrtracer_vars, tracerstate_vars, biogeo_state, betr_status)
   !
   !retrieve state variables for lsm mass balance check
   use tracer_varcon, only : catomw, natomw, patomw, c13atomw, c14atomw
   use BeTRTracerType           , only : BeTRTracer_Type
   use tracerstatetype          , only : tracerstate_type
   use BeTR_biogeoStateType     , only : betr_biogeo_state_type
   use BeTR_decompMod           , only : betr_bounds_type
   implicit none
   class(bgc_reaction_dioc_run_type) , intent(inout) :: this               !
   type(betr_bounds_type)               , intent(in)  :: bounds                      ! bounds
   integer                              , intent(in) :: lbj, ubj
   integer                              , intent(in) :: jtops(bounds%begc: )
   integer                              , intent(in)    :: num_soilc                   ! number of columns in column filter
   integer                              , intent(in)    :: filter_soilc(:)             ! column filter
   type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
   type(tracerstate_type)               , intent(inout) :: tracerstate_vars
   type(betr_biogeo_state_type)         , intent(inout) :: biogeo_state
   type(betr_status_type)               , intent(out):: betr_status

   call betr_status%reset()
   SHR_ASSERT_ALL((ubound(jtops)  == (/bounds%endc/)),   errMsg(mod_filename,__LINE__), betr_status)
   if(betr_status%check_status())return
   if (this%dummy_compiler_warning) continue
   if (bounds%begc > 0)             continue

   end subroutine retrieve_biostates
end module DIOCBGCReactionsType
