module BGCReactionsCenturyType

#include "shr_assert.h"

  !
  ! !DESCRIPTION:
  ! this code uses the operator automated down-regulation scheme described in Tang and Riley, 2015, BG
  ! HISTORY:
  ! Created by Jinyun Tang, Oct 2nd, 2014

  ! !USES
  !
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use clm_varctl            , only : iulog
  use abortutils            , only : endrun
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use clm_time_manager      , only : get_nstep
  use shr_infnan_mod        , only : nan => shr_infnan_nan, assignment(=)
  use decompMod             , only : bounds_type
  use BGCReactionsMod       , only : bgc_reaction_type
  use clm_varcon            , only : spval
  use clm_varctl            , only : spinup_state
  use tracer_varcon         , only : bndcond_as_conc, bndcond_as_flux
  use BGCCenturySubMod
  use BGCCenturySubCoreMod
  use LandunitType          , only : lun
  use ColumnType            , only : col
  use GridcellType          , only : grc
  use landunit_varcon       , only : istsoil, istcrop
  implicit none

  save
  private
  !
  ! !PUBLIC TYPES:
  public :: bgc_reaction_CENTURY_type
  type(centurybgc_type), private :: centurybgc_vars

  !integer, private :: lpr
  type, extends(bgc_reaction_type) :: &
       bgc_reaction_CENTURY_type
  private
contains
  procedure :: Init_betrbgc                    ! initialize betr bgc
  procedure :: set_boundary_conditions         ! set top/bottom boundary conditions for various tracers
  procedure :: calc_bgc_reaction               ! doing bgc calculation
  procedure :: init_boundary_condition_type    ! initialize type of top boundary conditions
  procedure :: do_tracer_equilibration         ! do equilibrium tracer chemistry
  procedure :: initCold
  procedure :: readParams
  procedure :: init_betr_alm_bgc_coupler       ! update state vars using other bgc parts in alm
  procedure :: betr_alm_flux_statevar_feedback !
end type bgc_reaction_CENTURY_type


type, private :: Extra_type
  real(r8), pointer :: cn_ratios(:)           !cn ratio of om pool
  real(r8), pointer :: cp_ratios(:)           !cp ratio of om pool
  real(r8), pointer :: k_decay(:)             !decay parameter for all reactions
  real(r8), pointer :: scal_f(:)              !scaling factor for first order sink
  real(r8), pointer :: conv_f(:)              !converting factor for first order sink
  real(r8), pointer :: conc_f(:)              !external forcing strength
  real(r8)          :: n2_n2o_ratio_denit     !ratio of n2 to n2o during denitrification
  real(r8)          :: cellsand               !sand content
  logical,  pointer :: is_zero_order(:)
  integer           :: nr                     !number of reactions involved
contains
  procedure, public :: Init_Allocate
  procedure, public :: DDeallocate
  procedure, public :: AAssign
end type Extra_type
type(Extra_type), private :: Extra_inst


  interface bgc_reaction_CENTURY_type
     module procedure constructor
  end interface bgc_reaction_CENTURY_type


contains

  subroutine Init_Allocate(this, nompools, nreacts, nprimstvars)

    !
    ! !DESCRIPTION:
    ! do memory allocation for the data type specified by this
    !
    ! !ARGUMENTS:
    class(Extra_type) :: this

    integer, intent(in) :: nompools
    integer, intent(in) :: nreacts
    integer, intent(in) :: nprimstvars     !number of primary state variables

    allocate(this%cn_ratios(nompools))
    allocate(this%cp_ratios(nompools))
    allocate(this%k_decay(nreacts))
    allocate(this%scal_f(nprimstvars    ));    this%scal_f(:) = 0._r8
    allocate(this%conv_f(nprimstvars    ));    this%conv_f(:) = 0._r8
    allocate(this%conc_f(nprimstvars    ));    this%conc_f(:) = 0._r8
    allocate(this%is_zero_order(nreacts )); this%is_zero_order(:) = .false.
    this%nr = nreacts

  end subroutine Init_Allocate

  !-------------------------------------------------------------------------------

  subroutine DDeallocate(this)
    !
    ! !DESCRIPTION:
    ! deallocate memory for the data type specified by this
    ! !ARGUMENTS:
    class(Extra_type) :: this


    deallocate(this%cn_ratios)
    deallocate(this%cp_ratios)
    deallocate(this%k_decay)
    deallocate(this%scal_f)
    deallocate(this%conv_f)
    deallocate(this%conc_f)

  end subroutine DDeallocate
  !-------------------------------------------------------------------------------

  subroutine AAssign(this, cn_r,cp_r, k_d,  n2_n2o_r_denit, cell_sand, &
       betrtracer_vars, gas2bulkcef, aere_cond, tracer_conc_atm)
    !
    ! !DESCRIPTION:
    ! assign memmber values for the data type specified by this
    ! !USES:
    use BeTRTracerType              , only : betrtracer_type
    ! !ARGUMENTS:
    class(Extra_type) :: this
    real(r8), dimension(:), intent(in) :: cn_r
    real(r8), dimension(:), intent(in) :: cp_r
    real(r8), dimension(:), intent(in) :: k_d
    real(r8)              , intent(in) :: n2_n2o_r_denit
    real(r8)              , intent(in) :: cell_sand
    type(BeTRtracer_type ), intent(in) :: betrtracer_vars
    real(r8)              , intent(in) :: gas2bulkcef(1:betrtracer_vars%nvolatile_tracers)
    real(r8)              , intent(in) :: aere_cond(1:betrtracer_vars%nvolatile_tracers)
    real(r8)              , intent(in) :: tracer_conc_atm(1:betrtracer_vars%nvolatile_tracers)

    ! !LOCAL VARIABLES:
    integer :: n1, n2, n3, j


    n1 = size(cn_r)
    n2 = size(cp_r)
    n3 = size(k_d)

    SHR_ASSERT_ALL((n1 == n2),        errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((n3 == this%nr),   errMsg(__FILE__,__LINE__))

    this%cn_ratios(1:n1)    = cn_r
    this%cp_ratios(1:n2)    = cp_r
    this%n2_n2o_ratio_denit = n2_n2o_r_denit
    this%cellsand           = cell_sand
    this%k_decay            = k_d


    do j = 1, betrtracer_vars%ngwmobile_tracers
       if(j == betrtracer_vars%id_trc_o2)then
          this%scal_f(centurybgc_vars%lid_o2)   = aere_cond(betrtracer_vars%volatileid(j))
          this%conc_f(centurybgc_vars%lid_o2)   = tracer_conc_atm(betrtracer_vars%volatileid(j))
          this%conv_f(centurybgc_vars%lid_o2)   = 1._r8/gas2bulkcef(betrtracer_vars%volatileid(j))

       elseif(j == betrtracer_vars%id_trc_n2)then
          this%scal_f(centurybgc_vars%lid_n2) = aere_cond(betrtracer_vars%volatileid(j))
          this%conc_f(centurybgc_vars%lid_n2) = tracer_conc_atm(betrtracer_vars%volatileid(j))
          this%conv_f(centurybgc_vars%lid_n2)   = 1._r8/gas2bulkcef(betrtracer_vars%volatileid(j))

       elseif(j == betrtracer_vars%id_trc_ar)then
          this%scal_f(centurybgc_vars%lid_ar) = aere_cond(betrtracer_vars%volatileid(j))
          this%conc_f(centurybgc_vars%lid_ar) = tracer_conc_atm(betrtracer_vars%volatileid(j))
          this%conv_f(centurybgc_vars%lid_ar)   = 1._r8/gas2bulkcef(betrtracer_vars%volatileid(j))

       elseif(j==betrtracer_vars%id_trc_co2x)then
          this%scal_f(centurybgc_vars%lid_co2)  = aere_cond(betrtracer_vars%volatileid(j))
          this%conc_f(centurybgc_vars%lid_co2) = tracer_conc_atm(betrtracer_vars%volatileid(j))
          this%conv_f(centurybgc_vars%lid_co2)   = 1._r8/gas2bulkcef(betrtracer_vars%volatileid(j))

       elseif(j==betrtracer_vars%id_trc_ch4) then
          this%scal_f(centurybgc_vars%lid_ch4) = aere_cond(betrtracer_vars%volatileid(j))
          this%conc_f(centurybgc_vars%lid_ch4) = tracer_conc_atm(betrtracer_vars%volatileid(j))
          this%conv_f(centurybgc_vars%lid_ch4)   = 1._r8/gas2bulkcef(betrtracer_vars%volatileid(j))

       elseif(j==betrtracer_vars%id_trc_n2o) then
          this%scal_f(centurybgc_vars%lid_n2o) = aere_cond(betrtracer_vars%volatileid(j))
          this%conc_f(centurybgc_vars%lid_n2o) = tracer_conc_atm(betrtracer_vars%volatileid(j))
          this%conv_f(centurybgc_vars%lid_n2o) = 1._r8/gas2bulkcef(betrtracer_vars%volatileid(j))
       endif
    enddo
  end subroutine AAssign

  !-------------------------------------------------------------------------------
  type(bgc_reaction_CENTURY_type) function constructor()
    !
    ! !DESCRIPTION:
    !
    ! create an object of type bgc_reaction_CENTURY_type.
    ! Right now it is purposely empty

  end function constructor


  !-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
    !
    ! !DESCRIPTION:
    ! initialize boundary condition types
    ! !USES:
    use TracerBoundaryCondType      , only : tracerboundarycond_type
    use BeTRTracerType              , only : betrtracer_type

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_type), intent(in) :: this
    type(bounds_type)               , intent(in) :: bounds
    type(BeTRtracer_type )          , intent(in) :: betrtracer_vars
    type(tracerboundarycond_type)   , intent(in) :: tracerboundarycond_vars

    ! !LOCAL VARIABLES:
    integer :: c

    associate(                               &
         groupid  => betrtracer_vars%groupid &
         )


      tracerboundarycond_vars%topbc_type(1:betrtracer_vars%ngwmobile_tracer_groups                                )  = bndcond_as_conc
      tracerboundarycond_vars%topbc_type(groupid(betrtracer_vars%id_trc_no3x)                                     ) = bndcond_as_flux
      tracerboundarycond_vars%topbc_type(betrtracer_vars%ngwmobile_tracer_groups+1:betrtracer_vars%ntracer_groups )  = bndcond_as_flux

    end associate

  end subroutine init_boundary_condition_type

  !-------------------------------------------------------------------------------
  
  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! initialize the betrbgc
    ! !USES:
    use CNSharedParamsMod     , only : CNParamsReadShared
    use ncdio_pio             , only : file_desc_t
    use BeTRTracerType        , only : betrtracer_type
    use MathfuncMod           , only : addone
    use clm_varctl            , only : cnallocate_carbon_only_set
    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_type), intent(in)    :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0
    type(BeTRtracer_type )          , intent(inout) :: betrtracer_vars

    ! !LOCAL VARIABLES:
    character(len=32)               , parameter     :: subname ='Init_betrbgc'
    integer                                         :: jj
    integer                                         :: nelm, itemp_mem
    integer                                         :: itemp, itemp_vgrp, itemp_v, itemp_grp
    integer                                         :: c_loc, n_loc, trcid
    logical                                         :: carbon_only = .false.

    call cnallocate_carbon_only_set(carbon_only)
    call centurybgc_vars%Init(bounds, lbj, ubj)

    nelm =centurybgc_vars%nelms
    c_loc=centurybgc_vars%c_loc
    n_loc=centurybgc_vars%n_loc

    itemp = 0
    betrtracer_vars%id_trc_n2   = addone(itemp)
    betrtracer_vars%id_trc_o2   = addone(itemp)
    betrtracer_vars%id_trc_ar   = addone(itemp)
    betrtracer_vars%id_trc_co2x = addone(itemp)
    betrtracer_vars%id_trc_ch4  = addone(itemp)
    betrtracer_vars%id_trc_nh3x = addone(itemp)
    betrtracer_vars%id_trc_no3x = addone(itemp)
    betrtracer_vars%id_trc_n2o  = addone(itemp)

    betrtracer_vars%ngwmobile_tracer_groups      = itemp                          ! n2, o2, ar, co2, ch4, n2o, nh3x and no3x
    betrtracer_vars%ngwmobile_tracers            = itemp
    betrtracer_vars%nvolatile_tracers            = itemp-2                        ! n2, o2, ar, co2, ch4 and n2o
    betrtracer_vars%nvolatile_tracer_groups      = itemp-2                        !
    betrtracer_vars%nsolid_passive_tracer_groups = 4                              ! som1, som2, som3 and others (lit1, lit2, lit3, cwd)
    betrtracer_vars%nsolid_passive_tracers       = centurybgc_vars%nom_pools*nelm !

    betrtracer_vars%nmem_max                     = nelm*4                         ! total number of elemnts, and 4 sub members (lit1, lit2, lit3, cwd)

    call betrtracer_vars%Init()

    betrtracer_vars%is_mobile(:) = .true.

    jj         = itemp
    itemp_vgrp = 0  !counter for volatile groups
    itemp_v    = 0  !counter for volatile tracers
    itemp_grp  = 0  !counter for tracer groups

    trcid = betrtracer_vars%id_trc_n2

    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_n2, trc_name='N2'   ,     &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = addone(itemp_grp),  &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)      ,  &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_o2, trc_name='O2'   ,     &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = addone(itemp_grp),  &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)      ,  &
         trc_volatile_group_id = addone(itemp_vgrp))


    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_ar, trc_name='AR'    ,    &
         is_trc_mobile=.false., is_trc_advective = .false., trc_group_id = addone(itemp_grp), &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)       , &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_co2x, trc_name='CO2x',    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = addone(itemp_grp) , &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)       , &
         trc_volatile_group_id = addone(itemp_vgrp))

    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_ch4, trc_name='CH4'  ,    &
         is_trc_mobile=.false., is_trc_advective = .false., trc_group_id = addone(itemp_grp), &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)       , &
         trc_volatile_group_id = addone(itemp_vgrp))


    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_nh3x, trc_name='NH3x',    &
         is_trc_mobile=.false., is_trc_advective = .false., trc_group_id = addone(itemp_grp), &
         trc_group_mem = 1, is_trc_volatile=.false.)


    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_no3x, trc_name='NO3x',    &
         is_trc_mobile=.true., is_trc_advective = .true., trc_group_id = addone(itemp_grp),   &
         trc_group_mem = 1, is_trc_volatile=.false.,trc_vtrans_scal=1._r8)


    call betrtracer_vars%set_tracer(trc_id = betrtracer_vars%id_trc_n2o, trc_name='N2O' ,     &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = addone(itemp_grp),  &
         trc_group_mem = 1, is_trc_volatile=.true., trc_volatile_id = addone(itemp_v)      ,  &
         trc_volatile_group_id = addone(itemp_vgrp))



    !------------------------------------------------------------------------------------
    itemp_mem=0
    itemp_grp=addone(itemp_grp)          !only one group passive solid litter tracers
    trcid = jj+(centurybgc_vars%lit1-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT1C'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%lit1-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT1N'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))

    !------------------------------------------------------------------------------------

    trcid = jj+(centurybgc_vars%lit2-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT2C'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%lit2-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT2N'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))

    !------------------------------------------------------------------------------------

    trcid = jj+(centurybgc_vars%lit3-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT3C'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%lit3-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT3N'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))
    !------------------------------------------------------------------------------------

    trcid = jj+(centurybgc_vars%cwd-1 )*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='CWDC'              ,            &
         is_trc_mobile=.false., is_trc_advective = .false., trc_group_id = itemp_grp,         &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%cwd-1 )*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='CWDN'              ,            &
         is_trc_mobile=.false., is_trc_advective = .false., trc_group_id = itemp_grp,         &
         trc_group_mem= addone(itemp_mem))
    !==========================================================================================
    !new group
    itemp_mem = 0
    itemp_grp = addone(itemp_grp)
    trcid = jj+(centurybgc_vars%som1-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM1C'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp ,         &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%som1-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM1N'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))
    !------------------------------------------------------------------------------------
    !new group
    itemp_mem=0
    itemp_grp = addone(itemp_grp)
    trcid = jj+(centurybgc_vars%som2-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM2C'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp ,         &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%som2-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM2N'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))
    !------------------------------------------------------------------------------------
    !new group
    itemp_mem=0
    itemp_grp = addone(itemp_grp)
    trcid = jj+(centurybgc_vars%som3-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM3C'            ,             &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%som3-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM3N'             ,            &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,          &
         trc_group_mem= addone(itemp_mem))

  end subroutine Init_betrbgc

  !-------------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, dz_top, betrtracer_vars, &
       waterflux_vars, tracerboundarycond_vars)
    !
    ! !DESCRIPTION:
    ! set up boundary conditions for tracer movement
    !
    ! !USES:
    use TracerBoundaryCondType, only : tracerboundarycond_type
    use abortutils            , only : endrun
    use shr_log_mod           , only : errMsg => shr_log_errMsg
    use BeTRTracerType        , only : betrtracer_type
    use WaterfluxType         , only : waterflux_type

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_type), intent(in)    :: this
    type(bounds_type)               , intent(in)    :: bounds
    integer                         , intent(in)    :: num_soilc                  ! number of columns in column filter
    integer                         , intent(in)    :: filter_soilc(:)            ! column filter
    type(betrtracer_type)           , intent(in)    :: betrtracer_vars
    real(r8)                        , intent(in)    :: dz_top(bounds%begc: )
    type(waterflux_type)            , intent(in)    :: waterflux_vars
    type(tracerboundarycond_type)   , intent(inout) :: tracerboundarycond_vars

    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'set_boundary_conditions'
    integer :: fc, c

    SHR_ASSERT_ALL((ubound(dz_top) == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))

    associate(                               &
         groupid  => betrtracer_vars%groupid &
         )

      do fc = 1, num_soilc
         c = filter_soilc(fc)

         !values below will be updated with datastream
         !eventually, the following code will be implemented using polymorphism
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%ngwmobile_tracers+1:betrtracer_vars%ntracers) =0._r8                        !zero incoming flux
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_n2)                                    =32.8_r8                      !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_o2)                                    =8.78_r8                      !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ar)                                    =0.3924_r8                    !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_co2x)                                  =0.0168_r8                    !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_ch4)                                   =6.939e-5_r8                  !mol m-3, contant boundary condition
         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_n2o)                                   =1.195e-5_r8                  !mol m-3, contant boundary condition

         tracerboundarycond_vars%tracer_gwdif_concflux_top_col(c,1:2,betrtracer_vars%id_trc_no3x)                                  = 0._r8
         tracerboundarycond_vars%bot_concflux_col(c,1,:)                                                                           = 0._r8                       !zero flux boundary condition
                                                                                                                                                                 !those will be updated with snow resistance and hydraulic wicking resistance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_n2))                                            = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_o2))                                            = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_ar))                                            = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_co2x))                                          = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_ch4))                                           = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
         tracerboundarycond_vars%condc_toplay_col(c,groupid(betrtracer_vars%id_trc_n2o))                                           = 2._r8*1.267e-5_r8/dz_top(c) !m/s surface conductance
      enddo
    end associate
  end subroutine set_boundary_conditions

  !-------------------------------------------------------------------------------
  subroutine calc_bgc_reaction(this, bounds, lbj, ubj, num_soilc, filter_soilc, num_soilp, filter_soilp, jtops, dtime, &
       betrtracer_vars, tracercoeff_vars, waterstate_vars, temperature_vars, soilstate_vars, chemstate_vars,           &
       cnstate_vars, carbonstate_vars, carbonflux_vars, nitrogenstate_vars, nitrogenflux_vars, tracerstate_vars,       &
       tracerflux_vars, plantsoilnutrientflux_vars)
    !
    ! !DESCRIPTION:
    ! do bgc reaction
    ! this returns net carbon fluxes from decay and translocation
    ! and also update the related carbon/nitrogen/phosphorus(potentially) pools of OM
    ! note it is assumed the stoichiometry of the om pools are not changed during decomposition
    !
    ! !USES:
    !
    use tracerfluxType           , only : tracerflux_type
    use tracerstatetype          , only : tracerstate_type
    use tracercoeffType          , only : tracercoeff_type
    use BetrTracerType           , only : betrtracer_type
    use WaterStateType           , only : Waterstate_Type
    use TemperatureType          , only : temperature_type
    use ChemStateType            , only : chemstate_type
    use SoilStatetype            , only : soilstate_type
    use ODEMod                   , only : ode_ebbks1
    use CNStateType              , only : cnstate_type
    use PlantSoilnutrientFluxType, only : plantsoilnutrientflux_type
    use CNVerticalProfileMod     , only : decomp_vertprofiles
    use CNCarbonStateType        , only : carbonstate_type
    use CNCarbonFluxType         , only : carbonflux_type
    use CNNitrogenFluxType       , only : nitrogenflux_type
    use CNNitrogenStateType      , only : nitrogenstate_type
    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_type) , intent(in)    :: this
    type(bounds_type)                , intent(in)    :: bounds                     ! bounds
    integer                          , intent(in)    :: num_soilc                  ! number of columns in column filter
    integer                          , intent(in)    :: filter_soilc(:)            ! column filter
    integer                          , intent(in)    :: num_soilp
    integer                          , intent(in)    :: filter_soilp(:)            ! pft filter
    integer                          , intent(in)    :: jtops(bounds%begc: )       ! top index of each column
    integer                          , intent(in)    :: lbj, ubj                   ! lower and upper bounds, make sure they are > 0
    real(r8)                         , intent(in)    :: dtime                      ! model time step
    type(Waterstate_Type)            , intent(in)    :: waterstate_vars            ! water state variables
    type(temperature_type)           , intent(in)    :: temperature_vars           ! energy state variable
    type(soilstate_type)             , intent(in)    :: soilstate_vars
    type(chemstate_type)             , intent(in)    :: chemstate_vars
    type(betrtracer_type)            , intent(in)    :: betrtracer_vars            ! betr configuration information
    type(tracercoeff_type)           , intent(in)    :: tracercoeff_vars
    type(carbonstate_type)           , intent(in)    :: carbonstate_vars
    type(cnstate_type)               , intent(inout) :: cnstate_vars
    type(carbonflux_type)            , intent(inout) :: carbonflux_vars
    type(nitrogenstate_type)         , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)          , intent(inout) :: nitrogenflux_vars
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    type(tracerflux_type)            , intent(inout) :: tracerflux_vars
    type(plantsoilnutrientflux_type) , intent(inout) :: plantsoilnutrientflux_vars !
    !
    ! !LOCAL VARIABLES:
    integer  :: fc, c, j, k
    real(r8) :: time
    real(r8) :: y0(centurybgc_vars%nstvars, bounds%begc:bounds%endc, lbj:ubj)
    real(r8) :: yf(centurybgc_vars%nstvars, bounds%begc:bounds%endc, lbj:ubj)
    real(r8) :: cn_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    real(r8) :: cp_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    real(r8) :: k_decay(centurybgc_vars%nreactions, bounds%begc:bounds%endc, lbj:ubj)
    real(r8) :: pot_decay_rates(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj) ![mol C/m3/s] potential decay rates for different om pools without nutrient limitation
    real(r8) :: pot_co2_hr(bounds%begc:bounds%endc, lbj:ubj)                                 ![mol C/m3/s], potential co2 respiration rate
    real(r8) :: pot_nh3_immob(bounds%begc:bounds%endc,lbj:ubj)
    real(r8) :: anaerobic_frac(bounds%begc:bounds%endc,lbj:ubj)
    real(r8) :: n2_n2o_ratio_denit(bounds%begc:bounds%endc, lbj:ubj)
    real(r8) :: nh4_no3_ratio(bounds%begc:bounds%endc, lbj:ubj)
    real(r8) :: nuptake_prof(bounds%begc:bounds%endc,1:ubj)
    real(r8) :: pscal
    character(len=32), parameter :: subname ='calc_bgc_reaction'

    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))

    call Extra_inst%Init_Allocate(centurybgc_vars%nom_pools, &
         centurybgc_vars%nreactions, centurybgc_vars%nprimvars)

    call set_reaction_order(centurybgc_vars%nreactions, &
         centurybgc_vars, Extra_inst%is_zero_order)

    !initialize local variables
    y0(:, :, :) = spval
    yf(:, :, :) = spval
    cn_ratios(:,:,:) = nan
    cp_ratios(:,:,:) = nan

    !initialize the state vector
    call init_state_vector(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, &
         centurybgc_vars%nstvars, tracerstate_vars, betrtracer_vars, centurybgc_vars, y0)

    !update the initial vector from external input
    !calculate elemental stoichiometry for different om pools and add mineral nutrient input from other than decaying process

    call bgcstate_ext_update_bfdecomp(bounds, 1, ubj, num_soilc, filter_soilc, carbonflux_vars, nitrogenflux_vars, &
         centurybgc_vars, betrtracer_vars, tracerflux_vars, y0, cn_ratios, cp_ratios)

    !calculate nitrogen uptake profile
    call calc_nuptake_prof(bounds, ubj, num_soilc, filter_soilc,                                                &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, 1:ubj, betrtracer_vars%id_trc_nh3x),  &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, 1:ubj, betrtracer_vars%id_trc_no3x),  &
         col%dz(bounds%begc:bounds%endc,1:ubj), cnstate_vars%nfixation_prof_col(bounds%begc:bounds%endc,1:ubj), &
         nuptake_prof(bounds%begc:bounds%endc,1:ubj))

    !update plant nitrogen uptake potential

    call plantsoilnutrientflux_vars%calc_nutrient_uptake_potential(bounds, num_soilc, filter_soilc, num_soilp, &
         filter_soilp, carbonstate_vars%frootc_patch)

    !calculate multiplicative scalars for decay parameters
    call calc_decompK_multiply_scalar(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, &
         waterstate_vars%finundated_col(bounds%begc:bounds%endc), col%z(bounds%begc:bounds%endc, lbj:ubj),&
         temperature_vars%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj), &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_o2), &
         tracercoeff_vars%aqu2bulkcef_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_o2), &
         soilstate_vars, centurybgc_vars, carbonflux_vars)

    !calculate decay coefficients
    call calc_som_deacyK(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%nom_pools, tracercoeff_vars, tracerstate_vars, &
         betrtracer_vars, centurybgc_vars, carbonflux_vars,dtime, k_decay(1:centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj))

    !calculate potential decay rates, without nutrient constraint
    call calc_sompool_decay(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars, &
         k_decay(1:centurybgc_vars%nom_pools,  bounds%begc:bounds%endc, lbj:ubj), y0(1:centurybgc_vars%nom_totelms, bounds%begc:bounds%endc, lbj:ubj), &
         pot_decay_rates)

    !calculate potential respiration rates by summarizing all om decomposition pathways
    call calc_potential_aerobic_hr(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, cn_ratios, cp_ratios, centurybgc_vars, pot_decay_rates, &
         soilstate_vars%cellsand_col(bounds%begc:bounds%endc,lbj:ubj), pot_co2_hr, pot_nh3_immob)

    !calculate fraction of anerobic environment
    call calc_anaerobic_frac(bounds, lbj, ubj, num_soilc, filter_soilc, jtops,                                 &
         temperature_vars%t_soisno_col(bounds%begc:bounds%endc,lbj:ubj),                                       &
         soilstate_vars, waterstate_vars%h2osoi_vol_col(bounds%begc:bounds%endc,lbj:ubj), pot_co2_hr,          &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_o2), &
         anaerobic_frac(bounds%begc:bounds%endc, lbj:ubj))
    
    !calculate normalized rate for nitrification and denitrification
    call calc_nitrif_denitrif_rate(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, col%dz(bounds%begc:bounds%endc, lbj:ubj), &
         temperature_vars%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj),                                                      &
         chemstate_vars%soil_pH(bounds%begc:bounds%endc, lbj:ubj),  pot_co2_hr, anaerobic_frac,                                &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_nh3x),               &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_no3x),               &
         soilstate_vars, waterstate_vars, carbonflux_vars, n2_n2o_ratio_denit, nh4_no3_ratio,                                  &
         k_decay(centurybgc_vars%lid_nh4_nit_reac, bounds%begc:bounds%endc, lbj:ubj),                                          &
         k_decay(centurybgc_vars%lid_no3_den_reac, bounds%begc:bounds%endc, lbj:ubj))

    !now there is no plant nitrogen uptake, I tend to create a new structure to indicate plant nutrient demand when it is hooked
    !back with CLM

    call calc_plant_nitrogen_uptake_prof(bounds, ubj, num_soilc, filter_soilc, col%dz(bounds%begc:bounds%endc, lbj:ubj), &
         plantsoilnutrientflux_vars%plant_minn_uptake_potential_col(bounds%begc:bounds%endc),                            &
         nuptake_prof(bounds%begc:bounds%endc,1:ubj),                                                                    &
         k_decay(centurybgc_vars%lid_plant_minn_up_reac, bounds%begc:bounds%endc ,1:ubj))

    !apply root distribution here
    call apply_plant_root_respiration_prof(bounds, ubj, num_soilc, filter_soilc, &
         carbonflux_vars%rr_col(bounds%begc:bounds%endc),                        &
         cnstate_vars%nfixation_prof_col(bounds%begc:bounds%endc,1:ubj),         &
         k_decay(centurybgc_vars%lid_at_rt_reac, bounds%begc:bounds%endc, 1:ubj))

    !do ode integration and update state variables for each layer
    !lpr = .true.
    do j = lbj, ubj
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          if(j<jtops(c))cycle
          !assign parameters for stoichiometric matrix calculation
          call Extra_inst%AAssign(cn_r=cn_ratios(:,c,j),cp_r=cp_ratios(:,c,j), &
               k_d=k_decay(:,c,j),                                             &
               n2_n2o_r_denit=n2_n2o_ratio_denit(c,j),                         &
               cell_sand=soilstate_vars%cellsand_col(c,j),                     &
               betrtracer_vars=betrtracer_vars,                                &
               gas2bulkcef=tracercoeff_vars%gas2bulkcef_mobile_col(c,j,:),     &
               aere_cond=tracercoeff_vars%aere_cond_col(c,:),                  &
               tracer_conc_atm=tracerstate_vars%tracer_conc_atm_col(c,:))
          !update state variables
          time = 0._r8

          yf(:,c,j)=y0(:,c,j) !this will allow to turn off the bgc reaction for debugging purpose

          call ode_ebbks1(one_box_century_bgc, y0(:,c,j),         &
               centurybgc_vars%nprimvars,centurybgc_vars%nstvars, &
               time, dtime, yf(:,c,j), pscal)

          if(pscal<5.e-1_r8)then
             write(iulog,*)'lat, lon=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
             write(iulog,*)'col, lev, pscal=',c, j, pscal
             write(iulog,*)'nstep =',get_nstep()
             call endrun()
          endif

       enddo
    enddo
    call bgcstate_ext_update_afdecomp(bounds, 1, ubj, num_soilc, filter_soilc, &
         carbonflux_vars, nitrogenflux_vars,                                   &
         centurybgc_vars, betrtracer_vars, tracerflux_vars, yf)

    !retrieve the flux variable
    call retrieve_flux_vars(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, &
         centurybgc_vars%nstvars, dtime, yf, y0,                              &
         centurybgc_vars, betrtracer_vars, tracerflux_vars,                   &
         carbonflux_vars, nitrogenflux_vars, plantsoilnutrientflux_vars)


    !retrieve the state variable, state variable will be updated later
    call retrieve_state_vars(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, &
         centurybgc_vars%nstvars,  yf, centurybgc_vars, betrtracer_vars, tracerstate_vars)

    call Extra_inst%DDeallocate()

  end subroutine calc_bgc_reaction


  !-------------------------------------------------------------------------------
  subroutine do_tracer_equilibration(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       betrtracer_vars, tracercoeff_vars, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! equilibrate tracers that has solid and mobile phases
    ! using the theory of mass action. When the redox-ladder is on, this
    ! subroutine will update the change of pH due to tracer transport, or
    ! !USES:
    !
    use tracerstatetype       , only : tracerstate_type
    use tracercoeffType       , only : tracercoeff_type
    use BeTRTracerType        , only : betrtracer_type

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_type), intent(in)    :: this
    type(bounds_type),                intent(in)    :: bounds
    integer,                          intent(in)    :: lbj, ubj
    integer,                          intent(in)    :: jtops(bounds%begc: )        ! top label of each column
    integer,                          intent(in)    :: num_soilc
    integer,                          intent(in)    :: filter_soilc(:)
    type(betrtracer_type),            intent(in)    :: betrtracer_vars
    type(tracercoeff_type),           intent(in)    :: tracercoeff_vars
    type(tracerstate_type),           intent(inout) :: tracerstate_vars
    ! !LOCAL VARIABLES:
    character(len=32) :: subname = 'do_tracer_equilibration'

    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))
    
    !depending on the simulation type, an implementation of aqueous chemistry will be
    !employed to separate out the adsorbed phase
    !It should be noted that this formulation excludes the use of linear isotherm, which
    !can be integrated through the retardation factor    
    
  end subroutine do_tracer_equilibration

  !-----------------------------------------------------------------------
  subroutine readParams(this, ncid, betrtracer_vars )
    !
    ! !DESCRIPTION:
    ! read in model parameters
    !
    ! !USES:
    use BeTRTracerType   , only : BeTRTracer_Type
    use ncdio_pio        , only : file_desc_t
    use BGCCenturyParMod , only : readCentDecompBgcParams
    use BGCCenturyParMod , only : readCentNitrifDenitrifParams
    use BGCCenturyParMod , only : readCentCNAllocParams
    !
    class(bgc_reaction_CENTURY_type) , intent(in)    :: this
    type(BeTRTracer_Type)            , intent(inout) :: betrtracer_vars
    type(file_desc_t)                , intent(inout) :: ncid  ! pio netCDF file id

    call readCentDecompBgcParams (ncid, centurybgc_vars%nelms, betrtracer_vars )

    call readCentNitrifDenitrifParams ( ncid )

    call readCentCNAllocParams(ncid)

  end subroutine readParams

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, betrtracer_vars,  waterstate_vars, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! cold initialization
    ! !USES:
    !
    use BeTRTracerType  , only : BeTRTracer_Type
    use tracerstatetype , only : tracerstate_type
    use WaterstateType  , only : waterstate_type
    use LandunitType    , only : lun
    use ColumnType      , only : col
    use PatchType       , only : pft
    use clm_varcon      , only : spval, ispval

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_type) , intent(in)    :: this
    type(bounds_type)                , intent(in)    :: bounds
    type(BeTRTracer_Type)            , intent(in)    :: betrtracer_vars
    type(waterstate_type)            , intent(in)    :: waterstate_vars
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l, k, j
    integer :: fc                                        ! filter index
    integer :: begc, endc
    integer :: begg, endg

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg

    associate(                                    &
         volatileid => betrtracer_vars%volatileid &
         )

      do c = bounds%begc, bounds%endc
         l = col%landunit(c)
         if (lun%ifspecial(l)) then
            if(betrtracer_vars%ngwmobile_tracers>0)then
               tracerstate_vars%tracer_conc_mobile_col(c,:,:)        = spval
               tracerstate_vars%tracer_conc_surfwater_col(c,:)       = spval
               tracerstate_vars%tracer_conc_aquifer_col(c,:)         = spval
               tracerstate_vars%tracer_conc_grndwater_col(c,:)       = spval
               tracerstate_vars%tracer_conc_atm_col(c,:)             = spval
            endif
            if(betrtracer_vars%ntracers > betrtracer_vars%ngwmobile_tracers)then
               tracerstate_vars%tracer_conc_solid_passive_col(c,:,:) = spval
            endif
            if(betrtracer_vars%nsolid_equil_tracers>0)then
               tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = spval
            endif
         endif
         tracerstate_vars%tracer_soi_molarmass_col(c,:)            = spval

         if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
            !dual phase tracers

             tracerstate_vars%tracer_conc_mobile_col    (c,:, :                                    )  = 0._r8
             tracerstate_vars%tracer_conc_surfwater_col (c,:                                       )  = 0._r8
             tracerstate_vars%tracer_conc_aquifer_col   (c,:                                       )  = 0._r8
             tracerstate_vars%tracer_conc_grndwater_col (c,:                                       )  = 0._r8
             tracerstate_vars%tracer_conc_atm_col       (c,volatileid (betrtracer_vars%id_trc_n2   )) = 32.8_r8
             tracerstate_vars%tracer_conc_atm_col       (c,volatileid (betrtracer_vars%id_trc_o2   )) = 8.78_r8
             tracerstate_vars%tracer_conc_atm_col       (c,volatileid (betrtracer_vars%id_trc_ar   )) = 0.3924_r8
             tracerstate_vars%tracer_conc_atm_col       (c,volatileid (betrtracer_vars%id_trc_co2x )) = 0.0168_r8
             tracerstate_vars%tracer_conc_atm_col       (c,volatileid (betrtracer_vars%id_trc_ch4  )) = 6.939e-5_r8
             tracerstate_vars%tracer_conc_atm_col       (c,volatileid (betrtracer_vars%id_trc_n2o  )) = 1.195e-5_r8

            !solid tracers
            if(betrtracer_vars%ngwmobile_tracers < betrtracer_vars%ntracers)then
               tracerstate_vars%tracer_conc_solid_passive_col(c,:,:) = 0._r8
            endif

            if(betrtracer_vars%nsolid_equil_tracers>0)then
               tracerstate_vars%tracer_conc_solid_equil_col(c, :, :) = 0._r8
            endif
            tracerstate_vars%tracer_soi_molarmass_col(c,:)          = 0._r8
         endif
      enddo
    end associate
  end subroutine InitCold

  !--------------------------------------------------------------------------------------------------------------------
  subroutine betr_alm_flux_statevar_feedback(this, bounds, num_soilc, filter_soilc, &
       carbonstate_vars, nitrogenstate_vars, nitrogenflux_vars, tracerstate_vars,   &
       tracerflux_vars,  betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! do state and flux variable exchange between betr and alm
    !
    ! !USES:
    use shr_kind_mod             , only : r8 => shr_kind_r8
    use tracerfluxType           , only : tracerflux_type
    use tracerstatetype          , only : tracerstate_type
    use decompMod                , only : bounds_type
    use BeTRTracerType           , only : BeTRTracer_Type
    use CNCarbonStateType        , only : carbonstate_type
    use CNNitrogenStateType      , only : nitrogenstate_type
    use CNNitrogenFluxType       , only : nitrogenflux_type


    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_type) , intent(in)    :: this
    type(bounds_type)                , intent(in)    :: bounds             ! bounds
    integer                          , intent(in)    :: num_soilc          ! number of columns in column filter
    integer                          , intent(in)    :: filter_soilc(:)    ! column filter
    type(betrtracer_type)            , intent(in)    :: betrtracer_vars    ! betr configuration information
    type(tracerstate_type)           , intent(in)    :: tracerstate_vars
    type(tracerflux_type)            , intent(in)    :: tracerflux_vars
    type(carbonstate_type)           , intent(inout) :: carbonstate_vars
    type(nitrogenflux_type)          , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)         , intent(inout) :: nitrogenstate_vars !

    call assign_nitrogen_hydroloss(bounds, num_soilc, filter_soilc, tracerflux_vars, &
         nitrogenflux_vars, betrtracer_vars)

    call assign_OM_CNpools(bounds, num_soilc, filter_soilc,  carbonstate_vars, &
         nitrogenstate_vars, tracerstate_vars, betrtracer_vars, centurybgc_vars)

  end subroutine betr_alm_flux_statevar_feedback

  !---------------------------------------------------------------
  subroutine init_betr_alm_bgc_coupler(this, bounds, carbonstate_vars, &
       nitrogenstate_vars, betrtracer_vars, tracerstate_vars)

    !
    ! !DESCRIPTION:
    ! state variable exchange between betr and alm
    ! !USES:
    use clm_varcon               , only : natomw, catomw
    use clm_varpar               , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
    use CNCarbonStateType        , only : carbonstate_type
    use CNNitrogenStateType      , only : nitrogenstate_type
    use tracerstatetype          , only : tracerstate_type
    use BetrTracerType           , only : betrtracer_type
    use clm_varpar               , only : nlevtrc_soil
    use landunit_varcon          , only : istsoil, istcrop

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_type) , intent(in)    :: this
    type(bounds_type)                , intent(in)    :: bounds
    type(tracerstate_type)           , intent(inout) :: tracerstate_vars
    type(betrtracer_type)            , intent(in)    :: betrtracer_vars   ! betr configuration information
    type(carbonstate_type)           , intent(in)    :: carbonstate_vars
    type(nitrogenstate_type)         , intent(in)    :: nitrogenstate_vars

    ! !LOCAL VARIABLES:
    integer, parameter :: i_soil1 = 5
    integer, parameter :: i_soil2 = 6
    integer, parameter :: i_soil3 = 7
    integer            :: c, j, k, l
    character(len=255)   :: subname = 'init_betr_alm_bgc_coupler'

    associate(                                                                         &
         id_trc_no3x               => betrtracer_vars%id_trc_no3x                    , &
         id_trc_nh3x               => betrtracer_vars%id_trc_nh3x                    , &
         decomp_cpools_vr          => carbonstate_vars%decomp_cpools_vr_col          , &
         decomp_npools_vr          => nitrogenstate_vars%decomp_npools_vr_col        , &
         smin_no3_vr_col           => nitrogenstate_vars%smin_no3_vr_col             , &
         smin_nh4_vr_col           => nitrogenstate_vars%smin_nh4_vr_col             , &
         tracer_conc_mobile        => tracerstate_vars%tracer_conc_mobile_col        , &
         tracer_conc_solid_passive => tracerstate_vars%tracer_conc_solid_passive_col , &
         c_loc                     => centurybgc_vars%c_loc                          , &
         n_loc                     => centurybgc_vars%n_loc                          , &
         lit1                      => centurybgc_vars%lit1                           , &
         lit2                      => centurybgc_vars%lit2                           , &
         lit3                      => centurybgc_vars%lit3                           , &
         som1                      => centurybgc_vars%som1                           , &
         som2                      => centurybgc_vars%som2                           , &
         som3                      => centurybgc_vars%som3                           , &
         cwd                       => centurybgc_vars%cwd                            , &
         nelms                     => centurybgc_vars%nelms                            &
         )

      !initialize tracer based on carbon/nitrogen pools
      !eventually, this will replace the century bgc
      do j = 1, nlevtrc_soil
         do c = bounds%begc, bounds%endc
            l = col%landunit(c)
            if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
               tracer_conc_mobile(c,j,id_trc_no3x)=smin_no3_vr_col(c,j) /natomw
               tracer_conc_mobile(c,j,id_trc_nh3x)=smin_nh4_vr_col(c,j) /natomw
               k = lit1; tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) = decomp_cpools_vr(c,j,i_met_lit) / catomw
               k = lit2; tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) = decomp_cpools_vr(c,j,i_cel_lit) / catomw
               k = lit3; tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) = decomp_cpools_vr(c,j,i_lig_lit) / catomw
               k = cwd ; tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) = decomp_cpools_vr(c,j,i_cwd    ) / catomw
               k = som1; tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) = decomp_cpools_vr(c,j,i_soil1  ) / catomw
               k = som2; tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) = decomp_cpools_vr(c,j,i_soil2  ) / catomw
               k = som3; tracer_conc_solid_passive(c,j,(k-1)*nelms+c_loc) = decomp_cpools_vr(c,j,i_soil3  ) / catomw

               k = lit1; tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) = decomp_npools_vr(c,j,i_met_lit) / natomw
               k = lit2; tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) = decomp_npools_vr(c,j,i_cel_lit) / natomw
               k = lit3; tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) = decomp_npools_vr(c,j,i_lig_lit) / natomw
               k = cwd ; tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) = decomp_npools_vr(c,j,i_cwd    ) / natomw
               k = som1; tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) = decomp_npools_vr(c,j,i_soil1  ) / natomw
               k = som2; tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) = decomp_npools_vr(c,j,i_soil2  ) / natomw
               k = som3; tracer_conc_solid_passive(c,j,(k-1)*nelms+n_loc) = decomp_npools_vr(c,j,i_soil3  ) / natomw
            endif
         enddo
      enddo
    end associate
  end subroutine init_betr_alm_bgc_coupler

  !-------------------------------------------------------------------------------
  
  
  subroutine one_box_century_bgc(ystate, dtime, time, nprimvars, nstvars, dydt)
    !
    ! !DESCRIPTION:
    ! do single box bgc
    !
    !the equations to be solved are in the form
    !
    ! dx/dt=I+A*R, where I is the input, A is the stoichiometric matrix, and R is the reaction vector
    ! the input only contains litter input and mineral nutrient, som is assumed to be of fixed stoichiometry
    ! !USES:
    use SOMStateVarUpdateMod  , only : calc_dtrend_som_bgc
    use BGCCenturySubMod      , only : calc_cascade_matrix
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in)  :: nstvars
    integer,  intent(in)  :: nprimvars
    real(r8), intent(in)  :: dtime
    real(r8), intent(in)  :: time
    real(r8), intent(in)  :: ystate(nstvars)
    real(r8), intent(out) :: dydt(nstvars)

    ! !LOCAL VARIABLES:
    integer  :: lk, jj
    real(r8) :: cascade_matrix(nstvars, Extra_inst%nr)
    logical  :: nitrogen_limit_flag(centurybgc_vars%nom_pools)
    real(r8) :: reaction_rates(Extra_inst%nr)
    real(r8) :: o2_consump, o2_limit

    !calculate cascade matrix, which contains the stoichiometry for all reactions
    call calc_cascade_matrix(nstvars, Extra_inst%nr, Extra_inst%cn_ratios, Extra_inst%cp_ratios, &
         Extra_inst%n2_n2o_ratio_denit, Extra_inst%cellsand, centurybgc_vars, nitrogen_limit_flag, cascade_matrix)

    !do pool degradation
    do lk = 1, Extra_inst%nr
       if(Extra_inst%is_zero_order(lk))then

          if ( spinup_state .eq. 1 ) then
             !spinup stage
             if(lk == centurybgc_vars%lid_o2_aere_reac)then
                jj = centurybgc_vars%lid_o2
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                !I add the following line to disconnect the nitrogen and oxygen interaction
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)
             else
                reaction_rates(lk) = Extra_inst%k_decay(lk)            !this effective defines the plant nitrogen demand
             endif
          else
             ! normal run stage
             if(lk == centurybgc_vars%lid_o2_aere_reac)then
                jj = centurybgc_vars%lid_o2
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)
             elseif(lk == centurybgc_vars%lid_ch4_aere_reac)then
                jj = centurybgc_vars%lid_ch4
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)

             elseif(lk == centurybgc_vars%lid_ar_aere_reac)then
                jj = centurybgc_vars%lid_ar
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)

             elseif(lk == centurybgc_vars%lid_n2_aere_reac)then
                jj = centurybgc_vars%lid_n2
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)

             elseif(lk == centurybgc_vars%lid_co2_aere_reac)then
                jj = centurybgc_vars%lid_co2
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)

             elseif(lk == centurybgc_vars%lid_n2o_aere_reac)then
                jj = centurybgc_vars%lid_n2o
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)

             else
                reaction_rates(lk) = Extra_inst%k_decay(lk)            !this effective defines the plant nitrogen demand
             endif
          endif
       else
          reaction_rates(lk)=ystate(centurybgc_vars%primvarid(lk))*Extra_inst%k_decay(lk)
       endif
    enddo

    !obtain total oxygen consumption rate
    o2_consump = DOT_PRODUCT(cascade_matrix(centurybgc_vars%lid_o2,1:Extra_inst%nr),reaction_rates(1:Extra_inst%nr))

    if(-o2_consump*dtime > ystate(centurybgc_vars%lid_o2))then
       o2_limit=-ystate(centurybgc_vars%lid_o2)/(o2_consump*dtime)
       do lk = 1, Extra_inst%nr
          if(centurybgc_vars%is_aerobic_reac(lk))then
             reaction_rates(lk) = reaction_rates(lk)*o2_limit
          endif
       enddo
    endif

    call apply_nutrient_down_regulation(nstvars, Extra_inst%nr, nitrogen_limit_flag, &
         ystate(centurybgc_vars%lid_nh4), ystate(centurybgc_vars%lid_no3),           &
         dtime, cascade_matrix, reaction_rates)

    call calc_dtrend_som_bgc(nstvars, Extra_inst%nr, &
         cascade_matrix(1:nstvars, 1:Extra_inst%nr), &
         reaction_rates(1:Extra_inst%nr), dydt)

  end subroutine one_box_century_bgc

  !-------------------------------------------------------------------------------
  subroutine apply_nutrient_down_regulation(nstvars, nreactions, nitrogen_limit_flag, &
       smin_nh4, smin_no3, dtime, cascade_matrix, reaction_rates)

    ! !DESCRIPTION:
    !
    ! this down-regulation considers nitrogen made available from gross mineralization
    ! this implements is corresponding to the CLM-2 approach as described in Tang and Riley (2015), BG, tehcnique note.
    ! !USES:
    use clm_varctl,   only : CNAllocate_Carbon_only
    use MathfuncMod,  only : safe_div

    ! !ARGUMENTS:
    integer , intent(in)    :: nstvars
    integer , intent(in)    :: nreactions
    logical , intent(in)    :: nitrogen_limit_flag(centurybgc_vars%nom_pools)
    real(r8), intent(in)    :: smin_nh4
    real(r8), intent(in)    :: smin_no3
    real(r8), intent(in)    :: dtime
    real(r8), intent(inout) :: cascade_matrix(nstvars, nreactions)
    real(r8), intent(inout) :: reaction_rates(nreactions)
    ! !LOCAL VARIABLES:
    real(r8)                :: decomp_plant_minn_demand_flx
    real(r8)                :: tot_nh4_demand_flx
    real(r8)                :: tot_no3_demand_flx
    real(r8)                :: decomp_plant_residual_minn_demand_flx
    real(r8)                :: smin_nh4_to_decomp_plant_flx
    real(r8)                :: smin_no3_to_decomp_plant_flx
    real(r8)                :: tot_sminn_to_decomp_plant_flx
    real(r8)                :: frac_nh4_to_decomp_plant
    real(r8)                :: supp_nh4_to_decomp_plant_flx
    real(r8)                :: frac_supp_nh4_to_decomp_plant
    real(r8)                :: gross_min_nh4_flx
    real(r8)                :: alpha
    real(r8)                :: frac_gross_immob=1.0_r8
    integer                 :: reac

    associate(                                                            & !
         nom_pools             => centurybgc_vars%nom_pools             , & !
         lid_nh4               => centurybgc_vars%lid_nh4               , & !
         lid_no3               => centurybgc_vars%lid_no3               , & !
         lid_plant_minn        => centurybgc_vars%lid_plant_minn        , & !
         lid_minn_nh4_immob    => centurybgc_vars%lid_minn_nh4_immob    , & !
         lid_minn_no3_immob    => centurybgc_vars%lid_minn_no3_immob    , & !
         lid_minn_nh4_plant    => centurybgc_vars%lid_minn_nh4_plant    , & !
         lid_minn_no3_plant    => centurybgc_vars%lid_minn_no3_plant    , & !
         lid_nh4_supp          => centurybgc_vars%lid_nh4_supp          , & !
         lid_nh4_nit           => centurybgc_vars%lid_nh4_nit           , & !
         lid_plant_minn_up_reac=> centurybgc_vars%lid_plant_minn_up_reac, & !
         lid_nh4_nit_reac      => centurybgc_vars%lid_nh4_nit_reac      , & !
         lid_no3_den_reac      => centurybgc_vars%lid_no3_den_reac        & !
         )

      decomp_plant_minn_demand_flx = 0._r8
      gross_min_nh4_flx = 0._r8
      do reac = 1,  nom_pools
         if(nitrogen_limit_flag(reac))then
            decomp_plant_minn_demand_flx = decomp_plant_minn_demand_flx - &
                 reaction_rates(reac) * cascade_matrix(lid_nh4, reac)
         else
            gross_min_nh4_flx = gross_min_nh4_flx + &
                 reaction_rates(reac) * cascade_matrix(lid_nh4, reac)
         endif
      enddo

      !add nitrogen demand from plant
      reac = lid_plant_minn_up_reac
      decomp_plant_minn_demand_flx = decomp_plant_minn_demand_flx - &
           reaction_rates(reac) * cascade_matrix(lid_nh4, reac)

      !in clm-century, nh4 is first competed between decomposer immobilization, plant and nitrification
      !
      reac = lid_nh4_nit_reac
      tot_nh4_demand_flx = decomp_plant_minn_demand_flx -       &
           reaction_rates(reac) * cascade_matrix(lid_nh4 ,reac) &
           - gross_min_nh4_flx*frac_gross_immob

      if(tot_nh4_demand_flx*dtime>smin_nh4)then
         if(CNAllocate_Carbon_only())then

            !nitrifier uses what it is provided
            !plant use the remaining nh4 and request external supply from supp nh4
            if(reaction_rates(reac)<1.e-40_r8)then
               alpha = 0._r8
            else
               alpha = -smin_nh4/(reaction_rates(reac)*cascade_matrix(lid_nh4,reac)*dtime)
            endif
            reaction_rates(reac) = reaction_rates(reac) * min(alpha, 1._r8)

            smin_nh4_to_decomp_plant_flx = smin_nh4/dtime+reaction_rates(reac)*cascade_matrix(lid_nh4,reac)
            decomp_plant_residual_minn_demand_flx = decomp_plant_minn_demand_flx - smin_nh4_to_decomp_plant_flx

         else
            !nitrifiers, decomposers and plants are nh4 limited
            alpha = smin_nh4/(tot_nh4_demand_flx*dtime)

            !downregulate nitrification
            reaction_rates(lid_nh4_nit_reac)      = reaction_rates(lid_nh4_nit_reac)*alpha
            smin_nh4_to_decomp_plant_flx          = smin_nh4/dtime + reaction_rates(lid_nh4_nit_reac) * cascade_matrix(lid_nh4, reac)
            decomp_plant_residual_minn_demand_flx = decomp_plant_minn_demand_flx - smin_nh4_to_decomp_plant_flx
         endif
      else
         !none is nh4 limited
         smin_nh4_to_decomp_plant_flx = decomp_plant_minn_demand_flx
         decomp_plant_residual_minn_demand_flx = 0._r8
      endif
      !avoid negative smin_nh4 due to roundoff
      smin_nh4_to_decomp_plant_flx = max(smin_nh4_to_decomp_plant_flx -1.e-21_r8, &
                                         smin_nh4_to_decomp_plant_flx)

      reac = lid_no3_den_reac
      tot_no3_demand_flx = decomp_plant_residual_minn_demand_flx - &
           reaction_rates(reac) * cascade_matrix(lid_no3 ,reac)

      !then no3 is competed between denitrification and residual request from decomposer immobilization and plant demand
      if(tot_no3_demand_flx * dtime>smin_no3)then

         if(CNAllocate_Carbon_only())then
            !denitrifiers is given what is available
            if(abs(reaction_rates(reac))<1.e-40_r8)then
               alpha = 0._r8
            else
               alpha = -smin_no3/(dtime*reaction_rates(reac)*cascade_matrix(lid_no3,reac))
            endif
            reaction_rates(lid_no3_den_reac ) = reaction_rates(lid_no3_den_reac )*min(alpha,1._r8)
            smin_no3_to_decomp_plant_flx = smin_no3/dtime + reaction_rates(lid_no3_den_reac ) * cascade_matrix(lid_no3 ,reac)
         else
            !denitrifiers, decomposers and plants are no3 limited
            alpha = smin_no3/(tot_no3_demand_flx*dtime)
            reaction_rates(lid_no3_den_reac ) = reaction_rates(lid_no3_den_reac )*alpha

            smin_no3_to_decomp_plant_flx = smin_no3/dtime + reaction_rates(lid_no3_den_reac ) * cascade_matrix(lid_no3 ,reac)
         endif
      else
         smin_no3_to_decomp_plant_flx = tot_no3_demand_flx
      endif

      !avoid negative smin_no3 due to roundoff
      smin_no3_to_decomp_plant_flx  = max(smin_no3_to_decomp_plant_flx-1.e-21_r8,0._r8)
      tot_sminn_to_decomp_plant_flx = smin_nh4_to_decomp_plant_flx + smin_no3_to_decomp_plant_flx

      if(CNAllocate_Carbon_only())then
         supp_nh4_to_decomp_plant_flx  = decomp_plant_minn_demand_flx - tot_sminn_to_decomp_plant_flx
         tot_sminn_to_decomp_plant_flx = decomp_plant_minn_demand_flx
      else
         supp_nh4_to_decomp_plant_flx = 0._r8
      endif

      if(tot_sminn_to_decomp_plant_flx < decomp_plant_minn_demand_flx)then
         !plant & decomp are nitrogen limited
         alpha = tot_sminn_to_decomp_plant_flx/decomp_plant_minn_demand_flx
      else
         alpha = 1._r8
      endif

      if(smin_nh4_to_decomp_plant_flx>=tot_sminn_to_decomp_plant_flx)then
         frac_nh4_to_decomp_plant = 1._r8
      else
         frac_nh4_to_decomp_plant = smin_nh4_to_decomp_plant_flx/tot_sminn_to_decomp_plant_flx

         if(supp_nh4_to_decomp_plant_flx>0._r8)then
            frac_supp_nh4_to_decomp_plant=1._r8-frac_nh4_to_decomp_plant
         else
            frac_supp_nh4_to_decomp_plant = 0._r8
         endif
      endif
      !revise the stoichiometry matix elements
      !for decomposers

      do reac = 1,  nom_pools
         if(nitrogen_limit_flag(reac))then

            reaction_rates(reac) = reaction_rates(reac) * alpha
            cascade_matrix(lid_no3, reac) = cascade_matrix(lid_nh4, reac)*(1._r8-frac_nh4_to_decomp_plant-frac_supp_nh4_to_decomp_plant)

            if(lid_nh4_supp>0)then
               cascade_matrix (lid_nh4_supp       , reac) =  cascade_matrix(lid_nh4, reac) * frac_supp_nh4_to_decomp_plant
               cascade_matrix (lid_nh4            , reac) =  cascade_matrix(lid_nh4, reac) - cascade_matrix(lid_no3, reac) - cascade_matrix(lid_nh4_supp, reac)
               cascade_matrix (lid_minn_nh4_immob , reac) = -cascade_matrix(lid_nh4, reac) - cascade_matrix(lid_nh4_supp, reac)
            else
               cascade_matrix(lid_nh4            , reac) =  cascade_matrix(lid_nh4, reac) - cascade_matrix(lid_no3, reac)
               cascade_matrix(lid_minn_nh4_immob , reac) = -cascade_matrix(lid_nh4, reac)
            endif
            cascade_matrix(lid_minn_no3_immob, reac) = -cascade_matrix(lid_no3, reac)
         endif
      enddo

      !for plant
      reac                          = lid_plant_minn_up_reac
      reaction_rates(reac)          = reaction_rates(reac) * alpha
      cascade_matrix(lid_nh4, reac) = -frac_nh4_to_decomp_plant

      if(lid_nh4_supp>0)then
         cascade_matrix(lid_nh4_supp       , reac) = -frac_supp_nh4_to_decomp_plant
         cascade_matrix(lid_no3            , reac) = -(1._r8 - frac_nh4_to_decomp_plant - frac_supp_nh4_to_decomp_plant)
         cascade_matrix(lid_minn_nh4_plant , reac) = -cascade_matrix(lid_nh4, reac)-cascade_matrix(lid_nh4_supp, reac)
      else
         cascade_matrix(lid_no3            , reac) = -(1._r8-frac_nh4_to_decomp_plant)
         cascade_matrix(lid_minn_nh4_plant , reac) = -cascade_matrix(lid_nh4, reac)
      endif
      cascade_matrix(lid_minn_no3_plant, reac) = -cascade_matrix(lid_no3, reac)
    end associate
  end subroutine apply_nutrient_down_regulation

end module BGCReactionsCenturyType

