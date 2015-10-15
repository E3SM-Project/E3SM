module BGCReactionsCenturyECAType

#include "shr_assert.h"

  !
  ! !DESCRIPTION:
  ! do ECA based nitrogen competition in betr.
  ! this code uses the operator automated down-regulation scheme
  ! HISTORY:
  ! Created by Jinyun Tang, Oct 2nd, 2014
  ! Note: ECA parameters are note tuned.
  !
  ! !USES:
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
  public :: bgc_reaction_CENTURY_ECA_type
  type(centurybgc_type), private :: centurybgc_vars
  logical :: ldebug
  !integer, private :: lpr
  type, extends(bgc_reaction_type) :: &
       bgc_reaction_CENTURY_ECA_type
  private

contains
    procedure :: Init_betrbgc                 ! initialize betr bgc
    procedure :: set_boundary_conditions      ! set top/bottom boundary conditions for various tracers
    procedure :: calc_bgc_reaction            ! doing bgc calculation
    procedure :: init_boundary_condition_type ! initialize type of top boundary conditions
    procedure :: do_tracer_equilibration      ! do equilibrium tracer chemistry
    procedure :: initCold
    procedure :: readParams
    procedure :: init_betr_alm_bgc_coupler    ! update state vars using other bgc parts in alm
    procedure :: betr_alm_flux_statevar_feedback
  end type bgc_reaction_CENTURY_ECA_type

  type, private :: Extra_type
     real(r8), pointer :: cn_ratios(:)           !cn ratio of om pool
     real(r8), pointer :: cp_ratios(:)           !cp ratio of om pool
     real(r8), pointer :: k_decay(:)             !decay parameter for all reactions
     real(r8), pointer :: scal_f(:)              !scaling factor for first order sink
     real(r8), pointer :: conv_f(:)              !converting factor for first order sink
     real(r8), pointer :: conc_f(:)              !external forcing strength
     real(r8)          :: n2_n2o_ratio_denit     !ratio of n2 to n2o during denitrification
     real(r8)          :: pct_sand               !sand content [0-100]
     real(r8)          :: pct_clay               !clay content [0-100]
     real(r8)          :: plant_frts             !fine roots for nutrient uptake
     logical,  pointer :: is_zero_order(:)
     integer           :: nr                     !number of reactions involved
   contains
     procedure, public :: Init_Allocate
     procedure, public :: DDeallocate
     procedure, public :: AAssign
  end type Extra_type
  type(Extra_type), private :: Extra_inst


  interface bgc_reaction_CENTURY_ECA_type
     module procedure constructor

  end interface bgc_reaction_CENTURY_ECA_type

contains

  subroutine Init_Allocate(this, nompools, nreacts, nprimstvars)
    !
    ! !DESCRIPTION:
    ! memory allocation for the data type specified by this
    !
    ! !ARGUMENTS:
    class(Extra_type) :: this

    integer, intent(in) :: nompools
    integer, intent(in) :: nreacts
    integer, intent(in) :: nprimstvars     !number of primary state variables

    allocate(this%cn_ratios(nompools))
    allocate(this%cp_ratios(nompools))
    allocate(this%k_decay(nreacts))
    allocate(this%scal_f(nprimstvars));    this%scal_f(:) = 0._r8
    allocate(this%conv_f(nprimstvars));    this%conv_f(:) = 0._r8
    allocate(this%conc_f(nprimstvars));    this%conc_f(:) = 0._r8
    allocate(this%is_zero_order(nreacts)); this%is_zero_order(:) = .false.
    this%nr = nreacts

  end subroutine Init_Allocate

  !-------------------------------------------------------------------------------

  subroutine DDeallocate(this)
    !
    ! !DESCRIPTION:
    ! deallocate memory for the data type specified by this
    !
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

  subroutine AAssign(this, cn_r,cp_r, k_d,  n2_n2o_r_denit, cell_sand, cell_clay, &
       plant_froots, betrtracer_vars, gas2bulkcef, aere_cond, tracer_conc_atm)
    !
    ! !DESCRIPTION:
    ! assign member values for the data type specified by this
    ! !USES:
    use BeTRTracerType              , only : betrtracer_type
    ! !ARGUMENTS:
    class(Extra_type) :: this
    real(r8), dimension(:), intent(in) :: cn_r
    real(r8), dimension(:), intent(in) :: cp_r
    real(r8), dimension(:), intent(in) :: k_d
    real(r8)              , intent(in) :: n2_n2o_r_denit
    real(r8)              , intent(in) :: cell_sand
    real(r8)              , intent(in) :: cell_clay
    real(r8)              , intent(in) :: plant_froots
    type(BeTRtracer_type ), intent(in) :: betrtracer_vars
    real(r8)              , intent(in) :: gas2bulkcef(1:betrtracer_vars%nvolatile_tracers)
    real(r8)              , intent(in) :: aere_cond(1:betrtracer_vars%nvolatile_tracers)
    real(r8)              , intent(in) :: tracer_conc_atm(1:betrtracer_vars%nvolatile_tracers)

    ! !LOCAL VARIABLES:
    integer :: n1, n2, n3, j

    n1 = size(cn_r)
    n2 = size(cp_r)
    n3 = size(k_d)
    SHR_ASSERT_ALL((n1              == n2),        errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((n3              == this%nr),   errMsg(__FILE__,__LINE__))
    this%cn_ratios(1:n1) = cn_r
    this%cp_ratios(1:n2) = cp_r

    this%n2_n2o_ratio_denit = n2_n2o_r_denit
    this%pct_sand           = cell_sand
    this%pct_clay           = cell_clay
    this%k_decay            = k_d
    this%plant_frts         = plant_froots

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
  type(bgc_reaction_CENTURY_ECA_type) function constructor()
    !
    ! ! DESCRIPTION:
    !
    ! create an object of type bgc_reaction_CENTURY_ECA_type.
    ! Right now it is purposely empty

  end function constructor


  !-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
    !
    ! DESCRIPTION:
    ! initialize boundary condition types
    ! !USES:
    use TracerBoundaryCondType      , only : tracerboundarycond_type
    use BeTRTracerType              , only : betrtracer_type

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECA_type), intent(in) :: this
    type(bounds_type)               , intent(in) :: bounds
    type(BeTRtracer_type )          ,  intent(in) :: betrtracer_vars
    type(tracerboundarycond_type)   ,  intent(in) :: tracerboundarycond_vars


    ! !LOCAL VARIABLES:
    integer :: c


    associate(                               &
         groupid  => betrtracer_vars%groupid &
         )

      tracerboundarycond_vars%topbc_type(1:betrtracer_vars%ngwmobile_tracer_groups) = bndcond_as_conc
      tracerboundarycond_vars%topbc_type(groupid(betrtracer_vars%id_trc_no3x)) = bndcond_as_flux

      tracerboundarycond_vars%topbc_type(betrtracer_vars%ngwmobile_tracer_groups+1:betrtracer_vars%ntracer_groups) = bndcond_as_flux

    end associate
  end subroutine init_boundary_condition_type

  !-------------------------------------------------------------------------------
  
  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars)
    !
    ! DESCRIPTION:
    ! initialize the betrbgc
    !
    ! !USES:
    use CNSharedParamsMod                , only : CNParamsReadShared
    use ncdio_pio                        , only : file_desc_t
    use BeTRTracerType                   , only : betrtracer_type
    use MathfuncMod                      , only : addone
    use clm_varctl                       , only : cnallocate_carbon_only_set

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECA_type) , intent(in)    :: this
    type(bounds_type)                    , intent(in)    :: bounds
    integer                              , intent(in)    :: lbj, ubj        ! lower and upper bounds, make sure they are > 0
    type(BeTRtracer_type )               , intent(inout) :: betrtracer_vars !

    ! !LOCAL VARIABLES:
    character(len=32), parameter                         :: subname ='Init_betrbgc'
    integer                                              :: jj
    integer                                              :: nelm, itemp_mem
    integer                                              :: itemp, itemp_vgrp, itemp_v, itemp_grp
    integer                                              :: c_loc, n_loc, trcid
    logical                                              :: carbon_only = .false.

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
    betrtracer_vars%nsolid_passive_tracer_groups =  4                             ! som1, som2, som3 and others (lit1, lit2, lit3, cwd)
    betrtracer_vars%nsolid_passive_tracers       = centurybgc_vars%nom_pools*nelm !

    betrtracer_vars%nmem_max                     = nelm*4                         ! total number of elemnts, and 4 sub members (lit1, lit2, lit3, cwd)

    call betrtracer_vars%Init()

    betrtracer_vars%is_mobile(:) = .true.

    jj = itemp
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
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT1C'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%lit1-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT1N'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
         trc_group_mem= addone(itemp_mem))

    !------------------------------------------------------------------------------------

    trcid = jj+(centurybgc_vars%lit2-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT2C'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%lit2-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT2N'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
         trc_group_mem= addone(itemp_mem))

    !------------------------------------------------------------------------------------

    trcid = jj+(centurybgc_vars%lit3-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT3C'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%lit3-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='LIT3N'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
         trc_group_mem= addone(itemp_mem))
    !------------------------------------------------------------------------------------

    trcid = jj+(centurybgc_vars%cwd-1 )*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='CWDC'              ,    &
         is_trc_mobile=.false., is_trc_advective = .false., trc_group_id = itemp_grp, &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%cwd-1 )*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='CWDN'              ,    &
         is_trc_mobile=.false., is_trc_advective = .false., trc_group_id = itemp_grp, &
         trc_group_mem= addone(itemp_mem))
    !==========================================================================================
    !new group
    itemp_mem = 0
    itemp_grp = addone(itemp_grp)
    trcid = jj+(centurybgc_vars%som1-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM1C'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp , &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%som1-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM1N'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
         trc_group_mem= addone(itemp_mem))
    !------------------------------------------------------------------------------------
    !new group
    itemp_mem=0
    itemp_grp = addone(itemp_grp)
    trcid = jj+(centurybgc_vars%som2-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM2C'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp , &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%som2-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM2N'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
         trc_group_mem= addone(itemp_mem))
    !------------------------------------------------------------------------------------
    !new group
    itemp_mem=0
    itemp_grp = addone(itemp_grp)
    trcid = jj+(centurybgc_vars%som3-1)*nelm+c_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM3C'            ,     &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
         trc_group_mem= addone(itemp_mem))

    trcid = jj+(centurybgc_vars%som3-1)*nelm+n_loc
    call betrtracer_vars%set_tracer(trc_id = trcid, trc_name='SOM3N'             ,    &
         is_trc_mobile=.true., is_trc_advective = .false., trc_group_id = itemp_grp,  &
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
    class(bgc_reaction_CENTURY_ECA_type) , intent(in)    :: this
    type(bounds_type)                    , intent(in)    :: bounds
    integer                              , intent(in)    :: num_soilc               ! number of columns in column filter
    integer                              , intent(in)    :: filter_soilc(:)         ! column filter
    type(betrtracer_type)                , intent(in)    :: betrtracer_vars
    real(r8)                             , intent(in)    :: dz_top(bounds%begc: )
    type(waterflux_type)                 , intent(in)    :: waterflux_vars
    type(tracerboundarycond_type)        , intent(inout) :: tracerboundarycond_vars !

    ! !LOCAL VARIABLES:
    character(len=255) :: subname = 'set_boundary_conditions'
    integer :: fc, c

    SHR_ASSERT_ALL((ubound(dz_top)                == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))

    associate(                                     &
         groupid  => betrtracer_vars%groupid          &
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
    ! !ARGUMENTS
    class(bgc_reaction_CENTURY_ECA_type) , intent(in) :: this
    type(bounds_type)                    , intent(in) :: bounds                        ! bounds
    integer                              , intent(in) :: num_soilc                     ! number of columns in column filter
    integer                              , intent(in) :: filter_soilc(:)               ! column filter
    integer                              , intent(in) :: num_soilp
    integer                              , intent(in) :: filter_soilp(:)               ! pft filter
    integer                              , intent(in) :: jtops(bounds%begc: )          ! top index of each column
    integer                              , intent(in) :: lbj, ubj                      ! lower and upper bounds, make sure they are > 0
    real(r8)                             , intent(in) :: dtime                         ! model time step
    type(Waterstate_Type)                , intent(in) :: waterstate_vars               ! water state variables
    type(temperature_type)               , intent(in) :: temperature_vars              ! energy state variable
    type(soilstate_type)                 , intent(in) :: soilstate_vars
    type(chemstate_type)                 , intent(in) :: chemstate_vars
    type(betrtracer_type)                , intent(in) :: betrtracer_vars               ! betr configuration information
    type(tracercoeff_type)               , intent(in) :: tracercoeff_vars
    type(carbonstate_type)               , intent(in) :: carbonstate_vars
    type(cnstate_type)                   , intent(inout) :: cnstate_vars
    type(carbonflux_type)                , intent(inout) :: carbonflux_vars
    type(nitrogenstate_type)             , intent(inout) :: nitrogenstate_vars
    type(nitrogenflux_type)              , intent(inout) :: nitrogenflux_vars
    type(tracerstate_type)               , intent(inout) :: tracerstate_vars
    type(tracerflux_type)                , intent(inout) :: tracerflux_vars
    type(plantsoilnutrientflux_type)     , intent(inout) :: plantsoilnutrientflux_vars !

    ! !LOCAL VARIABLES:
    character(len=32), parameter :: subname ='calc_bgc_reaction'
    integer                      :: fc, c, j, k
    real(r8)                     :: time
    real(r8)                     :: y0(centurybgc_vars%nstvars, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                     :: yf(centurybgc_vars%nstvars, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                     :: cn_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                     :: cp_ratios(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                     :: k_decay(centurybgc_vars%nreactions, bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                     :: pot_decay_rates(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj) ![mol C/m3/s] potential decay rates for different om pools without nutrient limitation
    real(r8)                     :: pot_co2_hr(bounds%begc:bounds%endc, lbj:ubj)                                 ![mol C/m3/s], potential co2 respiration rate
    real(r8)                     :: pot_nh3_immob(bounds%begc:bounds%endc,lbj:ubj)
    real(r8)                     :: anaerobic_frac(bounds%begc:bounds%endc,lbj:ubj)
    real(r8)                     :: n2_n2o_ratio_denit(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                     :: nh4_no3_ratio(bounds%begc:bounds%endc, lbj:ubj)
    real(r8)                     :: nuptake_prof(bounds%begc:bounds%endc,1:ubj)
    real(r8)                     :: pscal

    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))

    call Extra_inst%Init_Allocate(centurybgc_vars%nom_pools, centurybgc_vars%nreactions, centurybgc_vars%nprimvars)

    call set_reaction_order( centurybgc_vars%nreactions, centurybgc_vars, Extra_inst%is_zero_order)

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


    call bgcstate_ext_update_bfdecomp(bounds, 1, ubj, num_soilc, filter_soilc, &
         carbonflux_vars, nitrogenflux_vars, centurybgc_vars, betrtracer_vars, tracerflux_vars, y0, cn_ratios, cp_ratios)

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
    call calc_decompK_multiply_scalar(bounds, lbj, ubj, num_soilc, filter_soilc, jtops,                        &
         waterstate_vars%finundated_col(bounds%begc:bounds%endc), col%z(bounds%begc:bounds%endc, lbj:ubj),     &
         temperature_vars%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj),                                      &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_o2), &
         tracercoeff_vars%aqu2bulkcef_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_o2), &
         soilstate_vars, centurybgc_vars, carbonflux_vars)

    !calculate decay coefficients
    call calc_som_deacyK(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%nom_pools, &
         tracercoeff_vars, tracerstate_vars, betrtracer_vars, centurybgc_vars, carbonflux_vars,dtime, &
         k_decay(1:centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj))

    !calculate potential decay rates, without nutrient constraint
    call calc_sompool_decay(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars, &
         k_decay(1:centurybgc_vars%nom_pools,  bounds%begc:bounds%endc, lbj:ubj),              &
         y0(1:centurybgc_vars%nom_totelms, bounds%begc:bounds%endc, lbj:ubj),                  &
         pot_decay_rates)

    !calculate potential respiration rates by summarizing all om decomposition pathways
    call calc_potential_aerobic_hr(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, cn_ratios, cp_ratios, &
         centurybgc_vars, pot_decay_rates, soilstate_vars%cellsand_col(bounds%begc:bounds%endc,lbj:ubj),   &
         pot_co2_hr, pot_nh3_immob)

    !calculate fraction of anerobic environment
    call calc_anaerobic_frac(bounds, lbj, ubj, num_soilc, filter_soilc, jtops,                                 &
         temperature_vars%t_soisno_col(bounds%begc:bounds%endc,lbj:ubj),                                       &
         soilstate_vars, waterstate_vars%h2osoi_vol_col(bounds%begc:bounds%endc,lbj:ubj),                      &
         pot_co2_hr,                                                                                           &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_o2), &
         anaerobic_frac(bounds%begc:bounds%endc, lbj:ubj))

    !calculate normalized rate for nitrification and denitrification
    call calc_nitrif_denitrif_rate(bounds, lbj, ubj, num_soilc, filter_soilc, jtops,                             &
         col%dz(bounds%begc:bounds%endc, lbj:ubj),                                                               &
         temperature_vars%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj),                                        &
         chemstate_vars%soil_pH(bounds%begc:bounds%endc, lbj:ubj),                                               &
         pot_co2_hr,                                                                                             &
         anaerobic_frac,                                                                                         &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_nh3x), &
         tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_no3x), &
         soilstate_vars,                                                                                         &
         waterstate_vars,                                                                                        &
         carbonflux_vars,                                                                                        &
         n2_n2o_ratio_denit,                                                                                     &
         nh4_no3_ratio,                                                                                          &
         k_decay(centurybgc_vars%lid_nh4_nit_reac, bounds%begc:bounds%endc, lbj:ubj),                            &
         k_decay(centurybgc_vars%lid_no3_den_reac, bounds%begc:bounds%endc, lbj:ubj))

    !now there is no plant nitrogen uptake, I tend to create a new structure to indicate plant nutrient demand when it is hooked
    !back with CLM

    call calc_plant_nitrogen_uptake_prof(bounds, ubj, num_soilc, filter_soilc, col%dz(bounds%begc:bounds%endc, lbj:ubj), &
         plantsoilnutrientflux_vars%plant_minn_uptake_potential_col(bounds%begc:bounds%endc),                            &
         nuptake_prof(bounds%begc:bounds%endc,1:ubj),                                                                    &
         k_decay(centurybgc_vars%lid_plant_minn_up_reac, bounds%begc:bounds%endc ,1:ubj))

    !apply root distribution here
    call apply_plant_root_respiration_prof(bounds, ubj, num_soilc, filter_soilc,                                          &
         carbonflux_vars%rr_col(bounds%begc:bounds%endc), cnstate_vars%nfixation_prof_col(bounds%begc:bounds%endc,1:ubj), &
         k_decay(centurybgc_vars%lid_at_rt_reac, bounds%begc:bounds%endc, 1:ubj))

    call  apply_plant_root_nuptake_prof(bounds, ubj, num_soilc, filter_soilc    ,     &
         cnstate_vars%nfixation_prof_col(bounds%begc:bounds%endc,1:ubj)             , &
         plantsoilnutrientflux_vars)

    !do ode integration and update state variables for each layer

    do j = lbj, ubj
       do fc = 1, num_soilc
          c = filter_soilc(fc)
          if(j<jtops(c))cycle
          !assign parameters for stoichiometric matrix calculation
          call Extra_inst%AAssign(cn_r=cn_ratios(:,c,j),cp_r=cp_ratios(:,c,j), k_d=k_decay(:,c,j) ,    &
               n2_n2o_r_denit=n2_n2o_ratio_denit(c,j)                                                , &
               cell_sand=soilstate_vars%cellsand_col(c,j), cell_clay=soilstate_vars%cellclay_col(c,j), &
               plant_froots= plantsoilnutrientflux_vars%plant_frootsc_vr_col(c,j)                    , &
               betrtracer_vars=betrtracer_vars                                                       , &
               gas2bulkcef=tracercoeff_vars%gas2bulkcef_mobile_col(c,j,:)                            , &
               aere_cond=tracercoeff_vars%aere_cond_col(c,:), tracer_conc_atm=tracerstate_vars%tracer_conc_atm_col(c,:))

          !update state variables
          time = 0._r8

          yf(:,c,j)=y0(:,c,j) !this will allow to turn off the bgc reaction for debugging purpose
          call ode_ebbks1(one_box_century_bgc, y0(:,c,j), centurybgc_vars%nprimvars,centurybgc_vars%nstvars, &
               time, dtime, yf(:,c,j), pscal)

          if(pscal<5.e-1_r8)then
             write(iulog,*)'lat, lon=',grc%latdeg(col%gridcell(c)),grc%londeg(col%gridcell(c))
             write(iulog,*)'col, lev, pscal=',c, j, pscal
             write(iulog,*)'nstep =',get_nstep()
             call endrun()
          endif

       enddo
    enddo

    call bgcstate_ext_update_afdecomp(bounds, 1, ubj, num_soilc, filter_soilc, carbonflux_vars, nitrogenflux_vars, &
         centurybgc_vars, betrtracer_vars, tracerflux_vars, yf)

    !retrieve the flux variable
    call retrieve_flux_vars(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%nstvars, dtime, yf, y0, &
         centurybgc_vars, betrtracer_vars, tracerflux_vars, carbonflux_vars, nitrogenflux_vars, plantsoilnutrientflux_vars)

    !retrieve the state variable, state variable will be updated later
    call retrieve_state_vars(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%nstvars,  yf, &
         centurybgc_vars, betrtracer_vars, tracerstate_vars)

    call Extra_inst%DDeallocate()

  end subroutine calc_bgc_reaction

  !-------------------------------------------------------------------------------
  subroutine do_tracer_equilibration(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
       betrtracer_vars, tracercoeff_vars, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! equilibrate tracers that has solid and mobile phases
    ! using the theory of mass action.
    !
    ! !USES:
    !
    use tracerstatetype       , only : tracerstate_type
    use tracercoeffType       , only : tracercoeff_type
    use BeTRTracerType        , only : betrtracer_type

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECA_type), intent(in)    :: this
    type(bounds_type),                    intent(in)    :: bounds
    integer,                              intent(in)    :: lbj, ubj
    integer,                              intent(in)    :: jtops(bounds%begc: )        ! top label of each column
    integer,                              intent(in)    :: num_soilc
    integer,                              intent(in)    :: filter_soilc(:)
    type(betrtracer_type),                intent(in)    :: betrtracer_vars
    type(tracercoeff_type),               intent(in)    :: tracercoeff_vars
    type(tracerstate_type),               intent(inout) :: tracerstate_vars
    !
    ! !LOCAL VARIABLES:
    character(len=255) :: subname = 'do_tracer_equilibration'

    SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))

  end subroutine do_tracer_equilibration

  !-----------------------------------------------------------------------
  subroutine readParams(this, ncid, betrtracer_vars )
    !
    ! !DESCRIPTION:
    ! read model parameters
    ! !USES:
    use BeTRTracerType   , only : BeTRTracer_Type
    use ncdio_pio        , only : file_desc_t
    use BGCCenturyParMod , only : readCentDecompBgcParams, readCentNitrifDenitrifParams, readCentCNAllocParams

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECA_type) , intent(in)    :: this
    type(BeTRTracer_Type)                , intent(inout) :: betrtracer_vars
    type(file_desc_t)                    , intent(inout)  :: ncid  ! pio netCDF file id

    call readCentDecompBgcParams ( ncid, centurybgc_vars%nelms, betrtracer_vars )

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
    use BeTRTracerType           , only : BeTRTracer_Type
    use tracerstatetype          , only : tracerstate_type
    use WaterstateType           , only : waterstate_type
    use LandunitType             , only : lun
    use ColumnType               , only : col
    use PatchType                , only : pft
    use clm_varcon               , only : spval, ispval

    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_ECA_type) , intent(in)    :: this
    type(bounds_type)                    , intent(in)    :: bounds
    type(BeTRTracer_Type)                , intent(in)    :: betrtracer_vars
    type(waterstate_type)                , intent(in)    :: waterstate_vars
    type(tracerstate_type)               , intent(inout) :: tracerstate_vars

    !
    ! !LOCAL VARIABLES:
    integer :: p, c, l, k, j
    integer :: fc                                        ! filter index
    integer               :: begc, endc
    integer               :: begg, endg
    !-----------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begg = bounds%begg; endg= bounds%endg
    !-----------------------------------------------------------------------

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

            tracerstate_vars%tracer_conc_mobile_col    (c,:, :                                   )  = 0._r8
            tracerstate_vars%tracer_conc_surfwater_col (c,:                                      )  = 0._r8
            tracerstate_vars%tracer_conc_aquifer_col   (c,:                                      )  = 0._r8
            tracerstate_vars%tracer_conc_grndwater_col (c,:                                      )  = 0._r8
            tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_n2   )) = 32.8_r8
            tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_o2   )) = 8.78_r8
            tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_ar   )) = 0.3924_r8
            tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_co2x )) = 0.0168_r8
            tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_ch4  )) = 6.939e-5_r8
            tracerstate_vars%tracer_conc_atm_col       (c,volatileid(betrtracer_vars%id_trc_n2o  )) = 1.195e-5_r8

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
    class(bgc_reaction_CENTURY_ECA_type) , intent(in)    :: this
    type(bounds_type)                    , intent(in)    :: bounds          ! bounds
    integer                              , intent(in)    :: num_soilc       ! number of columns in column filter
    integer                              , intent(in)    :: filter_soilc(:) ! column filter
    type(betrtracer_type)                , intent(in)    :: betrtracer_vars ! betr configuration information
    type(tracerstate_type)               , intent(in)    :: tracerstate_vars
    type(tracerflux_type)                , intent(in)    :: tracerflux_vars
    type(carbonstate_type)               , intent(inout) :: carbonstate_vars
    type(nitrogenflux_type)              , intent(inout) :: nitrogenflux_vars
    type(nitrogenstate_type)             , intent(inout) :: nitrogenstate_vars


    call assign_nitrogen_hydroloss(bounds, num_soilc, filter_soilc, &
         tracerflux_vars, nitrogenflux_vars, betrtracer_vars)

    call assign_OM_CNpools(bounds, num_soilc, filter_soilc, &
         carbonstate_vars, nitrogenstate_vars, tracerstate_vars, betrtracer_vars, centurybgc_vars)

  end subroutine betr_alm_flux_statevar_feedback

  !---------------------------------------------------------------
  subroutine init_betr_alm_bgc_coupler(this, bounds, carbonstate_vars, &
       nitrogenstate_vars, betrtracer_vars, tracerstate_vars)
    !
    ! !DESCRIPTION:
    ! do state variable exchange between betr and alm
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
    class(bgc_reaction_CENTURY_ECA_type) , intent(in)    :: this
    type(bounds_type)                    , intent(in)    :: bounds
    type(tracerstate_type)               , intent(inout) :: tracerstate_vars
    type(betrtracer_type)                , intent(in)    :: betrtracer_vars ! betr configuration information
    type(carbonstate_type)               , intent(in)    :: carbonstate_vars
    type(nitrogenstate_type)             , intent(in)    :: nitrogenstate_vars

    ! !LOCAL VARIABLES:
    integer, parameter :: i_soil1 = 5
    integer, parameter :: i_soil2 = 6
    integer, parameter :: i_soil3 = 7
    character(len=255) :: subname = 'init_betr_alm_bgc_coupler'
    integer            :: c, j, k, l

    associate(                                                                         &
         id_trc_no3x        => betrtracer_vars%id_trc_no3x                           , &
         id_trc_nh3x        => betrtracer_vars%id_trc_nh3x                           , &
         decomp_cpools_vr   => carbonstate_vars%decomp_cpools_vr_col                 , &
         decomp_npools_vr   => nitrogenstate_vars%decomp_npools_vr_col               , &
         smin_no3_vr_col    => nitrogenstate_vars%smin_no3_vr_col                    , &
         smin_nh4_vr_col    => nitrogenstate_vars%smin_nh4_vr_col                    , &
         tracer_conc_mobile => tracerstate_vars%tracer_conc_mobile_col               , &
         tracer_conc_solid_passive => tracerstate_vars%tracer_conc_solid_passive_col , &
         c_loc              => centurybgc_vars%c_loc                                 , &
         n_loc              => centurybgc_vars%n_loc                                 , &
         lit1               => centurybgc_vars%lit1                                  , &
         lit2               => centurybgc_vars%lit2                                  , &
         lit3               => centurybgc_vars%lit3                                  , &
         som1               => centurybgc_vars%som1                                  , &
         som2               => centurybgc_vars%som2                                  , &
         som3               => centurybgc_vars%som3                                  , &
         cwd                => centurybgc_vars%cwd                                   , &
         nelms              => centurybgc_vars%nelms                                   &
         )

      !initialize tracer based on carbon/nitrogen pools
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
    ! the equations to be solved are in the form
    !
    ! dx/dt=I+A*R, where I is the input, A is the stoichiometric matrix, and R is the reaction vector
    !
    ! the input only contains litter input and mineral nutrient, som is assumed to be of fixed stoichiometry
    ! !USES:
    use SOMStateVarUpdateMod   , only : calc_dtrend_som_bgc
    use BGCCenturySubMod       , only : calc_cascade_matrix
    use MathfuncMod            , only : pd_decomp
    implicit none
    !
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
    real(r8) :: cascade_matrixp(nprimvars, Extra_inst%nr)
    real(r8) :: cascade_matrixd(nprimvars, Extra_inst%nr)
    logical  :: nitrogen_limit_flag(centurybgc_vars%nom_pools)
    real(r8) :: reaction_rates(Extra_inst%nr)
    real(r8) :: pscal(1:nprimvars)
    real(r8) :: rscal(1:Extra_inst%nr)
    real(r8) :: p_dt(1:nprimvars)
    real(r8) :: d_dt(1:nprimvars)
    integer  :: it
    logical  :: lneg

    !calculate cascade matrix, which contains the stoichiometry for all reactions
    call calc_cascade_matrix(nstvars, Extra_inst%nr, Extra_inst%cn_ratios, Extra_inst%cp_ratios,   &
         Extra_inst%n2_n2o_ratio_denit, Extra_inst%pct_sand, centurybgc_vars, nitrogen_limit_flag, &
         cascade_matrix)


    !obtain reaction rates
    do lk = 1, Extra_inst%nr
       if(Extra_inst%is_zero_order(lk))then

          if ( spinup_state .eq. 1 ) then
             !spinup stage
             if(lk == centurybgc_vars%lid_o2_aere_reac)then
                jj = centurybgc_vars%lid_o2
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                !avoid excessive arenchyma o2 transport into the atmosphere
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)
             else
                reaction_rates(lk) = Extra_inst%k_decay(lk)            !this effective defines the plant nitrogen demand
             endif
          else
             ! normal run stage
             if(lk == centurybgc_vars%lid_o2_aere_reac)then
                jj = centurybgc_vars%lid_o2
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                !avoid excessive arenchyma o2 transport into the atmosphere
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)
             elseif(lk == centurybgc_vars%lid_ch4_aere_reac)then
                jj = centurybgc_vars%lid_ch4
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                !avoid excessive arenchyma ch4 transport into the atmosphere
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)
             elseif(lk == centurybgc_vars%lid_ar_aere_reac)then
                jj = centurybgc_vars%lid_ar
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                !avoid excessive arenchyma ar transport into the atmosphere
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)
             elseif(lk == centurybgc_vars%lid_n2_aere_reac)then
                jj = centurybgc_vars%lid_n2
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                !avoid excessive arenchyma n2 transport into the atmosphere
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)
             elseif(lk == centurybgc_vars%lid_co2_aere_reac)then
                jj = centurybgc_vars%lid_co2
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                !avoid excessive arenchyma co2 transport into the atmosphere
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)
             elseif(lk == centurybgc_vars%lid_n2o_aere_reac)then
                jj = centurybgc_vars%lid_n2o
                reaction_rates(lk) = Extra_inst%scal_f(jj) *(Extra_inst%conv_f(jj)*ystate(jj) - Extra_inst%conc_f(jj))
                !avoid excessive arenchyma n2o transport into the atmosphere
                reaction_rates(lk) = min(reaction_rates(lk),ystate(jj)/dtime)
             else
                reaction_rates(lk) = Extra_inst%k_decay(lk)            !this effective defines the plant nitrogen demand
             endif
          endif

       else
          reaction_rates(lk)=ystate(centurybgc_vars%primvarid(lk))*Extra_inst%k_decay(lk)
       endif
    enddo

    call apply_ECA_nutrient_regulation(nprimvars, Extra_inst%nr, Extra_inst%pct_clay, nitrogen_limit_flag,  ystate(1:nprimvars), &
         Extra_inst%plant_frts, reaction_rates(1:Extra_inst%nr), cascade_matrix(1:nprimvars, 1:Extra_inst%nr))

    call pd_decomp(nprimvars, Extra_inst%nr, cascade_matrix(1:nprimvars, 1:Extra_inst%nr), &
         cascade_matrixp(1:nprimvars, 1:Extra_inst%nr),  cascade_matrixd(1:nprimvars, 1:Extra_inst%nr))
    it=0
    do
       call calc_dtrend_som_bgc(nprimvars, Extra_inst%nr, cascade_matrixp(1:nprimvars, 1:Extra_inst%nr), reaction_rates(1:Extra_inst%nr), p_dt)

       call calc_dtrend_som_bgc(nprimvars, Extra_inst%nr, cascade_matrixd(1:nprimvars, 1:Extra_inst%nr), reaction_rates(1:Extra_inst%nr), d_dt)


       !update the state variables
       call calc_pscal(nprimvars, dtime, ystate(1:nprimvars), p_dt(1:nprimvars), d_dt(1:nprimvars), pscal(1:nprimvars), lneg)

       if(lneg)then

          call calc_rscal(nprimvars, Extra_inst%nr, pscal, cascade_matrixd(1:nprimvars, 1:Extra_inst%nr), rscal)

          call reduce_reaction_rates(Extra_inst%nr, rscal(1:Extra_inst%nr), reaction_rates(1:Extra_inst%nr))
       else
          exit
       endif
       it = it + 1
       if(it>100)then
          write(iulog,*)'it',it
          call endrun('too many iterations')
       endif
    enddo

    call calc_dtrend_som_bgc(nstvars, Extra_inst%nr, cascade_matrix(1:nstvars, 1:Extra_inst%nr), &
         reaction_rates(1:Extra_inst%nr), dydt)

  end subroutine one_box_century_bgc
  !-------------------------------------------------------------------------------

  subroutine calc_pscal(nprimvars, dtime, ystate, p_dt,  d_dt, pscal, lneg)
    !
    ! !DESCRIPTION:
    ! calcualte limiting factor from each primary state variable
    !
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in)  :: nprimvars
    real(r8), intent(in)  :: dtime
    real(r8), intent(in)  :: ystate(1:nprimvars)
    real(r8), intent(in)  :: p_dt(1:nprimvars)
    real(r8), intent(in)  :: d_dt(1:nprimvars)
    real(r8), intent(out) :: pscal(1:nprimvars)
    logical,  intent(out) :: lneg

    ! !LOCAL VARIABLES:
    real(r8) :: yt
    real(r8) :: bb=0.999_r8
    integer  :: j
    lneg =.false.

    do j = 1, nprimvars
       yt = ystate(j) + (p_dt(j)+d_dt(j))*dtime
       if(yt<0._r8)then
          pscal(j) = -(p_dt(j)*dtime+ystate(j))/(dtime*d_dt(j))*bb
          lneg=.true.
          if(pscal(j)<0._r8)then
             call endrun('ngeative p in calc_pscal')
          endif
       else
          pscal(j) = 1._r8
       endif
    enddo
  end subroutine calc_pscal



  !-------------------------------------------------------------------------------
  subroutine calc_rscal(nprimvars, nr, pscal, cascade_matrixd, rscal)
    !
    ! !DESCRIPTION:
    ! calcualte limiting factor for each reaction
    ! !USES:
    use MathfuncMod , only : minp
    implicit none
    ! !ARGUMENTS:
    integer , intent(in) :: nprimvars
    integer , intent(in) :: nr
    real(r8), intent(in) :: pscal(1:nprimvars)
    real(r8), intent(in) :: cascade_matrixd(1:nprimvars, 1:nr)
    real(r8), intent(out):: rscal(1:nr)

    ! !LOCAL VARIABLES:
    integer :: j

    do j = 1, nr
       rscal(j) = minp(pscal,cascade_matrixd(1:nprimvars, j))
    enddo

  end subroutine calc_rscal

  !-------------------------------------------------------------------------------
  subroutine  reduce_reaction_rates(nr, rscal, reaction_rates)
    !
    ! !DESCRIPTION:
    ! reduce reaction rates using input scalar
    !
    implicit none
    ! !ARGUMENTS:
    integer , intent(in)    :: nr
    real(r8), intent(in)    :: rscal(1:nr)
    real(r8), intent(inout) :: reaction_rates(1:nr)
    ! !LOCAL VARIABLES:
    integer :: j

    do j = 1, nr
       reaction_rates(j) = reaction_rates(j)*rscal(j)
    enddo
  end subroutine  reduce_reaction_rates

  !-------------------------------------------------------------------------------
  subroutine apply_ECA_nutrient_regulation(nprimvars, nr, pct_clay, nitrogen_limit_flag, &
       ystate, plant_frts, reaction_rates, cascade_matrix)
    !
    ! !DESCRIPTION:
    ! do ECA competition
    !
    ! !USES:
    use KineticsMod, only : kd_infty, ecacomplex_cell_norm
    use MathfuncMod, only : safe_div
    implicit none
    ! !ARGUMENTS:
    integer , intent(in) :: nprimvars
    integer , intent(in) :: nr
    real(r8), intent(in) :: pct_clay
    logical , intent(in) :: nitrogen_limit_flag(centurybgc_vars%nom_pools)
    real(r8), intent(in) :: ystate(1:nprimvars)
    real(r8), intent(in) :: plant_frts
    real(r8), intent(inout) :: reaction_rates(1:nr)
    real(r8), intent(inout) :: cascade_matrix(1:nprimvars, 1:nr)

    ! !LOCAL VARIABLES:
    real(r8) :: k_mat(2,1:centurybgc_vars%ncompets)
    real(r8) :: vcompet(1:centurybgc_vars%ncompets)
    real(r8) :: siej_cell_norm(2, 1:centurybgc_vars%ncompets)
    real(r8) :: eca_nh4, eca_no3
    integer  :: j

    ! the following parameters are arbitrary
    real(r8), parameter :: kd_nh4_nit   = 1._r8
    real(r8), parameter :: kd_nh4_plant = 1._r8
    real(r8), parameter :: kd_nh4_clay  = 1._r8
    real(r8), parameter :: kd_no3_denit = 1._r8
    real(r8), parameter :: kd_no3_plant = 1._r8
    real(r8), parameter :: kd_nh4_decomp= 1._r8
    real(r8), parameter :: kd_no3_decomp= 1._r8


    !assume microbial biomass are 1% of the respective som pool
    !assume the conversion factor between clay (%) and NH4 adorsption capacity is gamma
    !also assume no competition between NH4 adsoprtion and other chemical adsorption
    !nh4 adsoprtion follows linear isotherm

    associate(                                                             &
         lid_nitri_compet   => centurybgc_vars%lid_nitri_compet          , &
         lid_denit_compet   => centurybgc_vars%lid_denit_compet          , &
         lid_plant_compet   => centurybgc_vars%lid_plant_compet          , &
         lid_clay_compet    => centurybgc_vars%lid_clay_compet           , &
         lid_lit1_compet    => centurybgc_vars%lid_lit1_compet           , &
         lid_lit2_compet    => centurybgc_vars%lid_lit2_compet           , &
         lid_lit3_compet    => centurybgc_vars%lid_lit3_compet           , &
         lid_cwd_compet     => centurybgc_vars%lid_cwd_compet            , &
         lid_som1_compet    => centurybgc_vars%lid_som1_compet           , &
         lid_som2_compet    => centurybgc_vars%lid_som2_compet           , &
         lid_som3_compet    => centurybgc_vars%lid_som3_compet           , &
         lit1               => centurybgc_vars%lit1                      , &
         lit2               => centurybgc_vars%lit2                      , &
         lit3               => centurybgc_vars%lit3                      , &
         cwd                => centurybgc_vars%cwd                       , &
         som1               => centurybgc_vars%som1                      , &
         som2               => centurybgc_vars%som2                      , &
         som3               => centurybgc_vars%som3                      , &
         lid_nh4            => centurybgc_vars%lid_nh4                   , &
         lid_no3            => centurybgc_vars%lid_no3                   , &
         lid_nh4_nit_reac   => centurybgc_vars%lid_nh4_nit_reac          , &
         lid_no3_den_reac   => centurybgc_vars%lid_no3_den_reac          , &
         lid_plant_minn_up_reac=> centurybgc_vars%lid_plant_minn_up_reac , &
         nelms              => centurybgc_vars%nelms                     , &
         c_loc              => centurybgc_vars%c_loc                       &
         )

      !form the K matrix
      
      !nh4
      k_mat(1,:) = kd_infty
      !no3
      k_mat(2,:) = kd_infty
      vcompet=0._r8
      do j = 1,  centurybgc_vars%nom_pools
         if(nitrogen_limit_flag(j))then
            k_mat(1,j) = kd_nh4_decomp
            k_mat(2,j) = kd_no3_decomp
         endif
      enddo

      k_mat(1,lid_nitri_compet) = kd_nh4_nit
      k_mat(1,lid_plant_compet) = kd_nh4_plant
      k_mat(1,lid_clay_compet)  = kd_nh4_clay
      !
      k_mat(2,lid_denit_compet) = kd_no3_denit
      k_mat(2,lid_plant_compet) = kd_no3_plant

      !form the competitor vector
      vcompet(lid_lit1_compet) = ystate((lit1-1)*nelms+c_loc) * 0.01_r8
      vcompet(lid_lit2_compet) = ystate((lit2-1)*nelms+c_loc) * 0.01_r8
      vcompet(lid_lit3_compet) = ystate((lit3-1)*nelms+c_loc) * 0.01_r8
      vcompet(lid_cwd_compet)  = ystate((cwd-1)*nelms+c_loc) * 0.01_r8

      ! by default som decomposition releases mineral nutrient, but I include them as
      ! subject to potential change
      vcompet(lid_som1_compet) = ystate((som1-1)*nelms+c_loc) * 0.01_r8
      vcompet(lid_som2_compet) = ystate((som2-1)*nelms+c_loc) * 0.01_r8
      vcompet(lid_som3_compet) = ystate((som3-1)*nelms+c_loc) * 0.01_r8

      vcompet(lid_nitri_compet)= ystate(lid_nh4) * 1.e-3_r8    !this number is arbitrary
      vcompet(lid_plant_compet)= plant_frts
      vcompet(lid_denit_compet)= ystate(lid_no3) * 1.e-3_r8    !this number is arbitrary
      vcompet(lid_clay_compet) = 1._r8
      !form the resource vector
      call ecacomplex_cell_norm(k_mat,(/ystate(lid_nh4),ystate(lid_no3)/),vcompet,siej_cell_norm)

      !now modify the reaction rates
      do j = 1,  centurybgc_vars%nom_pools
         if(nitrogen_limit_flag(j))then
            eca_nh4 = siej_cell_norm(1,j)/kd_nh4_decomp
            eca_no3 = siej_cell_norm(2,j)/kd_no3_decomp

            reaction_rates(j) = reaction_rates(j) * (eca_nh4 + eca_no3)
            cascade_matrix(lid_no3, j) = cascade_matrix(lid_nh4,j) * safe_div(eca_no3,eca_nh4+eca_no3)
            cascade_matrix(lid_nh4, j) = cascade_matrix(lid_nh4, j) - cascade_matrix(lid_no3,j)

         endif
      enddo

      !adjust for nitrification
      reaction_rates(lid_nh4_nit_reac) = reaction_rates(lid_nh4_nit_reac) * siej_cell_norm(1,lid_nitri_compet)
      !adjust for denitrification
      reaction_rates(lid_no3_den_reac) = reaction_rates(lid_no3_den_reac) * siej_cell_norm(2,lid_denit_compet)

      !adjust for plant mineral nitrogen uptake
      eca_nh4 = siej_cell_norm(1,lid_plant_compet)/kd_nh4_decomp
      eca_no3 = siej_cell_norm(2,lid_plant_compet)/kd_no3_decomp


      reaction_rates(lid_plant_compet) = reaction_rates(lid_plant_compet) * (eca_nh4+eca_no3) * vcompet(lid_plant_compet)
      cascade_matrix(lid_no3, lid_plant_minn_up_reac) = cascade_matrix(lid_nh4, lid_plant_minn_up_reac) * &
           safe_div(eca_no3, eca_nh4+eca_no3)
      cascade_matrix(lid_nh4,lid_plant_minn_up_reac) = cascade_matrix(lid_nh4,lid_plant_minn_up_reac) - &
           cascade_matrix(lid_no3, lid_plant_minn_up_reac)

    end associate
  end subroutine apply_ECA_nutrient_regulation

end module BGCReactionsCenturyECAType
