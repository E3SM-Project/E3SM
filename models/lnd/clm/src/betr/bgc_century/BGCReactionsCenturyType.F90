module BGCReactionsCenturyType

#include "shr_assert.h"

!
! !DESCRIPTION
! This is an example on how to use polymorphism to create your own bgc modules that will be run with BeTR
!
! HISTORY:
! Created by Jinyun Tang, Oct 2nd, 2014
! Questions for thought/things to do: Should I consider redox flucutation here (specifically redox lag)? It also
! seems necessary to create a separate subroutine to account for that, rather than using the ch4 code.

! !USES

  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use decompMod             , only : bounds_type
  use BGCReactionsMod       , only : bgc_reaction_type
  use clm_varcon            , only : spval
  use tracer_varcon         , only : bndcond_as_conc, bndcond_as_flux
  use BGCCenturySubMod     
implicit none

  save
  private
  !
  ! !PUBLIC TYPES:
  public :: bgc_reaction_CENTURY_type

  type(centurybgc_type), private :: centurybgc_vars


  type, extends(bgc_reaction_type) :: &
  bgc_reaction_CENTURY_type
  private
  contains
   procedure :: Init_betrbgc                  ! initialize betr bgc
   procedure :: set_boundary_conditions       ! set top/bottom boundary conditions for various tracers
   procedure :: calc_bgc_reaction             ! doing bgc calculation
   procedure :: init_boundary_condition_type  ! initialize type of top boundary conditions
   procedure :: do_tracer_equilibration       ! do equilibrium tracer chemistry
   procedure :: initCold
   procedure :: readParams
  end type bgc_reaction_CENTURY_type

  
   type, private :: Extra_type
     real(r8), pointer :: cn_ratios(:)           !cn ratio of om pool
     real(r8), pointer :: cp_ratios(:)           !cp ratio of om pool
     real(r8), pointer :: k_decay(:)             !decay parameter for all reactions
     real(r8)          :: n2_n2o_ratio_denit     !ratio of n2 to n2o during denitrification
     real(r8)          :: nh4_no3_ratio          !ratio of available nh4 to no3
     real(r8)          :: cellsand               !sand content
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

   subroutine Init_Allocate(this, nompools, nreacts)
   
   class(Extra_type) :: this
   
   integer, intent(in) :: nompools
   integer, intent(in) :: nreacts
   
   allocate(this%cn_ratios(nompools))
   allocate(this%cp_ratios(nompools))
   allocate(this%k_decay(nreacts))
   this%nr = nreacts
   end subroutine Init_Allocate

!-------------------------------------------------------------------------------   

   subroutine DDeallocate(this)
   
   class(Extra_type) :: this
   
   
   deallocate(this%cn_ratios)
   deallocate(this%cp_ratios)
   deallocate(this%k_decay)
   end subroutine DDeallocate
!-------------------------------------------------------------------------------   

   subroutine AAssign(this, cn_r,cp_r, k_d,  n2_n2o_r_denit, nh4_no3_r, cell_sand)
   
   class(Extra_type) :: this
   real(r8), dimension(:), intent(in) :: cn_r
   real(r8), dimension(:), intent(in) :: cp_r
   real(r8), dimension(:), intent(in) :: k_d
   real(r8)              , intent(in) :: n2_n2o_r_denit
   real(r8)              , intent(in) :: nh4_no3_r
   real(r8)              , intent(in) :: cell_sand
   
   integer :: n1, n2, n3
   
   
   n1 = size(cn_r)
   n2 = size(cp_r)
   n3 = size(k_d)
   SHR_ASSERT_ALL((n1              == n2),        errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((n3              == this%nr),   errMsg(__FILE__,__LINE__))
   this%cn_ratios(1:n1) = cn_r
   this%cp_ratios(1:n2) = cp_r
   
   this%n2_n2o_ratio_denit = n2_n2o_r_denit
   this%nh4_no3_ratio      = nh4_no3_r
   this%cellsand           = cell_sand
   this%k_decay            = k_d  
   end subroutine AAssign
   
!-------------------------------------------------------------------------------
  type(bgc_reaction_CENTURY_type) function constructor()
  !
  ! ! DESCRIPTION
  ! create an object of type bgc_reaction_CENTURY_type.
  ! Right now it is purposely left empty

  end function constructor


!-------------------------------------------------------------------------------
  subroutine init_boundary_condition_type(this, bounds, betrtracer_vars, tracerboundarycond_vars )
  !
  ! DESCRIPTIONS
  ! initialize boundary condition types

  use TracerBoundaryCondType      , only : tracerboundarycond_type
  use BeTRTracerType              , only : betrtracer_type
  
  
  class(bgc_reaction_CENTURY_type), intent(in) :: this
  type(bounds_type)               , intent(in) :: bounds  
  type(BeTRtracer_type )          ,  intent(in) :: betrtracer_vars
  type(tracerboundarycond_type)   ,  intent(in) :: tracerboundarycond_vars



  integer :: c


  tracerboundarycond_vars%topbc_type(1:betrtracer_vars%ngwmobile_tracers) = bndcond_as_conc
    
  tracerboundarycond_vars%topbc_type(betrtracer_vars%ngwmobile_tracers+1:betrtracer_vars%ntracers) = bndcond_as_flux
  end subroutine init_boundary_condition_type

!-------------------------------------------------------------------------------

  subroutine Init_betrbgc(this, bounds, lbj, ubj, betrtracer_vars)
  !
  ! DESCRIPTION
  ! initialize the betrbgc
  use CNSharedParamsMod     , only : CNParamsReadShared
  use ncdio_pio             , only : file_desc_t
  use BeTRTracerType        , only : betrtracer_type
  
  
  class(bgc_reaction_CENTURY_type), intent(in) :: this
  type(bounds_type)               , intent(in) :: bounds
  integer                         , intent(in) :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0  
  type(BeTRtracer_type )          , intent(inout) :: betrtracer_vars
  
  !local variables
  character(len=*)                 , parameter :: subname ='Init_betrbgc'

  integer :: jj
  !type(file_desc_t) :: ncid
  
  !ncid%fh=10

  call centurybgc_vars%Init(bounds, lbj, ubj)

  betrtracer_vars%ngwmobile_tracers=7                                 ! n2, o2, ar, co2, ch4, nh3x and no3x
  betrtracer_vars%nsolid_passive_tracers=centurybgc_vars%nom_pools    !
  betrtracer_vars%nvolatile_tracers=6                                 ! n2, o2, ar, co2, ch4 and nh3x 
  betrtracer_vars%ntracers=betrtracer_vars%ngwmobile_tracers+betrtracer_vars%nsolid_passive_tracers
  

  call betrtracer_vars%Init()
  
  betrtracer_vars%id_trc_n2  = 1
  betrtracer_vars%id_trc_o2  = 2
  betrtracer_vars%id_trc_ar  = 3
  betrtracer_vars%id_trc_co2x= 4
  betrtracer_vars%id_trc_ch4 = 5
  betrtracer_vars%id_trc_nh3x = 6
  betrtracer_vars%id_trc_no3x = 7
  
  jj = 7
  betrtracer_vars%tracernames(betrtracer_vars%id_trc_n2)='N2'
  betrtracer_vars%tracernames(betrtracer_vars%id_trc_o2)='O2'
  betrtracer_vars%tracernames(betrtracer_vars%id_trc_ar)='AR'
  betrtracer_vars%tracernames(betrtracer_vars%id_trc_co2x)='CO2x'
  betrtracer_vars%tracernames(betrtracer_vars%id_trc_ch4)='CH4'
  betrtracer_vars%tracernames(betrtracer_vars%id_trc_nh3x)='NH3X'
  betrtracer_vars%tracernames(betrtracer_vars%id_trc_no3x)='NO3X'
 
  betrtracer_vars%tracernames(jj+1)  = 'LIT1'
  betrtracer_vars%tracernames(jj+2)  = 'LIT2'
  betrtracer_vars%tracernames(jj+3)  = 'LIT3'
  betrtracer_vars%tracernames(jj+4)  = 'SOM1'
  betrtracer_vars%tracernames(jj+5) = 'SOM2'
  betrtracer_vars%tracernames(jj+6) = 'SOM3'
  betrtracer_vars%tracernames(jj+7) = 'CWD'
  
  betrtracer_vars%is_volatile(betrtracer_vars%id_trc_n2)  =.true.
  betrtracer_vars%is_volatile(betrtracer_vars%id_trc_o2)  =.true.
  betrtracer_vars%is_volatile(betrtracer_vars%id_trc_ar)  =.true.
  betrtracer_vars%is_volatile(betrtracer_vars%id_trc_co2x)=.true.
  betrtracer_vars%is_volatile(betrtracer_vars%id_trc_ch4) =.true.
  
  betrtracer_vars%volatileid(betrtracer_vars%id_trc_n2)   = 1
  betrtracer_vars%volatileid(betrtracer_vars%id_trc_o2)   = 2
  betrtracer_vars%volatileid(betrtracer_vars%id_trc_ar)   = 3
  betrtracer_vars%volatileid(betrtracer_vars%id_trc_co2x) = 4
  betrtracer_vars%volatileid(betrtracer_vars%id_trc_ch4)  = 5
  
  betrtracer_vars%is_mobile(betrtracer_vars%id_trc_n2)   = .true.
  betrtracer_vars%is_mobile(betrtracer_vars%id_trc_o2)   = .true.
  betrtracer_vars%is_mobile(betrtracer_vars%id_trc_ar)   = .true.
  betrtracer_vars%is_mobile(betrtracer_vars%id_trc_co2x) = .true.
  betrtracer_vars%is_mobile(betrtracer_vars%id_trc_ch4)  = .true.
  betrtracer_vars%is_mobile(betrtracer_vars%id_trc_no3x)  = .true.
  
  betrtracer_vars%is_advective(betrtracer_vars%id_trc_n2)   = .true.
  betrtracer_vars%is_advective(betrtracer_vars%id_trc_o2)   = .true.
  betrtracer_vars%is_advective(betrtracer_vars%id_trc_ar)   = .true.
  betrtracer_vars%is_advective(betrtracer_vars%id_trc_co2x) = .true.
  betrtracer_vars%is_advective(betrtracer_vars%id_trc_ch4)  = .true.
  betrtracer_vars%is_advective(betrtracer_vars%id_trc_no3x)  = .true.
  
  
  
  !comment following lines out when it is hooked to CLM at the moment
  !call CNParamsReadShared(ncid)
  
  !call readCNNitrifDenitrifParams ( ncid )
  
  end subroutine Init_betrbgc

!-------------------------------------------------------------------------------
  subroutine set_boundary_conditions(this, bounds, num_soilc, filter_soilc, dz_top, betrtracer_vars, &
    waterflux_vars, tracerboundarycond_vars)
  !
  ! DESCRIPTION
  ! set up boundary conditions for tracer movement
  !

  use clm_varctl            , only : iulog
  use TracerBoundaryCondType, only : tracerboundarycond_type
  use abortutils            , only : endrun
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use BeTRTracerType        , only : betrtracer_type
  use WaterfluxType         , only : waterflux_type  

  
  class(bgc_reaction_CENTURY_type), intent(in) :: this
  type(bounds_type)               , intent(in) :: bounds
  integer                         , intent(in) :: num_soilc                 ! number of columns in column filter
  integer                         , intent(in) :: filter_soilc(:)            ! column filter
  type(betrtracer_type)           , intent(in) :: betrtracer_vars
  real(r8)                        , intent(in) :: dz_top(bounds%begc: )
  type(waterflux_type)            , intent(in) :: waterflux_vars
  type(tracerboundarycond_type)   , intent(inout) :: tracerboundarycond_vars


  !local variables
  character(len=255) :: subname = 'set_boundary_conditions'

  SHR_ASSERT_ALL((ubound(dz_top)                == (/bounds%endc/)),   errMsg(__FILE__,__LINE__))


  !eventually, the following code will be implemented using polymorphism
  tracerboundarycond_vars%tracer_gwdif_concflux_top_col(1,1:2,betrtracer_vars%ngwmobile_tracers+1:betrtracer_vars%ntracers)=0._r8       !zero incoming flux
  tracerboundarycond_vars%tracer_gwdif_concflux_top_col(1,1:2,betrtracer_vars%id_trc_n2)=32.8_r8                         !mol m-3, contant boundary condition
  tracerboundarycond_vars%tracer_gwdif_concflux_top_col(1,1:2,betrtracer_vars%id_trc_o2)=8.78_r8                         !mol m-3, contant boundary condition
  tracerboundarycond_vars%tracer_gwdif_concflux_top_col(1,1:2,betrtracer_vars%id_trc_ar)=0.3924_r8                       !mol m-3, contant boundary condition
  tracerboundarycond_vars%tracer_gwdif_concflux_top_col(1,1:2,betrtracer_vars%id_trc_co2x)=0.0168_r8                     !mol m-3, contant boundary condition  
  tracerboundarycond_vars%tracer_gwdif_concflux_top_col(1,1:2,betrtracer_vars%id_trc_ch4)=6.939e-5_r8                    !mol m-3, contant boundary condition

    
  tracerboundarycond_vars%bot_concflux_col(1,1,:) = 0._r8                                  !zero flux boundary condition
  tracerboundarycond_vars%condc_toplay_col(1,betrtracer_vars%id_trc_n2) = 2._r8*1.267e-5_r8/dz_top(1)     !m/s surface conductance
  tracerboundarycond_vars%condc_toplay_col(1,betrtracer_vars%id_trc_o2) = 2._r8*1.267e-5_r8/dz_top(1)     !m/s surface conductance
  tracerboundarycond_vars%condc_toplay_col(1,betrtracer_vars%id_trc_ar) = 2._r8*1.267e-5_r8/dz_top(1)     !m/s surface conductance
  tracerboundarycond_vars%condc_toplay_col(1,betrtracer_vars%id_trc_co2x) = 2._r8*1.267e-5_r8/dz_top(1)   !m/s surface conductance
  tracerboundarycond_vars%condc_toplay_col(1,betrtracer_vars%id_trc_ch4) = 2._r8*1.267e-5_r8/dz_top(1)    !m/s surface conductance  
  
  end subroutine set_boundary_conditions
!-------------------------------------------------------------------------------

  subroutine calc_bgc_reaction(this, bounds, lbj, ubj, num_soilc, filter_soilc, num_soilp, filter_soilp, jtops, dtime, col, &
    betrtracer_vars, tracercoeff_vars, waterstate_vars, temperature_vars, soilstate_vars, chemstate_vars, &
    cnstate_vars, canopystate_vars, tracerstate_vars, tracerflux_vars, plantsoilnutrientflux_vars)
  !
  ! do bgc reaction
  ! eventually this will be an abstract subroutine, but now I use the select case approach for a quick and dirty implementation.
  !USES
  !
  use tracerfluxType           , only : tracerflux_type
  use tracerstatetype          , only : tracerstate_type
  use tracercoeffType          , only : tracercoeff_type
  use BetrTracerType           , only : betrtracer_type
  use WaterStateType           , only : Waterstate_Type
  use TemperatureType          , only : temperature_type
  use ChemStateType            , only : chemstate_type
  use ColumnType               , only : column_type
  use SoilStatetype            , only : soilstate_type
  use ODEMod                   , only : ode_adapt_mbbks1
  use CanopyStateType          , only : canopystate_type
  use CNStateType              , only : cnstate_type  
  use PlantSoilnutrientFluxType, only : plantsoilnutrientflux_type
  use CNVerticalProfileMod     , only : decomp_vertprofiles
  !ARGUMENTS
  class(bgc_reaction_CENTURY_type)   , intent(in) :: this
  type(bounds_type)                  , intent(in) :: bounds                             ! bounds
  integer                            , intent(in) :: num_soilc                               ! number of columns in column filter
  integer                            , intent(in) :: filter_soilc(:)                          ! column filter
  integer                            , intent(in) :: num_soilp
  integer                            , intent(in) :: filter_soilp(:)                    ! pft filter
  integer                            , intent(in) :: jtops(bounds%begc: )               ! top index of each column
  integer                            , intent(in) :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0
  real(r8)                           , intent(in) :: dtime                              ! model time step
  type(column_type)                  , intent(in) :: col                                ! column type
  type(Waterstate_Type)              , intent(in) :: waterstate_vars                    ! water state variables
  type(temperature_type)             , intent(in) :: temperature_vars                   ! energy state variable
  type(soilstate_type)               , intent(in) :: soilstate_vars  
  type(chemstate_type)               , intent(in) :: chemstate_vars
  type(betrtracer_type)              , intent(in) :: betrtracer_vars                    ! betr configuration information
  type(tracercoeff_type)             , intent(in) :: tracercoeff_vars
  type(canopystate_type)             , intent(in)    :: canopystate_vars
  type(cnstate_type)                 , intent(inout) :: cnstate_vars
  type(tracerstate_type)             , intent(inout) :: tracerstate_vars
  type(tracerflux_type)              , intent(inout) :: tracerflux_vars
  type(plantsoilnutrientflux_type)   , intent(inout) :: plantsoilnutrientflux_vars
  
  character(len=*), parameter :: subname ='calc_bgc_reaction'

  integer :: fc, c, j

  real(r8) :: y0(centurybgc_vars%neqs, bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: yf(centurybgc_vars%neqs, bounds%begc:bounds%endc, lbj:ubj)  
  real(r8) :: cn_ratios(centurybgc_vars%nom_pools)
  real(r8) :: cp_ratios(centurybgc_vars%nom_pools)
  real(r8) :: time
  real(r8) :: k_decay(centurybgc_vars%nreactions, bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: pot_decay_rates(centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj)   ![mol C/m3/s] potential decay rates for different om pools without nutrient limitation
  real(r8) :: pot_co2_hr(bounds%begc:bounds%endc, lbj:ubj)                   ![mol C/m3/s], potential co2 respiration rate
  real(r8) :: anaerobic_frac(bounds%begc:bounds%endc,lbj:ubj)
  real(r8) :: n2_n2o_ratio_denit(bounds%begc:bounds%endc, lbj:ubj)
  real(r8) :: nh4_no3_ratio(bounds%begc:bounds%endc, lbj:ubj)
  
  SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))
  
  call Extra_inst%Init_Allocate(centurybgc_vars%nom_pools, centurybgc_vars%nreactions)
  
  !initialize local variables
  y0(:, :, :) = spval
  yf(:, :, :) = spval

  !calculate vertical profiles to destribute various variables
  call decomp_vertprofiles(bounds, &
           num_soilc, filter_soilc, num_soilp, filter_soilp, &
           soilstate_vars, canopystate_vars, cnstate_vars)
  
  !add new litter 
  call calc_om_input(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, dtime, col, betrtracer_vars, centurybgc_vars, cnstate_vars, tracerstate_vars)
      
  !initialize the state vector
  call init_state_vector(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%neqs, tracerstate_vars, betrtracer_vars, centurybgc_vars, y0)

  !calculate multiplicative scalars for decay parameters
  call calc_decompK_multiply_scalar(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, &
    waterstate_vars%finundated_col(bounds%begc:bounds%endc), col%z(bounds%begc:bounds%endc, lbj:ubj),&
    temperature_vars%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj), &
    tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_o2), &
    tracercoeff_vars%aqu2bulkcef_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_o2), &
     soilstate_vars, centurybgc_vars)
  
  !calculate decay coefficients
  call calc_som_deacyK(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%nom_pools, tracercoeff_vars, tracerstate_vars, &
    betrtracer_vars, centurybgc_vars, k_decay(1:centurybgc_vars%nom_pools, bounds%begc:bounds%endc, lbj:ubj))
  
  !calculate potential decay rates, without nutrient constraint
  call calc_sompool_decay(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%nom_pools, &
     k_decay(1:centurybgc_vars%nom_pools,  bounds%begc:bounds%endc, lbj:ubj), y0(1:centurybgc_vars%nom_pools, lbj:ubj, bounds%begc:bounds%endc),&
     pot_decay_rates)
      
  !calculate potential respiration rates by summarizing all om decomposition pathways
  call calc_potential_aerobic_hr(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%nom_pools, pot_decay_rates, &
    soilstate_vars%cellsand_col(bounds%begc:bounds%endc,lbj:ubj), pot_co2_hr)      
      
  !calculate fraction of anerobic environment
  call calc_anaerobic_frac(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, temperature_vars%t_soisno_col(bounds%begc:bounds%endc,lbj:ubj),&
     soilstate_vars, waterstate_vars%h2osoi_vol_col(bounds%begc:bounds%endc,lbj:ubj), pot_co2_hr, &
     tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_o2), &
     anaerobic_frac(bounds%begc:bounds%endc, lbj:ubj))
      
  !calculate normalized rate for nitrification and denitrification
  call calc_nitrif_denitrif_rate(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, col%dz(bounds%begc:bounds%endc, lbj:ubj), &
     temperature_vars%t_soisno_col(bounds%begc:bounds%endc, lbj:ubj), &
     chemstate_vars%soil_pH(bounds%begc:bounds%endc, lbj:ubj),  pot_co2_hr, anaerobic_frac, &
     tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_nh3x), &     
     tracerstate_vars%tracer_conc_mobile_col(bounds%begc:bounds%endc, lbj:ubj, betrtracer_vars%id_trc_no3x), &
     soilstate_vars, waterstate_vars, centurybgc_vars, n2_n2o_ratio_denit, nh4_no3_ratio, &
     k_decay(centurybgc_vars%lid_nh4, bounds%begc:bounds%endc, lbj:ubj), &
     k_decay(centurybgc_vars%lid_no3, bounds%begc:bounds%endc, lbj:ubj))
     
  !now there is no plant nitrogen uptake, I tend to create a new structure to indicate plant nutrient demand when it is hooked
  !back with CLM
  k_decay(centurybgc_vars%lid_plant_minn, : ,: )=0._r8    
  !do ode integration and update state variables for each layer
  do j = lbj, ubj
    do fc = 1, num_soilc
      c = filter_soilc(fc)
      if(j<jtops(c))cycle   
      !assign parameters for stoichiometric matrix calculation
      call Extra_inst%AAssign(cn_ratios,cp_ratios, k_decay(:,c,j), n2_n2o_ratio_denit(c,j), nh4_no3_ratio(c,j), soilstate_vars%cellsand_col(c,j))
      !update state variables
      time = 0._r8 
      call ode_adapt_mbbks1(one_box_century_bgc, y0(:,c,j), centurybgc_vars%neqs, time, dtime, yf(:,c,j))
    enddo
  enddo  

  
  !retrieve the flux variable
  call retrieve_flux(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%neqs, dtime, yf, y0, tracerflux_vars, centurybgc_vars, plantsoilnutrientflux_vars)
  
  !retrieve the state variable    
  call retrieve_state_vector(bounds, lbj, ubj, num_soilc, filter_soilc, jtops, centurybgc_vars%neqs,  yf, centurybgc_vars, betrtracer_vars, tracerstate_vars)      
  
  
  call Extra_inst%DDeallocate()
  end subroutine calc_bgc_reaction
  
    
  
!-------------------------------------------------------------------------------
  subroutine do_tracer_equilibration(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, betrtracer_vars, tracercoeff_vars, tracerstate_vars)
  !
  ! DESCRIPTIONS
  ! requilibrate tracers that has solid and mobile phases
  ! using the theory of mass action. When the redox-ladder is on, this
  ! subroutine will update the change of pH due to tracer transport, or
  ! USES
  !
  use tracerstatetype       , only : tracerstate_type
  use tracercoeffType       , only : tracercoeff_type
  use BeTRTracerType        , only : betrtracer_type
  
  class(bgc_reaction_CENTURY_type),    intent(in) :: this

  type(bounds_type),      intent(in) :: bounds
  integer,                intent(in) :: lbj, ubj
  integer,                intent(in) :: jtops(bounds%begc: )        ! top label of each column
  integer,                intent(in) :: num_soilc
  integer,                intent(in) :: filter_soilc(:)
  type(betrtracer_type),  intent(in) :: betrtracer_vars
  type(tracercoeff_type), intent(in) :: tracercoeff_vars
  type(tracerstate_type), intent(inout) :: tracerstate_vars
  character(len=255) :: subname = 'do_tracer_equilibration'


  SHR_ASSERT_ALL((ubound(jtops) == (/bounds%endc/)), errMsg(__FILE__,__LINE__))

  !depending on the simulation type, an implementation of aqueous chemistry will be
  !employed to separate out the adsorbed phase
  !It should be noted that this formulation excludes the use of linear isotherm, which
  !can be integrated through the retardation factor


  end subroutine do_tracer_equilibration
  
  !-----------------------------------------------------------------------
  subroutine readParams(this, ncid)

  use ncdio_pio               , only : file_desc_t
                                         
  class(bgc_reaction_CENTURY_type) , intent(in)    :: this  
  
  type(file_desc_t)  :: ncid  ! pio netCDF file id
  
  

  call readCentDecompBgcParams ( ncid )
  
  call readCentNitrifDenitrifParams ( ncid )
  
  end subroutine readParams
  
  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, betrtracer_vars,  waterstate_vars, tracerstate_vars)
    !
    ! !USES:
    !
    use BeTRTracerType           , only : BeTRTracer_Type
    use tracerstatetype          , only : tracerstate_type
    use WaterstateType           , only : waterstate_type        
    use LandunitType             , only : lun                
    use ColumnType               , only : col                
    use PatchType                , only : pft
    use clm_varcon               , only : spval, ispval
    use landunit_varcon          , only : istsoil, istcrop    
    
    ! !ARGUMENTS:
    class(bgc_reaction_CENTURY_type) , intent(in)    :: this    
    type(bounds_type)                 , intent(in)    :: bounds
    type(BeTRTracer_Type)             , intent(in)    :: betrtracer_vars
    type(waterstate_type)             , intent(in)    :: waterstate_vars    
    type(tracerstate_type)            , intent(inout) :: tracerstate_vars
    
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

    
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
         if(betrtracer_vars%ngwmobile_tracers>0)then
           tracerstate_vars%tracer_conc_mobile_col(c,:,:)        = spval
           tracerstate_vars%tracer_conc_surfwater_col(c,:)       = spval
           tracerstate_vars%tracer_conc_aquifer_col(c,:)         = spval
           tracerstate_vars%tracer_conc_grndwater_col(c,:)       = spval           
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

         tracerstate_vars%tracer_conc_mobile_col(c,:, :)          = 0._r8
         tracerstate_vars%tracer_conc_surfwater_col(c,:)          = 0._r8
         tracerstate_vars%tracer_conc_aquifer_col(c,:)            = 0._r8
         tracerstate_vars%tracer_conc_grndwater_col(c,:)          = 0._r8                      

      
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
    
  end subroutine InitCold


!-------------------------------------------------------------------------------


  subroutine one_box_century_bgc(ystate, dtime, time, neqs, dydt)
  !
  ! DESCRIPTIONS
  !do single box bgc
  !
  use SOMStateVarUpdateMod  , only : calc_dtrend_som_bgc
  use BGCCenturySubMod      , only : calc_cascade_matrix
  implicit none
  integer,  intent(in)  :: neqs
  real(r8), intent(in)  :: dtime
  real(r8), intent(in)  :: time
  real(r8), intent(in)  :: ystate(neqs)
  real(r8), intent(out) :: dydt(neqs)
  
  !local variables
  integer :: lk
  real(r8) :: cascade_matrix(neqs, Extra_inst%nr)
  real(r8) :: reaction_rates(Extra_inst%nr)
  
  
    !calculate cascade matrix, which contains the stoichiometry for all reactions
  call calc_cascade_matrix(neqs, Extra_inst%nr, Extra_inst%cn_ratios, Extra_inst%cp_ratios, &
      Extra_inst%n2_n2o_ratio_denit, Extra_inst%nh4_no3_ratio, Extra_inst%cellsand, centurybgc_vars, cascade_matrix)
    
    !do pool degradation
  do lk = 1, Extra_inst%nr
    reaction_rates(lk)=ystate(lk)*Extra_inst%k_decay(lk)
  enddo
    
  call calc_dtrend_som_bgc(neqs, Extra_inst%nr, cascade_matrix(1:neqs, 1:Extra_inst%nr), reaction_rates(1:Extra_inst%nr), dydt)

  end subroutine one_box_century_bgc
  
end module BGCReactionsCenturyType
