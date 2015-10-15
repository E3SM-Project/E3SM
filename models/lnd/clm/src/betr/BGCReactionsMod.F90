module BGCReactionsMod
  !
  ! !DESCRIPTION:
  ! template for doing bgc reaction in betr
  !
  ! !USES:
  use LandunitType , only : lun
  use ColumnType   , only : col
  implicit none
  save
  private
  public ::  bgc_reaction_type

  type, abstract :: bgc_reaction_type
     private
   contains
     !initialize betr bgc
     procedure(Init_betrbgc_interface)                    , deferred :: Init_betrbgc

     !doing bgc reaction
     procedure(calc_bgc_reaction_interface)               , deferred :: calc_bgc_reaction

     !set boundary condition for related tracer transport
     procedure(set_boundary_conditions_interface)         , deferred :: set_boundary_conditions

     procedure(init_boundary_condition_type_interface)    , deferred :: init_boundary_condition_type

     !do equilibrium tracer chemistry
     procedure(do_tracer_equilibration_interface )        , deferred :: do_tracer_equilibration

     !do cold initialization of different tracers
     procedure(initCold_interface)                        , deferred :: initCold

     !read in implementation specific parameters
     procedure(readParams_interface)                      , deferred :: readParams

     !send back state flux variables to other parts of alm
     procedure(betr_alm_flux_statevar_feedback_interface) , deferred :: betr_alm_flux_statevar_feedback

     !initialize betr state variable from other bgc components in alm
     procedure(init_betr_alm_bgc_coupler_interface)       , deferred :: init_betr_alm_bgc_coupler

  end type bgc_reaction_type

  abstract interface
     !----------------------------------------------------------------------
     subroutine Init_betrbgc_interface(this, bounds, lbj, ubj, betrtracer_vars)
       !
       ! !DESCRIPTION:
       ! template for init_betrbgc
       !
       ! !USES:
       use BeTRTracerType        , only : BeTRtracer_type
       use decompMod             , only : bounds_type
       !
       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type) , intent(in)    :: this
       type(bounds_type)        , intent(in)    :: bounds
       integer                  , intent(in)    :: lbj, ubj
       type(BeTRtracer_type )   , intent(inout) :: betrtracer_vars

     end subroutine Init_betrbgc_interface
     !----------------------------------------------------------------------
     subroutine calc_bgc_reaction_interface(this, bounds, lbj, ubj, num_soilc, filter_soilc, &
          num_soilp,filter_soilp, jtops, dtime,                                              &
          betrtracer_vars, tracercoeff_vars, waterstate_vars, temperature_vars,              &
          soilstate_vars, chemstate_vars, cnstate_vars, carbonstate_vars, carbonflux_vars,   &
          nitrogenstate_vars, nitrogenflux_vars,                                             &
          tracerstate_vars, tracerflux_vars, plantsoilnutrientflux_vars)
       !
       ! !DESCRIPTION:
       ! template for calc_bgc_reaction
       !
       ! !USES:
       use tracerfluxType           , only : tracerflux_type
       use tracerstatetype          , only : tracerstate_type
       use tracercoeffType          , only : tracercoeff_type
       use WaterStateType           , only : Waterstate_Type
       use TemperatureType          , only : temperature_type
       use ChemStateType            , only : chemstate_type
       use SoilStatetype            , only : soilstate_type
       use decompMod                , only : bounds_type
       use PlantSoilnutrientFluxType, only : plantsoilnutrientflux_type
       use BeTRTracerType           , only : BeTRTracer_Type
       use shr_kind_mod             , only : r8 => shr_kind_r8
       use CanopyStateType          , only : canopystate_type
       use CNStateType              , only : cnstate_type
       use CNCarbonStateType        , only : carbonstate_type
       use CNCarbonFluxType         , only : carbonflux_type
       use CNNitrogenFluxType       , only : nitrogenflux_type
       use CNNitrogenStateType      , only : nitrogenstate_type

       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)         , intent(in)    :: this
       type(bounds_type)                , intent(in)    :: bounds                      ! bounds
       integer                          , intent(in)    :: num_soilc                   ! number of columns in column filter
       integer                          , intent(in)    :: filter_soilc(:)             ! column filter
       integer                          , intent(in)    :: num_soilp
       integer                          , intent(in)    :: filter_soilp(:)
       integer                          , intent(in)    :: jtops( : )                  ! top index of each column
       integer                          , intent(in)    :: lbj, ubj                    ! lower and upper bounds, make sure they are > 0
       real(r8)                         , intent(in)    :: dtime                       ! model time step
       type(Waterstate_Type)            , intent(in)    :: waterstate_vars             ! water state variables
       type(temperature_type)           , intent(in)    :: temperature_vars            ! energy state variable
       type(chemstate_type)             , intent(in)    :: chemstate_vars
       type(betrtracer_type)            , intent(in)    :: betrtracer_vars             ! betr configuration information
       type(soilstate_type)             , intent(in)    :: soilstate_vars
       type(cnstate_type)               , intent(inout) :: cnstate_vars
       type(carbonstate_type)           , intent(in)    :: carbonstate_vars
       type(carbonflux_type)            , intent(inout) :: carbonflux_vars
       type(nitrogenstate_type)         , intent(inout) :: nitrogenstate_vars
       type(nitrogenflux_type)          , intent(inout) :: nitrogenflux_vars
       type(tracercoeff_type)           , intent(in)    :: tracercoeff_vars
       type(tracerstate_type)           , intent(inout) :: tracerstate_vars
       type(tracerflux_type)            , intent(inout) :: tracerflux_vars
       type(plantsoilnutrientflux_type) , intent(inout) ::  plantsoilnutrientflux_vars !total nitrogen yield to plant

     end subroutine calc_bgc_reaction_interface
     !----------------------------------------------------------------------

     subroutine set_boundary_conditions_interface(this, bounds, num_soilc, filter_soilc, dz_top, &
          betrtracer_vars, waterflux_vars, tracerboundarycond_vars)

       ! !DESCRIPTION:
       ! template for set_boundary_conditions
       !
       ! !USES:
       use TracerBoundaryCondType , only : tracerboundarycond_type
       use decompMod              , only : bounds_type
       use BeTRTracerType         , only : BeTRTracer_Type
       use WaterfluxType          , only : waterflux_type
       use shr_kind_mod           , only : r8 => shr_kind_r8

       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)      , intent(in) :: this                       !
       type(bounds_type)             , intent(in) :: bounds                     !
       integer                       , intent(in) :: num_soilc                  ! number of columns in column filter
       integer                       , intent(in) :: filter_soilc(:)            ! column filter
       type(betrtracer_type)         , intent(in) :: betrtracer_vars            !
       real(r8)                      , intent(in) :: dz_top( : )                !
       type(waterflux_type)          , intent(in) :: waterflux_vars             !
       type(tracerboundarycond_type) , intent(inout) :: tracerboundarycond_vars !

     end subroutine set_boundary_conditions_interface

     !----------------------------------------------------------------------

     subroutine init_boundary_condition_type_interface(this, bounds, &
          betrtracer_vars, tracerboundarycond_vars )
       !
       ! !DESCRIPTION:
       ! template for init_boundary_condition
       !
       ! !USES:
       use BeTRTracerType        , only : betrtracer_type
       use TracerBoundaryCondType, only : tracerboundarycond_type
       use decompMod             , only : bounds_type

       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)         , intent(in) :: this
       type(BeTRtracer_type )           , intent(in) :: betrtracer_vars
       type(bounds_type)                , intent(in) :: bounds
       type(tracerboundarycond_type)    , intent(in) :: tracerboundarycond_vars

     end subroutine init_boundary_condition_type_interface


     !-------------------------------------------------------------------------------
     subroutine do_tracer_equilibration_interface(this, bounds, lbj, ubj, jtops, num_soilc, filter_soilc, &
          betrtracer_vars, tracercoeff_vars, tracerstate_vars)
       !
       ! !DESCRIPTION:
       ! template for do_tracer_equilibration
       ! !USES:
       !
       use tracerstatetype       , only : tracerstate_type
       use tracercoeffType       , only : tracercoeff_type
       use BeTRTracerType        , only : BeTRTracer_Type
       use decompMod             , only : bounds_type

       ! !ARGUMENTS:
       import :: bgc_reaction_type

       class(bgc_reaction_type)   , intent(in)    :: this
       type(bounds_type)          , intent(in)    :: bounds
       integer                    , intent(in)    :: lbj, ubj
       integer                    , intent(in)    :: jtops( : )        ! top label of each column
       integer                    , intent(in)    :: num_soilc
       integer                    , intent(in)    :: filter_soilc(:)
       type(betrtracer_type)      , intent(in)    :: betrtracer_vars
       type(tracercoeff_type)     , intent(in)    :: tracercoeff_vars
       type(tracerstate_type)     , intent(inout) :: tracerstate_vars


     end subroutine do_tracer_equilibration_interface

     !-------------------------------------------------------------------------------
     subroutine InitCold_interface(this, bounds, betrtracer_vars, waterstate_vars, tracerstate_vars)
       !
       ! !DESCRIPTION:
       ! template for InitCold
       ! !USES:
       !
       use BeTRTracerType           , only : BeTRTracer_Type
       use tracerstatetype          , only : tracerstate_type
       use WaterstateType           , only : waterstate_type
       use LandunitType             , only : lun
       use ColumnType               , only : col
       use PatchType                , only : pft
       use decompMod                , only : bounds_type

       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)          , intent(in)    :: this
       type(bounds_type)                 , intent(in)    :: bounds
       type(BeTRTracer_Type)             , intent(in)    :: betrtracer_vars
       type(waterstate_type)             , intent(in)    :: waterstate_vars
       type(tracerstate_type)            , intent(inout) :: tracerstate_vars


     end subroutine InitCold_interface

     !-------------------------------------------------------------------------------
     subroutine readParams_interface(this, ncid, betrtracer_vars)
       !
       ! !DESCRIPTION:
       ! template for readParams
       ! !USES:
       use ncdio_pio                , only : file_desc_t
       use BeTRTracerType           , only : BeTRTracer_Type

       ! !ARGUMENTS:
       import :: bgc_reaction_type

       class(bgc_reaction_type)          , intent(in)    :: this
       type(file_desc_t)                 , intent(inout) :: ncid  ! pio netCDF file id
       type(BeTRTracer_Type)             , intent(inout) :: betrtracer_vars

     end subroutine readParams_interface

     !-------------------------------------------------------------------------------
     subroutine betr_alm_flux_statevar_feedback_interface(this, bounds, num_soilc, filter_soilc, &
          carbonstate_vars, nitrogenstate_vars, nitrogenflux_vars, &
          tracerstate_vars, tracerflux_vars,  betrtracer_vars)

       ! !DESCRIPTION:
       ! template for betr_alm_flux_statevar_feedback
       ! !USES:
       use decompMod                , only : bounds_type
       use shr_kind_mod             , only : r8 => shr_kind_r8
       use tracerfluxType           , only : tracerflux_type
       use tracerstatetype          , only : tracerstate_type
       use BeTRTracerType           , only : BeTRTracer_Type
       use CNCarbonStateType        , only : carbonstate_type
       use CNNitrogenStateType      , only : nitrogenstate_type
       use CNNitrogenFluxType       , only : nitrogenflux_type

       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)   , intent(in) :: this                  !
       type(bounds_type)          , intent(in) :: bounds                ! bounds
       integer                    , intent(in) :: num_soilc             ! number of columns in column filter
       integer                    , intent(in) :: filter_soilc(:)       ! column filter
       type(betrtracer_type)      , intent(in) :: betrtracer_vars       ! betr configuration information
       type(tracerstate_type)     , intent(in) :: tracerstate_vars      !
       type(tracerflux_type)      , intent(in) :: tracerflux_vars       !
       type(carbonstate_type)     , intent(inout) :: carbonstate_vars   !
       type(nitrogenflux_type)    , intent(inout) :: nitrogenflux_vars  !
       type(nitrogenstate_type)   , intent(inout) :: nitrogenstate_vars !

     end subroutine betr_alm_flux_statevar_feedback_interface

     !-------------------------------------------------------------------------------


     subroutine init_betr_alm_bgc_coupler_interface(this, bounds, carbonstate_vars, &
          nitrogenstate_vars, betrtracer_vars, tracerstate_vars)
       !
       ! !DESCRIPTION:
       ! template for init_betr_alm_bgc_coupler

       ! !USES:
       use decompMod                , only : bounds_type
       use clm_varcon               , only : natomw, catomw
       use clm_varpar               , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
       use CNCarbonStateType        , only : carbonstate_type
       use CNNitrogenStateType      , only : nitrogenstate_type
       use tracerstatetype          , only : tracerstate_type
       use BetrTracerType           , only : betrtracer_type
       use clm_varpar               , only : nlevtrc_soil
       use landunit_varcon          , only : istsoil, istcrop
       !
       ! !ARGUMENTS:
       import :: bgc_reaction_type
       class(bgc_reaction_type)           , intent(in)    :: this
       type(bounds_type)                  , intent(in)    :: bounds
       type(tracerstate_type)             , intent(inout) :: tracerstate_vars
       type(betrtracer_type)              , intent(in)    :: betrtracer_vars    ! betr configuration information
       type(carbonstate_type)             , intent(in)    :: carbonstate_vars
       type(nitrogenstate_type)           , intent(in)    :: nitrogenstate_vars


     end subroutine init_betr_alm_bgc_coupler_interface

  end interface
end module BGCReactionsMod
