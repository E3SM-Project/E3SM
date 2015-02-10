module BGCReactionsMod
!
! module doing bgc reaction
! created by Jinyun Tang
! This is dirty version just to make sure the model is running
! Eventually, I want to introduce polymorphism to make it
! consistent with other developments in soil hydrology and the clm4.5/clm5 bgc

implicit none
  save
  private
  public ::  bgc_reaction_type
  
  type, abstract :: bgc_reaction_type
      private
    contains
      !initialize betr bgc
      procedure(Init_betrbgc_interface)                 , deferred :: Init_betrbgc
      !doing bgc reaction
      procedure(calc_bgc_reaction_interface)            , deferred :: calc_bgc_reaction
      
      !set boundary condition for related tracer transport
      procedure(set_boundary_conditions_interface)      , deferred :: set_boundary_conditions
        
      procedure(init_boundary_condition_type_interface) , deferred :: init_boundary_condition_type
      
      !do equilibrium tracer chemistry
      procedure(do_tracer_equilibration_interface )     , deferred :: do_tracer_equilibration
      
      !do cold initialization of different tracers
      procedure(initCold_interface)                     , deferred :: initCold
  end type bgc_reaction_type  
  
  abstract interface
!----------------------------------------------------------------------      
    subroutine Init_betrbgc_interface(this, bounds, lbj, ubj, betrtracer_vars)
    !
    ! USES
    use BeTRTracerType        , only : BeTRtracer_type
    use decompMod             , only : bounds_type      
    !
    import :: bgc_reaction_type
    class(bgc_reaction_type),    intent(in) :: this    
    type(bounds_type),           intent(in) :: bounds
    integer,                     intent(in) :: lbj, ubj      
    type(BeTRtracer_type ),   intent(inout) :: betrtracer_vars
    
    
    end subroutine Init_betrbgc_interface
!----------------------------------------------------------------------    
    subroutine calc_bgc_reaction_interface(this, bounds, lbj, ubj, numf, filter, jtops, dtime, &
       col, betrtracer_vars, tracercoeff_vars, waterstate_vars, temperature_vars, &
       soilstate_vars, chemstate_vars, tracerstate_vars, tracerflux_vars, plantsoilnutrientflux_vars)
  !
  ! do bgc reaction
  ! eventually this will be an abstract subroutine, but now I use the select case approach for a quick and dirty implementation.
  !USES
  !
   use tracerfluxType           , only : tracerflux_type
   use tracerstatetype          , only : tracerstate_type
   use tracercoeffType          , only : tracercoeff_type  
   use WaterStateType           , only : Waterstate_Type  
   use TemperatureType          , only : temperature_type
   use ChemStateType            , only : chemstate_type
   use ColumnType               , only : column_type
   use SoilStatetype            , only : soilstate_type
   use decompMod                , only : bounds_type
   use PlantSoilnutrientFluxType, only : plantsoilnutrientflux_type
   use BeTRTracerType           , only : BeTRTracer_Type       
   use shr_kind_mod             , only : r8 => shr_kind_r8

   import :: bgc_reaction_type
   class(bgc_reaction_type),    intent(in) :: this
   type(bounds_type),           intent(in) :: bounds                             ! bounds   
   integer,                     intent(in) :: numf                               ! number of columns in column filter
   integer,                     intent(in) :: filter(:)                          ! column filter
   integer,                     intent(in) :: jtops( : )               ! top index of each column
   integer,                     intent(in) :: lbj, ubj                           ! lower and upper bounds, make sure they are > 0  
   real(r8),                    intent(in) :: dtime                              ! model time step
   type(column_type),           intent(in) :: col                                ! column type
   type(Waterstate_Type),       intent(in) :: waterstate_vars                    ! water state variables
   type(temperature_type),      intent(in) :: temperature_vars                   ! energy state variable
   type(chemstate_type),        intent(in) :: chemstate_vars
   type(betrtracer_type),       intent(in) :: betrtracer_vars                    ! betr configuration information
   type(soilstate_type)  ,      intent(in) :: soilstate_vars  
   type(tracercoeff_type),      intent(in) :: tracercoeff_vars
   type(tracerstate_type),   intent(inout) :: tracerstate_vars
   type(tracerflux_type),    intent(inout) :: tracerflux_vars
   type(plantsoilnutrientflux_type), intent(inout) ::  plantsoilnutrientflux_vars                    !total nitrogen yield to plant     
    
   end subroutine calc_bgc_reaction_interface
!----------------------------------------------------------------------    
   
   subroutine set_boundary_conditions_interface(this, bounds, numf, filter, dz_top, &
      betrtracer_vars, waterflux_vars, tracerboundarycond_vars)

   use TracerBoundaryCondType, only : tracerboundarycond_type   
   use decompMod                , only : bounds_type
   use BeTRTracerType           , only : BeTRTracer_Type    
   use WaterfluxType            , only : waterflux_type    
   use shr_kind_mod             , only : r8 => shr_kind_r8
   
   
   import :: bgc_reaction_type
   class(bgc_reaction_type)   , intent(in) :: this
   type(bounds_type)          , intent(in) :: bounds
   integer                    , intent(in) :: numf                 ! number of columns in column filter
   integer                    , intent(in) :: filter(:)            ! column filter   
   type(betrtracer_type)      , intent(in) :: betrtracer_vars
   real(r8)                   , intent(in) :: dz_top( : )
   type(waterflux_type)       , intent(in) :: waterflux_vars     
   type(tracerboundarycond_type), intent(inout) :: tracerboundarycond_vars

   
   
   end subroutine set_boundary_conditions_interface
   
!----------------------------------------------------------------------    

   subroutine init_boundary_condition_type_interface(this, bounds, betrtracer_vars, tracerboundarycond_vars )
   !
   ! DESCRIPTIONS
   ! initialize boundary condition types
   use BeTRTracerType        , only : betrtracer_type
   use TracerBoundaryCondType, only : tracerboundarycond_type
   use decompMod             , only : bounds_type   
   import :: bgc_reaction_type
   
   class(bgc_reaction_type),    intent(in) :: this    
   type(BeTRtracer_type ),            intent(in) :: betrtracer_vars
   type(bounds_type),                 intent(in) :: bounds
   type(tracerboundarycond_type),     intent(in) :: tracerboundarycond_vars
  
   end subroutine init_boundary_condition_type_interface
   
   
!-------------------------------------------------------------------------------  
  subroutine do_tracer_equilibration_interface(this, bounds, lbj, ubj, jtops, numf, filter, betrtracer_vars, &
    tracercoeff_vars, tracerstate_vars)
  !
  ! DESCRIPTIONS
  ! requilibrate tracers that has solid and mobile phases
  ! using the theory of mass action. When the redox-ladder is on, this
  ! subroutine will update the change of pH due to tracer transport, or 
  ! USES
  !
  use tracerstatetype       , only : tracerstate_type
  use tracercoeffType       , only : tracercoeff_type 
  use BeTRTracerType        , only : BeTRTracer_Type    
  use decompMod             , only : bounds_type   

  import :: bgc_reaction_type
  
  class(bgc_reaction_type)   , intent(in)    :: this  
  type(bounds_type)          , intent(in)    :: bounds
  integer                    , intent(in)    :: lbj, ubj
  integer                    , intent(in)    :: jtops( : )        ! top label of each column
  integer                    , intent(in)    :: numf
  integer                    , intent(in)    :: filter(:)
  type(betrtracer_type)      , intent(in)    :: betrtracer_vars
  type(tracercoeff_type)     , intent(in)    :: tracercoeff_vars
  type(tracerstate_type)     , intent(inout) :: tracerstate_vars
  
  
  end subroutine do_tracer_equilibration_interface 
  
!-------------------------------------------------------------------------------
  subroutine InitCold_interface(this, bounds, betrtracer_vars, waterstate_vars, tracerstate_vars)
  !
  ! !USES:
  !
  use BeTRTracerType           , only : BeTRTracer_Type
  use tracerstatetype          , only : tracerstate_type
  use WaterstateType           , only : waterstate_type  
  use LandunitType             , only : lun                
  use ColumnType               , only : col  
  use PatchType                , only : pft
  use decompMod                , only : bounds_type
    
  import :: bgc_reaction_type      

  ! !ARGUMENTS:
  class(bgc_reaction_type)          , intent(in)    :: this  
  type(bounds_type)                 , intent(in)    :: bounds
  type(BeTRTracer_Type)             , intent(in)    :: betrtracer_vars
  type(waterstate_type)             , intent(in)    :: waterstate_vars  
  type(tracerstate_type)            , intent(inout) :: tracerstate_vars
  
  
  end subroutine InitCold_interface  
  
  end interface
end module BGCReactionsMod
