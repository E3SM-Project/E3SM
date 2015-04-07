module NutrientCompetitionMethodMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Abstract base class for functions to calculate nutrient competition
  !
  ! Created by Jinyun Tang, following Bill Sack's implementation of polymorphism
  ! !USES:
  implicit none
  private
  !
  ! !PUBLIC TYPES:
  public :: nutrient_competition_method_type

  type, abstract :: nutrient_competition_method_type
     private
   contains

     ! read in nutrient competition kinetic parameters
     procedure(readParams_interface), deferred :: readParams

      ! compute plant nutrient demand
     procedure(calc_plant_nutrient_demand_interface), deferred :: calc_plant_nutrient_demand

     ! compute the nutrient yield for different components
     procedure(calc_plant_nutrient_competition_interface), deferred :: calc_plant_nutrient_competition

  end type nutrient_competition_method_type

  abstract interface

     ! Note: The following code is adapted based on what Bill Scaks has done for soil water retention curve
     ! polymorphism. Therefore, I also keep some suggestions he gave there.
     !
     ! - Make the interfaces contain all possible inputs that are needed by any
     !   implementation; each implementation will then ignore the inputs it doesn't need.
     !
     ! - For inputs that are needed only by particular implementations - and particularly
     !   for inputs that are constant in time 
     !   pass these into the constructor, and save pointers to these inputs as components
     !   of the child type that needs them. Then they aren't needed as inputs to the
     !   individual routines, allowing the interfaces for these routines to remain more
     !   consistent between different implementations.
     !
     !---------------------------------------------------------------------------
     subroutine readParams_interface(this, ncid)
       ! !DESCRIPTION:
       ! read in kinetic parameters that are needed for doing nutrient competition
       !
       ! !USES:
       use ncdio_pio, only : file_desc_t
       import :: nutrient_competition_method_type
       !
       ! !ARGUMENTS:
       class(nutrient_competition_method_type) , intent(in)    :: this
       type(file_desc_t)                       , intent(inout) :: ncid   ! pio netCDF file id

     end subroutine readParams_interface

     !---------------------------------------------------------------------------     
     subroutine calc_plant_nutrient_demand_interface (this, bounds, num_soilp, filter_soilp, &
          photosyns_inst, crop_inst, canopystate_inst,                             &
          cnveg_state_inst, cnveg_carbonstate_inst, cnveg_carbonflux_inst,         &
          c13_cnveg_carbonflux_inst, c14_cnveg_carbonflux_inst,                    &
          cnveg_nitrogenstate_inst, cnveg_nitrogenflux_inst, aroot, arepr)
       !
       ! DESCRIPTION
       ! calculate nutrient yield after considering competition between different components
       !
       ! USES
       use shr_kind_mod           , only : r8 => shr_kind_r8
       use decompMod              , only : bounds_type       
       use PhotosynthesisMod      , only : photosyns_type
       use CropType               , only : crop_type
       use CanopyStateType        , only : canopystate_type
       use CNVegStateType         , only : cnveg_state_type
       use CNVegCarbonStateType   , only : cnveg_carbonstate_type
       use CNVegNitrogenStateType , only : cnveg_nitrogenstate_type
       use CNVegCarbonFluxType    , only : cnveg_carbonflux_type
       use CNVegNitrogenFluxType  , only : cnveg_nitrogenflux_type
       import :: nutrient_competition_method_type
       !
       ! !ARGUMENTS:
       class(nutrient_competition_method_type) , intent(in)    :: this
       type(bounds_type)               , intent(in)    :: bounds
       integer                         , intent(in)    :: num_soilp        ! number of soil patches in filter
       integer                         , intent(in)    :: filter_soilp(:)  ! filter for soil patches
       type(photosyns_type)            , intent(in)    :: photosyns_inst
       type(crop_type)                 , intent(in)    :: crop_inst
       type(canopystate_type)          , intent(in)    :: canopystate_inst
       type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
       type(cnveg_carbonstate_type)    , intent(inout) :: cnveg_carbonstate_inst
       type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
       type(cnveg_carbonflux_type)     , intent(inout) :: c13_cnveg_carbonflux_inst
       type(cnveg_carbonflux_type)     , intent(inout) :: c14_cnveg_carbonflux_inst
       type(cnveg_nitrogenstate_type)  , intent(in)    :: cnveg_nitrogenstate_inst
       type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
       real(r8)                        , intent(out)   :: aroot(bounds%begp:)
       real(r8)                        , intent(out)   :: arepr(bounds%begp:)

     end subroutine calc_plant_nutrient_demand_interface

     !-----------------------------------------------------------------------
     subroutine calc_plant_nutrient_competition_interface (this, bounds, num_soilp, filter_soilp, &
          cnveg_state_inst, cnveg_carbonflux_inst, c13_cnveg_carbonflux_inst,                     &
          c14_cnveg_carbonflux_inst, cnveg_nitrogenflux_inst,                                     &
          aroot, arepr, fpg_col)                                              
       !
       ! !USES:
       use shr_kind_mod          , only : r8 => shr_kind_r8
       use decompMod             , only : bounds_type       
       use CNVegStateType        , only : cnveg_state_type
       use CNVegCarbonFluxType   , only : cnveg_carbonflux_type
       use CNVegNitrogenFluxType , only : cnveg_nitrogenflux_type
       import :: nutrient_competition_method_type
       !
       ! !ARGUMENTS:
       class(nutrient_competition_method_type), intent(in) :: this
       type(bounds_type)               , intent(in)    :: bounds
       integer                         , intent(in)    :: num_soilp        ! number of soil patches in filter
       integer                         , intent(in)    :: filter_soilp(:)  ! filter for soil patches
       type(cnveg_state_type)          , intent(inout) :: cnveg_state_inst
       type(cnveg_carbonflux_type)     , intent(inout) :: cnveg_carbonflux_inst
       type(cnveg_carbonflux_type)     , intent(inout) :: c13_cnveg_carbonflux_inst
       type(cnveg_carbonflux_type)     , intent(inout) :: c14_cnveg_carbonflux_inst
       type(cnveg_nitrogenflux_type)   , intent(inout) :: cnveg_nitrogenflux_inst
       real(r8)                        , intent(in)    :: aroot(bounds%begp:)
       real(r8)                        , intent(in)    :: arepr(bounds%begp:)
       real(r8)                        , intent(in)    :: fpg_col(bounds%begc:)

     end subroutine calc_plant_nutrient_competition_interface

  end interface

end module NutrientCompetitionMethodMod
