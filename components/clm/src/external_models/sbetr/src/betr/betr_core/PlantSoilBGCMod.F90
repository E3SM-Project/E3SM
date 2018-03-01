module PlantSoilBGCMod

!
! !DESCRIPTION:
! template for doing plant and soil bgc coupling
!
! !USES:
  use BeTR_decompMod         , only : betr_bounds_type
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
type, abstract :: plant_soilbgc_type
   private
   ! dummy var to remove compiler warnings
   logical, public :: dummy_compiler_warning
 contains

   !initialize Init_plant_soilbgc
   procedure(Init_plant_soilbgc_interface)                    , deferred :: Init_plant_soilbgc
   !send back plant nutrient yield
   procedure(lsm_betr_plant_soilbgc_recv_interface)           , deferred :: lsm_betr_plant_soilbgc_recv
   !summarize active+passive plant nutrient yield
   procedure(plant_soilbgc_summary_interface)                 , deferred :: plant_soilbgc_summary
   !do profile integration of fluxes
   procedure(integrate_vr_flux_interface)                     , deferred :: integrate_vr_flux
   !send in bgc inputs
   procedure(lsm_betr_plant_soilbgc_send_interface)          , deferred  :: lsm_betr_plant_soilbgc_send

end type plant_soilbgc_type


  abstract interface
  !----------------------------------------------------------------------
  subroutine Init_plant_soilbgc_interface(this, bounds, lbj, ubj, namelist_buffer)

  !
  ! !DESCRIPTION:
  ! template for init_betrbgc
  !
  ! !USES:
  use BeTR_decompMod         , only : betr_bounds_type

  ! !ARGUMENTS:
  import :: plant_soilbgc_type
  implicit none
  class(plant_soilbgc_type) , intent(inout) :: this
  type(betr_bounds_type)    , intent(in) :: bounds
  character(len=*)          , intent(in) :: namelist_buffer
  integer                   , intent(in) :: lbj, ubj

  end subroutine Init_plant_soilbgc_interface


  !----------------------------------------------------------------------
  subroutine plant_soilbgc_summary_interface(this,bounds, lbj, ubj, pft, numf, &
       filter, dtime, dz, betrtracer_vars, tracerflux_vars, biogeo_flux, betr_status)

  ! !USES:
  use BeTRTracerType , only : BeTRtracer_type
  use tracerfluxType , only : tracerflux_type
  use BeTR_decompMod , only : betr_bounds_type
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use BetrStatusType , only : betr_status_type
  use BeTR_PatchType , only : betr_patch_type
  use BeTR_biogeoFluxType  , only : betr_biogeo_flux_type
  ! !ARGUMENTS:
  import :: plant_soilbgc_type

  class(plant_soilbgc_type) , intent(inout) :: this
  type(betr_bounds_type)    , intent(in) :: bounds
  integer                   , intent(in) :: lbj, ubj
  integer                   , intent(in) :: numf
  type(betr_patch_type)     , intent(in) :: pft
  integer                   , intent(in) :: filter(:)
  real(r8)                  , intent(in) :: dtime
  real(r8)                  , intent(in) :: dz(bounds%begc: ,1: )
  type(BeTRtracer_type )    , intent(in) :: betrtracer_vars
  type(tracerflux_type)     , intent(in) :: tracerflux_vars
  type(betr_biogeo_flux_type)      , intent(inout) :: biogeo_flux
  type(betr_status_type)    , intent(out):: betr_status

  end subroutine plant_soilbgc_summary_interface


  !----------------------------------------------------------------------

  subroutine integrate_vr_flux_interface(this, bounds, numf, filter)
  !
  ! !DESCRIPTIONS
  ! integrate 3d fluxes into 2d fluxes
  use BeTR_decompMod         , only : betr_bounds_type
  ! !ARGUMENTS:
  import :: plant_soilbgc_type

  class(plant_soilbgc_type) , intent(inout) :: this
  type(betr_bounds_type)    , intent(in) :: bounds
  integer                   , intent(in) :: numf
  integer                   , intent(in) :: filter(:)


  end subroutine integrate_vr_flux_interface
  !----------------------------------------------------------------------

  subroutine lsm_betr_plant_soilbgc_recv_interface(this, bounds, numf, filter, &
     betr_pft, biogeo_fluxes)
  !DESCRIPTION
  !return plant nutrient yield
  !
  !USES
  use BeTR_decompMod       , only : betr_bounds_type
  use BeTR_biogeoStateType , only : betr_biogeo_state_type
  use BeTR_biogeoFluxType  , only : betr_biogeo_flux_type
  use BeTR_PatchType, only : betr_patch_type
  ! !ARGUMENTS:
  import :: plant_soilbgc_type
  implicit none
  class(plant_soilbgc_type)   , intent(inout)    :: this
  type(betr_bounds_type)      , intent(in)    :: bounds
  integer                     , intent(in)    :: numf
  integer                     , intent(in)    :: filter(:)
  type(betr_patch_type) , intent(in) :: betr_pft
  type(betr_biogeo_flux_type) , intent(inout) :: biogeo_fluxes


  end subroutine lsm_betr_plant_soilbgc_recv_interface
  !----------------------------------------------------------------------
  subroutine lsm_betr_plant_soilbgc_send_interface(this, bounds, numf, &
                 filter, betr_pft, biogeo_forc, biogeo_states, biogeo_fluxes)
  !DESCRIPTION
  ! initialize feedback variables for plant soil bgc interactions
  !
  !USES
  use BeTR_biogeoStateType , only : betr_biogeo_state_type
  use BeTR_biogeoFluxType  , only : betr_biogeo_flux_type
  use BeTR_decompMod       , only : betr_bounds_type
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use BeTR_PatchType, only : betr_patch_type
  ! !ARGUMENTS:
  import :: plant_soilbgc_type

  class(plant_soilbgc_type)    , intent(inout) :: this
  type(betr_bounds_type)       , intent(in) :: bounds
  integer                      , intent(in) :: numf
  integer                      , intent(in) :: filter(:)
  type(betr_patch_type)        , intent(in) :: betr_pft
  type(betr_biogeophys_input_type), intent(in):: biogeo_forc
  type(betr_biogeo_state_type) , intent(in) :: biogeo_states
  type(betr_biogeo_flux_type)  , intent(in) :: biogeo_fluxes

  end subroutine lsm_betr_plant_soilbgc_send_interface

  end interface
end module PlantSoilBGCMod
