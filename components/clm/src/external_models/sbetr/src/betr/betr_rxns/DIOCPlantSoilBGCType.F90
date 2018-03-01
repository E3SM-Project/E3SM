module DIOCPlantSoilBGCType
  !
  !DESCRIPTION
  ! mock interface for plant soil bgc coupling
#include "bshr_assert.h"  
  !USES
  use PlantSoilBGCMod , only : plant_soilbgc_type
  use betr_decompMod  , only : bounds_type => betr_bounds_type
  use bshr_log_mod    , only : errMsg => shr_log_errMsg
  implicit none

  private
  character(len=*), private, parameter :: filename = &
       __FILE__
  public :: plant_soilbgc_dioc_run_type

  type, extends(plant_soilbgc_type) :: &
    plant_soilbgc_dioc_run_type
  private
    contains
    procedure :: Init_plant_soilbgc
    procedure :: plant_soilbgc_summary
    procedure :: integrate_vr_flux
    procedure :: lsm_betr_plant_soilbgc_recv
    procedure :: lsm_betr_plant_soilbgc_send
  end type plant_soilbgc_dioc_run_type

  interface plant_soilbgc_dioc_run_type
    module procedure constructor
  end interface plant_soilbgc_dioc_run_type

  contains

  !-------------------------------------------------------------------------------
  type(plant_soilbgc_dioc_run_type) function constructor()
  !
  ! !DESCRIPTION:
  ! create an object of type plant_soilbgc_dioc_run_type.
  ! Right now it is purposely empty
    type(plant_soilbgc_dioc_run_type), allocatable :: plants
    allocate(plants)
    constructor = plants
  end function constructor

  !-------------------------------------------------------------------------------
  subroutine Init_plant_soilbgc(this, bounds, lbj, ubj, namelist_buffer)
  !
  ! !DESCRIPTION:
  ! template for init_betrbgc
  !
  ! !USES:
  use gbetrType      , only : gbetr_type
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_dioc_run_type) , intent(inout) :: this
  type(bounds_type)                  , intent(in) :: bounds
  integer                            , intent(in) :: lbj, ubj
  character(len=*)                   , intent(in) :: namelist_buffer

  ! remove compiler warnings for unused dummy args
  if (this%dummy_compiler_warning) continue
  if (bounds%begc > 0)             continue
  if (lbj > 0)                     continue
  if (ubj > 0)                     continue

  end subroutine Init_plant_soilbgc


  !----------------------------------------------------------------------
  subroutine plant_soilbgc_summary(this,bounds, lbj, ubj, pft, numf, &
       filter, dtime, dz, betrtracer_vars, tracerflux_vars, biogeo_flux, betr_status)
  !DESCRIPTION
  !summarize bgc coupling flux variables
  ! !USES:
  Use BeTRTracerType , only : BeTRtracer_type
  use tracerfluxType , only : tracerflux_type
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use BetrStatusType , only : betr_status_type
  use BeTR_PatchType , only : betr_patch_type
  use BeTR_biogeoFluxType  , only : betr_biogeo_flux_type
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_dioc_run_type) , intent(inout) :: this
  type(bounds_type)                  , intent(in) :: bounds
  integer                            , intent(in) :: lbj, ubj
  type(betr_patch_type)              , intent(in) :: pft
  integer                            , intent(in) :: numf
  integer                            , intent(in) :: filter(:)
  real(r8)                           , intent(in) :: dtime
  real(r8)                           , intent(in) :: dz(bounds%begc: ,1: )
  type(BeTRtracer_type )             , intent(in) :: betrtracer_vars
  type(tracerflux_type)              , intent(in) :: tracerflux_vars
  type(betr_biogeo_flux_type)        , intent(inout) :: biogeo_flux
  type(betr_status_type)             , intent(out):: betr_status

  call betr_status%reset()
  SHR_ASSERT_ALL((ubound(dz)==(/bounds%endc,ubj/)), errMsg(filename,__LINE__), betr_status)
  if(betr_status%check_status())return

  ! remove compiler warnings for unused dummy args
  if (this%dummy_compiler_warning)                       continue
  if (bounds%begc > 0)                                   continue
  if (numf > 0)                                          continue
  if (size(filter) > 0)                                  continue
  if (lbj > 0)                                           continue
  if (ubj > 0)                                           continue
  if (size(dz) > 0)                                      continue
  if (len(betrtracer_vars%betr_simname) > 0)             continue
  if (size(tracerflux_vars%tracer_flx_top_soil_col) > 0) continue

  end subroutine plant_soilbgc_summary


  !----------------------------------------------------------------------

  subroutine integrate_vr_flux(this, bounds, numf, filter)


  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_dioc_run_type) , intent(inout) :: this
  type(bounds_type)                  , intent(in) :: bounds
  integer                            , intent(in) :: numf
  integer                            , intent(in) :: filter(:)

  ! remove compiler warnings for unused dummy args
  if (this%dummy_compiler_warning) continue
  if (bounds%begc > 0)             continue
  if (numf > 0)                    continue
  if (size(filter) > 0)            continue

  end subroutine integrate_vr_flux

  !----------------------------------------------------------------------

  subroutine lsm_betr_plant_soilbgc_recv(this, bounds, numf, filter, betr_pft, biogeo_fluxes)

  !DESCRIPTION
  !return plant nutrient yield
  !
  !USES
  use BeTR_biogeoFluxType, only : betr_biogeo_flux_type
  use BeTR_PatchType, only : betr_patch_type
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_dioc_run_type) , intent(inout)    :: this
  type(bounds_type)                  , intent(in)    :: bounds
  integer                            , intent(in)    :: numf
  integer                            , intent(in)    :: filter(:)
  type(betr_patch_type) , intent(in) :: betr_pft
  type(betr_biogeo_flux_type)        , intent(inout) :: biogeo_fluxes

  ! remove compiler warnings for unused dummy args
  if (this%dummy_compiler_warning)        continue
  if (bounds%begc > 0)                    continue
  if (numf > 0)                           continue
  if (size(filter) > 0)                   continue
  if (size(biogeo_fluxes%qflx_adv_col)>0) continue
  end subroutine lsm_betr_plant_soilbgc_recv


  !----------------------------------------------------------------------

  subroutine lsm_betr_plant_soilbgc_send(this, bounds, numf, filter,  &
    betr_pft, biogeo_forc, biogeo_states, biogeo_fluxes)
  !
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
  class(plant_soilbgc_dioc_run_type) , intent(inout) :: this
  type(betr_bounds_type)             , intent(in) :: bounds
  integer                            , intent(in) :: numf
  integer                            , intent(in) :: filter(:)
  type(betr_patch_type)              , intent(in) :: betr_pft
  type(betr_biogeophys_input_type), intent(in):: biogeo_forc
  type(betr_biogeo_state_type)       , intent(in) :: biogeo_states
  type(betr_biogeo_flux_type)        , intent(in) :: biogeo_fluxes

  ! remove compiler warnings for unused dummy args
  if (this%dummy_compiler_warning)       continue
  if (bounds%begc > 0)                   continue
  if (numf > 0)                          continue
  if (size(filter) > 0)                  continue
  if (size(biogeo_states%zwts_col)>0)    continue
  if(size(biogeo_fluxes%qflx_adv_col)>0) continue

  end subroutine lsm_betr_plant_soilbgc_send
end module DIOCPlantSoilBGCType
