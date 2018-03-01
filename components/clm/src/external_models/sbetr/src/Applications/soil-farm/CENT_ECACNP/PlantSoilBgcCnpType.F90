module PlantSoilBgcCnpType

#include "bshr_assert.h"
  use PlantSoilBGCMod , only : plant_soilbgc_type
  use bshr_kind_mod           , only : r8 => shr_kind_r8
  use bshr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use bshr_log_mod            , only : errMsg => shr_log_errMsg
  implicit none

  private

  character(len=*), private, parameter :: filename = &
       __FILE__
  public :: plant_soilbgc_cnp_type

  type, extends(plant_soilbgc_type) :: &
    plant_soilbgc_cnp_type

    real(r8),pointer :: rt_vr_col(:,:) => null()
    real(r8),pointer :: plant_root_exudates_c(:) => null()
    real(r8),pointer :: plant_root_exudates_n(:) => null()
    real(r8),pointer :: plant_root_exudates_p(:) => null()


    real(r8), pointer :: plant_minn_nh4_active_yield_flx_vr_patch  (:,:)  => null() !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minn_no3_active_yield_flx_vr_patch  (:,:)  => null() !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minp_active_yield_flx_vr_patch   (:,:)  => null() !column level mineral phosphorus yeild from soil bgc calculation
    real(r8), pointer :: plant_minn_nh4_active_yield_flx_vr_col  (:,:)  => null() !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minn_no3_active_yield_flx_vr_col  (:,:)  => null() !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minp_active_yield_flx_vr_col(:,:) => null()
    real(r8), pointer :: plant_minn_active_yield_flx_col(:) => null()
    real(r8), pointer :: plant_minp_active_yield_flx_col(:) => null()
    real(r8), pointer :: plant_minn_nh4_passive_yield_flx_vr_patch  (:,:)  => null() !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minn_no3_passive_yield_flx_vr_patch  (:,:)  => null() !patch level mineral nitrogen yeild from soil bgc calculation
    real(r8), pointer :: plant_minp_passive_yield_flx_vr_patch   (:,:)  => null() !column level mineral phosphorus yeild from soil bgc calculation

  contains
    procedure :: Init_plant_soilbgc
    procedure :: plant_soilbgc_summary
    procedure :: integrate_vr_flux
    procedure :: lsm_betr_plant_soilbgc_recv
    procedure :: lsm_betr_plant_soilbgc_send
    procedure, private :: set_profiles_vars
    procedure, private :: InitAllocate
  end type plant_soilbgc_cnp_type

  interface plant_soilbgc_cnp_type
    module procedure constructor
  end interface plant_soilbgc_cnp_type

  contains

  !-------------------------------------------------------------------------------
  type(plant_soilbgc_cnp_type) function constructor()
  !
  ! !DESCRIPTION:
  ! create an object of type plant_soilbgc_cnp_type.
  ! Right now it is purposely empty
   type(plant_soilbgc_cnp_type), allocatable :: plants
   allocate(plants)
   constructor = plants
  end function constructor

  !-------------------------------------------------------------------------------
  subroutine Init_plant_soilbgc(this, bounds, lbj, ubj, namelist_buffer)

  !
  ! !DESCRIPTION:
  ! template for init_betrbgc
  !
  ! here I call alm instances directly?
  ! !USES:
  use BeTR_decompMod       , only : betr_bounds_type
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_cnp_type) , intent(inout) :: this
  type(betr_bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: lbj, ubj
  character(len=*)          , intent(in) :: namelist_buffer

  call this%InitAllocate(bounds, lbj, ubj)
  end subroutine Init_plant_soilbgc
  !----------------------------------------------------------------------

  subroutine InitAllocate(this, bounds, lbj, ubj)
  !
  !DESCRIPTION
  !allocate memories
  use BeTR_decompMod       , only : betr_bounds_type
  use betr_varcon         , only : betr_maxpatch_pft
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_cnp_type) , intent(inout) :: this
  type(betr_bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: lbj, ubj

  integer :: begc, endc
  integer :: begp, endp

  begc = bounds%begc; endc=bounds%endc
  begp = bounds%begp; endp=bounds%endp

  allocate(this%rt_vr_col(begc:endc,1:ubj)); this%rt_vr_col(:,:) = 0._r8


  allocate(this%plant_minn_nh4_active_yield_flx_vr_patch  (begp:endp,1:ubj)) !patch level mineral nitrogen yeild from soil bgc calculation
  allocate(this%plant_minn_no3_active_yield_flx_vr_patch  (begp:endp,1:ubj)) !patch level mineral nitrogen yeild from soil bgc calculation
  allocate(this%plant_minp_active_yield_flx_vr_patch  (begp:endp,1:ubj)) !column level mineral phosphorus yeild from soil bgc calculation

  allocate(this%plant_minn_nh4_passive_yield_flx_vr_patch  (begp:endp,1:ubj)) !patch level mineral nitrogen yeild from soil bgc calculation
  allocate(this%plant_minn_no3_passive_yield_flx_vr_patch  (begp:endp,1:ubj)) !patch level mineral nitrogen yeild from soil bgc calculation
  allocate(this%plant_minp_passive_yield_flx_vr_patch  (begp:endp,1:ubj)) !column level mineral phosphorus yeild from soil bgc calculation

  allocate(this%plant_minn_nh4_active_yield_flx_vr_col  (begc:endc,1:ubj)) !patch level mineral nitrogen yeild from soil bgc calculation
  allocate(this%plant_minn_no3_active_yield_flx_vr_col  (begc:endc,1:ubj)) !patch level mineral nitrogen yeild from soil bgc calculation
  allocate(this%plant_minp_active_yield_flx_vr_col  (begc:endc,1:ubj)) !column level mineral phosphorus yeild from soil bgc calculation

  allocate(this%plant_minn_active_yield_flx_col(begc:endc))
  allocate(this%plant_minp_active_yield_flx_col(begc:endc))
  end subroutine InitAllocate
  !----------------------------------------------------------------------
  subroutine plant_soilbgc_summary(this, bounds, lbj, ubj, pft, numf, &
       filter, dtime, dz, betrtracer_vars, tracerflux_vars,biogeo_flux, betr_status)

  ! !USES:
  use BeTRTracerType        , only : BeTRtracer_type
  use tracerfluxType        , only : tracerflux_type
  use BeTR_decompMod        , only : betr_bounds_type
  use bshr_kind_mod         , only : r8 => shr_kind_r8
  use BetrStatusType        , only : betr_status_type
  use BeTR_PatchType        , only : betr_patch_type
  use tracer_varcon         , only : natomw, patomw
  use BeTR_biogeoFluxType  , only : betr_biogeo_flux_type
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_cnp_type) , intent(inout) :: this
  type(betr_bounds_type)        , intent(in) :: bounds
  integer                   , intent(in) :: lbj, ubj
  type(betr_patch_type)     , intent(in) :: pft
  integer                   , intent(in) :: numf
  integer                   , intent(in) :: filter(:)
  real(r8)                  , intent(in) :: dz(bounds%begc: ,1: )
  real(r8)                  , intent(in) :: dtime
  type(BeTRtracer_type )    , intent(in) :: betrtracer_vars
  type(tracerflux_type)     , intent(in) :: tracerflux_vars
  type(betr_biogeo_flux_type) , intent(inout) :: biogeo_flux
  type(betr_status_type)    , intent(out):: betr_status

  !local variables
  integer :: p, c, fc

  call betr_status%reset()
  SHR_ASSERT_ALL((ubound(dz)==(/bounds%endc,ubj/)), errMsg(filename,__LINE__), betr_status)
  if(betr_status%check_status())return
  associate(                                                                &
    tracer_flx_vtrans_patch  => tracerflux_vars%tracer_flx_vtrans_patch   , &
    id_trc_no3x  => betrtracer_vars%id_trc_no3x, &
    id_trc_nh3x  => betrtracer_vars%id_trc_nh3x, &
    id_trc_p_sol => betrtracer_vars%id_trc_p_sol &

  )
  !now summarize all nutrient uptake for each plant patch
  do p = 1, pft%npfts
    c = pft%column(p)
    biogeo_flux%n14flux_vars%smin_nh4_to_plant_patch(p) =  &
      tracer_flx_vtrans_patch(p, id_trc_nh3x) * natomw / dtime + &
      dot_product(this%plant_minn_nh4_active_yield_flx_vr_patch(p,1:ubj), dz(c,1:ubj))

    biogeo_flux%n14flux_vars%smin_nh4_to_plant_patch(p) = biogeo_flux%n14flux_vars%smin_nh4_to_plant_patch(p)/pft%wtcol(p)

    biogeo_flux%n14flux_vars%smin_no3_to_plant_patch(p) = &
      tracer_flx_vtrans_patch(p, id_trc_no3x) * natomw / dtime + &
      dot_product(this%plant_minn_no3_active_yield_flx_vr_patch(p,1:ubj), dz(c,1:ubj))

    biogeo_flux%n14flux_vars%smin_no3_to_plant_patch(p) = biogeo_flux%n14flux_vars%smin_no3_to_plant_patch(p)/pft%wtcol(p)

    biogeo_flux%p31flux_vars%sminp_to_plant_patch(p) =  &
      tracer_flx_vtrans_patch(p, id_trc_p_sol) * patomw / dtime + &
      dot_product(this%plant_minp_active_yield_flx_vr_patch(p,1:ubj), dz(c,1:ubj))

    biogeo_flux%p31flux_vars%sminp_to_plant_patch(p) = biogeo_flux%p31flux_vars%sminp_to_plant_patch(p)/pft%wtcol(p)

  enddo

  end associate
  end subroutine plant_soilbgc_summary


  !----------------------------------------------------------------------

  subroutine integrate_vr_flux(this, bounds, numf, filter)

  !DESCRIPTION
  !integrates depth resolved fluxes into surface fluxes
  use BeTR_decompMod       , only : betr_bounds_type
  implicit none
  ! !ARGUMENTS:

  class(plant_soilbgc_cnp_type) , intent(inout) :: this
  type(betr_bounds_type)         , intent(in) :: bounds
  integer                   , intent(in) :: numf
  integer                   , intent(in) :: filter(:)


  end subroutine integrate_vr_flux

  !----------------------------------------------------------------------

  subroutine lsm_betr_plant_soilbgc_recv(this, bounds, numf, filter, betr_pft, biogeo_fluxes)
  !DESCRIPTION
  !return plant nutrient yield
  !
  !USES
  use BeTR_decompMod       , only : betr_bounds_type
  use BeTR_biogeoStateType , only : betr_biogeo_state_type
  use BeTR_biogeoFluxType  , only : betr_biogeo_flux_type
  use BeTR_PatchType, only : betr_patch_type
  implicit none
  ! !ARGUMENTS:

  class(plant_soilbgc_cnp_type) , intent(inout) :: this
  type(betr_bounds_type)      , intent(in)    :: bounds
  integer                     , intent(in)    :: numf
  integer                     , intent(in)    :: filter(:)
  type(betr_patch_type) , intent(in) :: betr_pft
  type(betr_biogeo_flux_type) , intent(inout) :: biogeo_fluxes


  end subroutine lsm_betr_plant_soilbgc_recv
  !----------------------------------------------------------------------
  subroutine lsm_betr_plant_soilbgc_send(this, bounds, numf, &
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
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_cnp_type) , intent(inout) :: this
  type(betr_bounds_type)       , intent(in) :: bounds
  integer                      , intent(in) :: numf
  integer                      , intent(in) :: filter(:)
  type(betr_patch_type) , intent(in) :: betr_pft
  type(betr_biogeophys_input_type), intent(in):: biogeo_forc
  type(betr_biogeo_state_type) , intent(in) :: biogeo_states
  type(betr_biogeo_flux_type)  , intent(in) :: biogeo_fluxes

  integer :: p, c

  !set flux profiles, e.g. root respiration
  call this%set_profiles_vars(bounds, numf, filter, betr_pft, biogeo_forc, biogeo_fluxes)

  !set root exudation, which will be added in the future.

  end subroutine lsm_betr_plant_soilbgc_send

  !----------------------------------------------------------------------
  subroutine set_profiles_vars(this, bounds, numf, filter, betr_pft, biogeo_forc, biogeo_fluxes)
  !
  !DESCRIPTION
  !this setup root respiration profiles
  use BeTR_biogeoFluxType  , only : betr_biogeo_flux_type
  use BeTR_decompMod       , only : betr_bounds_type
  use BeTR_biogeophysInputType , only : betr_biogeophys_input_type
  use BeTR_PatchType, only : betr_patch_type
  implicit none
  ! !ARGUMENTS:
  class(plant_soilbgc_cnp_type) , intent(inout) :: this
  type(betr_bounds_type)       , intent(in) :: bounds
  integer                      , intent(in) :: numf
  integer                      , intent(in) :: filter(:)
  type(betr_patch_type) , intent(in) :: betr_pft
  type(betr_biogeophys_input_type), intent(in):: biogeo_forc
  type(betr_biogeo_flux_type)  , intent(in) :: biogeo_fluxes

  integer :: j, p, c

  do j = 1, bounds%ubj
    this%rt_vr_col(:,j) = 0._r8
    do p = 1, betr_pft%npfts
      c = betr_pft%column(p)
      this%rt_vr_col(c,j)  = this%rt_vr_col(c,j) + biogeo_forc%rr_patch(p,j) * betr_pft%wtcol(p) !gC/m2/s
    enddo
  enddo

  end subroutine set_profiles_vars
  !----------------------------------------------------------------------
end module PlantSoilBgcCnpType
