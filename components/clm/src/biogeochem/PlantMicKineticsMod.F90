module PlantMicKineticsMod

!
! DESCRIPTION
! compute depth-dependent kinetic parameters used for nutrient competition
  use shr_kind_mod           , only : r8 => shr_kind_r8
  use shr_infnan_mod         , only : nan => shr_infnan_nan, assignment(=)
  use decompMod              , only : bounds_type
  use clm_varpar             , only : nlevdecomp_full, nlevgrnd, nlevdecomp
implicit none

  type, public :: PlantMicKinetics_type
    real(r8), pointer :: plant_nh4_vmax_vr_patch(:,:)
    real(r8), pointer :: plant_no3_vmax_vr_patch(:,:)
    real(r8), pointer :: plant_p_vmax_vr_patch(:,:)
    real(r8), pointer :: plant_nh4_km_vr_patch(:,:)
    real(r8), pointer :: plant_no3_km_vr_patch(:,:)
    real(r8), pointer :: plant_p_km_vr_patch(:,:)

    real(r8), pointer :: plant_eff_ncompet_b_vr_patch(:,:)
    real(r8), pointer :: plant_eff_pcompet_b_vr_patch(:,:)
    real(r8), pointer :: decomp_eff_ncompet_b_vr_col(:,:)
    real(r8), pointer :: decomp_eff_pcompet_b_vr_col(:,:)
    real(r8), pointer :: minsurf_p_compet_vr_col(:,:)

    real(r8), pointer :: vmax_minsurf_p_vr_col(:,:)
    real(r8), pointer :: km_minsurf_p_vr_col(:,:)
    real(r8), pointer :: km_decomp_nh4_vr_col(:,:)
    real(r8), pointer :: km_decomp_no3_vr_col(:,:)
    real(r8), pointer :: km_decomp_p_vr_col(:,:)
    real(r8), pointer :: km_nit_nh4_vr_col(:,:)
    real(r8), pointer :: km_den_no3_vr_col(:,:)

  contains
    procedure, public  :: Init
    procedure, public  :: InitAllocate
    procedure, public  :: InitCold

  end type PlantMicKinetics_type
  !------------------------------------------------------------------------
  contains
    !------------------------------------------------------------------------
    subroutine Init(this, bounds)

     class(PlantMicKinetics_type) :: this
     type(bounds_type), intent(in) :: bounds

     call this%InitAllocate ( bounds)
     call this%InitCold (bounds )
    end subroutine Init
    !------------------------------------------------------------------------
    subroutine InitAllocate(this, bounds)

     class(PlantMicKinetics_type) :: this
     type(bounds_type), intent(in) :: bounds
     integer :: begp, endp, begc, endc

     begp = bounds%begp; endp=bounds%endp
     begc = bounds%begc; endc=bounds%endc
     allocate(this%plant_nh4_vmax_vr_patch(begp:endp, 1:nlevdecomp_full)); this%plant_nh4_vmax_vr_patch(:,:) = nan
     allocate(this%plant_no3_vmax_vr_patch(begp:endp, 1:nlevdecomp_full)); this%plant_no3_vmax_vr_patch(:,:) = nan
     allocate(this%plant_p_vmax_vr_patch(begp:endp, 1:nlevdecomp_full)); this%plant_p_vmax_vr_patch(:,:) = nan

     allocate(this%plant_no3_km_vr_patch(begp:endp, 1:nlevdecomp_full)); this%plant_no3_km_vr_patch(:,:) = nan
     allocate(this%plant_nh4_km_vr_patch(begp:endp, 1:nlevdecomp_full)); this%plant_nh4_km_vr_patch(:,:) = nan
     allocate(this%plant_p_km_vr_patch(begp:endp, 1:nlevdecomp_full)); this%plant_p_km_vr_patch(:,:) = nan


     allocate(this%plant_eff_ncompet_b_vr_patch(begp:endp,1:nlevdecomp_full)); this%plant_eff_ncompet_b_vr_patch(:,:)=nan
     allocate(this%plant_eff_pcompet_b_vr_patch(begp:endp,1:nlevdecomp_full)); this%plant_eff_pcompet_b_vr_patch(:,:)=nan
     allocate(this%decomp_eff_ncompet_b_vr_col(begc:endc,1:nlevdecomp_full)); this%decomp_eff_ncompet_b_vr_col(:,:) = nan
     allocate(this%decomp_eff_pcompet_b_vr_col(begc:endc,1:nlevdecomp_full)); this%decomp_eff_pcompet_b_vr_col(:,:) = nan
     allocate(this%minsurf_p_compet_vr_col(begc:endc,1:nlevdecomp_full)); this%minsurf_p_compet_vr_col(:,:) = nan

     allocate(this%vmax_minsurf_p_vr_col(begc:endc, 1:nlevdecomp_full)); this%vmax_minsurf_p_vr_col(:,:) = nan
     allocate(this%km_minsurf_p_vr_col(begc:endc,1:nlevdecomp_full)); this%km_minsurf_p_vr_col(:,:) = nan
     allocate(this%km_decomp_nh4_vr_col(begc:endc, 1:nlevdecomp_full)); this%km_decomp_nh4_vr_col(:,:) = nan
     allocate(this%km_decomp_no3_vr_col(begc:endc,1:nlevdecomp_full)); this%km_decomp_no3_vr_col(:,:) = nan
     allocate(this%km_decomp_p_vr_col(begc:endc,1:nlevdecomp_full)); this%km_decomp_p_vr_col(:,:) = nan
     allocate(this%km_nit_nh4_vr_col(begc:endc,1:nlevdecomp_full)); this%km_nit_nh4_vr_col(:,:) = nan
     allocate(this%km_den_no3_vr_col(begc:endc,1:nlevdecomp_full)); this%km_den_no3_vr_col(:,:) = nan
    end subroutine InitAllocate
    !------------------------------------------------------------------------
    subroutine InitCold(this, bounds)

     class(PlantMicKinetics_type) :: this
     type(bounds_type), intent(in) :: bounds

    end subroutine InitCold

end module PlantMicKineticsMod
