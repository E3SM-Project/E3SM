module PlantNutKineticsMod

!
! DESCRIPTION
! transient plant kinetic parameters for nutrient compettion
!
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
  use betr_decompMod      , only : betr_bounds_type
implicit none

  type, public :: PlantNutKinetics_type
    real(r8), pointer :: plant_nh4_vmax_vr_patch(:,:)
    real(r8), pointer :: plant_no3_vmax_vr_patch(:,:)
    real(r8), pointer :: plant_p_vmax_vr_patch(:,:)
    real(r8), pointer :: plant_nh4_km_vr_patch(:,:)
    real(r8), pointer :: plant_no3_km_vr_patch(:,:)
    real(r8), pointer :: plant_p_km_vr_patch(:,:)
    real(r8), pointer :: plant_eff_ncompet_b_vr_patch(:,:)
    real(r8), pointer :: plant_eff_pcompet_b_vr_patch(:,:)
    real(r8), pointer :: minsurf_p_compet_vr_col(:,:)
    real(r8), pointer :: minsurf_nh4_compet_vr_col(:,:)

    !the following is only for eca-cnp bgc, for other applications
    !some of the parameters will be read through the interface
    real(r8), pointer :: km_minsurf_p_vr_col(:,:)
    real(r8), pointer :: km_decomp_nh4_vr_col(:,:)
    real(r8), pointer :: km_decomp_no3_vr_col(:,:)
    real(r8), pointer :: km_decomp_p_vr_col(:,:)
    real(r8), pointer :: km_nit_nh4_vr_col(:,:)
    real(r8), pointer :: km_den_no3_vr_col(:,:)
    real(r8), pointer :: decomp_eff_ncompet_b_vr_col(:,:)
    real(r8), pointer :: decomp_eff_pcompet_b_vr_col(:,:)
    real(r8), pointer :: den_eff_ncompet_b_vr_col(:,:)
    real(r8), pointer :: nit_eff_ncompet_b_vr_col(:,:)
  contains
    procedure, public  :: Init
    procedure, public  :: InitAllocate
    procedure, public  :: InitCold

  end type PlantNutKinetics_type

  !------------------------------------------------------------------------
  contains
    !------------------------------------------------------------------------
    subroutine Init(this, bounds)

     class(PlantNutKinetics_type) :: this
     type(betr_bounds_type), intent(in) :: bounds

     call this%InitAllocate ( bounds)
     call this%InitCold (bounds )
    end subroutine Init
    !------------------------------------------------------------------------
    subroutine InitAllocate(this, bounds)

     class(PlantNutKinetics_type) :: this
     type(betr_bounds_type), intent(in) :: bounds
     integer :: begp, endp, begc, endc
     integer :: lbj, ubj

     begp = bounds%begp; endp=bounds%endp
     begc = bounds%begc; endc=bounds%endc
     lbj = bounds%lbj; ubj=bounds%ubj
     allocate(this%plant_nh4_vmax_vr_patch(begp:endp, 1:ubj)); this%plant_nh4_vmax_vr_patch(:,:) = nan
     allocate(this%plant_no3_vmax_vr_patch(begp:endp, 1:ubj)); this%plant_no3_vmax_vr_patch(:,:) = nan
     allocate(this%plant_p_vmax_vr_patch(begp:endp, 1:ubj)); this%plant_p_vmax_vr_patch(:,:) = nan
     allocate(this%plant_no3_km_vr_patch(begp:endp, 1:ubj)); this%plant_no3_km_vr_patch(:,:) = nan
     allocate(this%plant_nh4_km_vr_patch(begp:endp, 1:ubj)); this%plant_nh4_km_vr_patch(:,:) = nan
     allocate(this%plant_p_km_vr_patch(begp:endp, 1:ubj)); this%plant_p_km_vr_patch(:,:) = nan

     allocate(this%plant_eff_ncompet_b_vr_patch(begp:endp,1:ubj)); this%plant_eff_ncompet_b_vr_patch(:,:) = nan
     allocate(this%plant_eff_pcompet_b_vr_patch(begp:endp,1:ubj)); this%plant_eff_pcompet_b_vr_patch(:,:) = nan

     allocate(this%minsurf_nh4_compet_vr_col(begc:endc,1:ubj)); this%minsurf_nh4_compet_vr_col(:,:) = nan
     allocate(this%minsurf_p_compet_vr_col(begc:endc,1:ubj)); this%minsurf_p_compet_vr_col(:,:) = nan
     allocate(this%km_minsurf_p_vr_col(begc:endc,1:ubj)); this%km_minsurf_p_vr_col(:,:) = nan
     allocate(this%km_decomp_nh4_vr_col(begc:endc, 1:ubj)); this%km_decomp_nh4_vr_col(:,:) = nan
     allocate(this%km_decomp_no3_vr_col(begc:endc,1:ubj)); this%km_decomp_no3_vr_col(:,:) = nan
     allocate(this%km_decomp_p_vr_col(begc:endc,1:ubj)); this%km_decomp_p_vr_col(:,:) = nan
     allocate(this%km_nit_nh4_vr_col(begc:endc,1:ubj)); this%km_nit_nh4_vr_col(:,:) = nan
     allocate(this%km_den_no3_vr_col(begc:endc,1:ubj)); this%km_den_no3_vr_col(:,:) = nan
     allocate(this%decomp_eff_ncompet_b_vr_col(begc:endc,1:ubj)); this%decomp_eff_ncompet_b_vr_col(:,:)=nan
     allocate(this%decomp_eff_pcompet_b_vr_col(begc:endc,1:ubj)); this%decomp_eff_pcompet_b_vr_col(:,:)=nan
     allocate(this%den_eff_ncompet_b_vr_col(begc:endc,1:ubj)); this%den_eff_ncompet_b_vr_col(:,:)=nan
     allocate(this%nit_eff_ncompet_b_vr_col(begc:endc,1:ubj)); this%nit_eff_ncompet_b_vr_col(:,:)=nan

    end subroutine InitAllocate
    !------------------------------------------------------------------------
    subroutine InitCold(this, bounds)

     class(PlantNutKinetics_type) :: this
     type(betr_bounds_type), intent(in) :: bounds


    end subroutine InitCold


end module PlantNutKineticsMod
