module BeTR_phosphorusfluxRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_phosphorusflux_recv_type
    real(r8), pointer :: sminp_leached_col(:) => null()
    real(r8), pointer :: sminp_qdrain_col(:) => null()
    real(r8), pointer :: sminp_runoff_col(:) => null()
    real(r8), pointer :: som_p_leached_col(:) => null()
    real(r8), pointer :: som_p_runoff_col(:) => null()
    real(r8), pointer :: som_p_qdrain_col(:) => null()
    real(r8), pointer :: sminp_to_plant_patch(:) => null()   !integrated phosphate goes to plant at patch (gN/m2/s), will be summarized within the bgc model
    real(r8), pointer :: fire_decomp_ploss_vr_col(:,:) => null()  !will be summarized within the bgc model
    real(r8), pointer :: fire_decomp_ploss_col(:) => null()  !will be summarized within the bgc model
    real(r8), pointer :: supplement_to_sminp_col(:) => null()
    real(r8), pointer :: secondp_to_occlp_col(:) => null()
    real(r8), pointer :: supplement_to_sminp_vr_col(:,:) => null()
    real(r8), pointer :: secondp_to_occlp_vr_col(:,:) => null()
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
    procedure, public  :: summary
  end type betr_phosphorusflux_recv_type

 contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_phosphorusflux_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_phosphorusflux_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%sminp_leached_col(begc:endc))
  allocate(this%sminp_qdrain_col(begc:endc))
  allocate(this%sminp_runoff_col(begc:endc))
  allocate(this%som_p_leached_col(begc:endc))
  allocate(this%som_p_qdrain_col(begc:endc))
  allocate(this%som_p_runoff_col(begc:endc))

  allocate(this%sminp_to_plant_patch(begp:endp))
  allocate(this%fire_decomp_ploss_col(begc:endc))
  allocate(this%supplement_to_sminp_col(begc:endc))
  allocate(this%secondp_to_occlp_col(begc:endc))
  allocate(this%supplement_to_sminp_vr_col(begc:endc, lbj:ubj))
  allocate(this%secondp_to_occlp_vr_col(begc:endc, lbj:ubj))
  allocate(this%fire_decomp_ploss_vr_col(begc:endc,lbj:ubj))
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_phosphorusflux_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%sminp_leached_col(:) = value_column
  this%sminp_runoff_col(:) = value_column
  this%sminp_qdrain_col(:) = value_column
  this%som_p_leached_col(:) = value_column
  this%som_p_runoff_col(:) = value_column
  this%som_p_qdrain_col(:) = value_column

  this%supplement_to_sminp_vr_col(:,:) = value_column
  this%secondp_to_occlp_vr_col(:,:) = value_column
  this%fire_decomp_ploss_vr_col(:,:) = value_column
  end subroutine reset

  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz)

  implicit none
  class(betr_phosphorusflux_recv_type),intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)
  integer :: c, j
  this%supplement_to_sminp_col(:) = 0._r8
  this%secondp_to_occlp_col(:) = 0._r8
  this%fire_decomp_ploss_col(:) = 0._r8
  do j = lbj, ubj
    do c = bounds%begc, bounds%endc
      this%supplement_to_sminp_col(c) = this%supplement_to_sminp_col(c) + dz(c,j) * this%supplement_to_sminp_vr_col(c,j)
      this%secondp_to_occlp_col(c) = this%secondp_to_occlp_col(c) + dz(c,j) * this%secondp_to_occlp_vr_col(c,j)
      this%fire_decomp_ploss_col(c) = this%fire_decomp_ploss_col(c) + dz(c,j) * this%fire_decomp_ploss_vr_col(c,j)
    enddo
  enddo

  end subroutine summary
end module BeTR_phosphorusfluxRecvType
