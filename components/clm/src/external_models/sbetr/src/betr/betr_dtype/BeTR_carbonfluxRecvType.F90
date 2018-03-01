module BeTR_carbonfluxRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_carbonflux_recv_type
    real(r8), pointer :: hr_col(:) => null()
    real(r8), pointer :: som_c_leached_col(:) => null()
    real(r8), pointer :: som_c_runoff_col(:) => null()
    real(r8), pointer :: som_c_qdrain_col(:) => null()
    real(r8), pointer :: hr_vr_col(:,:) => null()
    real(r8), pointer :: fire_decomp_closs_vr_col(:,:) => null()  !will be summarized from the specific bgc model
    real(r8), pointer :: fire_decomp_closs_col(:) => null()  !will be summarized from the specific bgc model
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
    procedure, public  :: summary
  end type betr_carbonflux_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_carbonflux_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  implicit none
  class(betr_carbonflux_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%hr_col(begc:endc))
  allocate(this%fire_decomp_closs_col(begc:endc))
  allocate(this%som_c_leached_col(begc:endc))
  allocate(this%som_c_runoff_col(begc:endc))
  allocate(this%som_c_qdrain_col(begc:endc))
  allocate(this%hr_vr_col(begc:endc,lbj:ubj))
  allocate(this%fire_decomp_closs_vr_col(begc:endc,lbj:ubj))
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_carbonflux_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%hr_vr_col(:,:) = value_column
  this%fire_decomp_closs_vr_col(:,:) = value_column
  this%som_c_leached_col(:) = value_column
  this%som_c_qdrain_col(:) = value_column
  this%som_c_runoff_col(:) = value_column

  end subroutine reset


  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz)

  implicit none
  class(betr_carbonflux_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)

  integer :: c, j

  this%hr_col(:) = 0._r8
  this%fire_decomp_closs_col(:) = 0._r8
  do j = lbj, ubj
    do c=bounds%begc, bounds%endc
      this%hr_col(c) = this%hr_col(c) + dz(c,j) * this%hr_vr_col(c,j)
      this%fire_decomp_closs_col(c) = this%fire_decomp_closs_col(c) + dz(c,j) * &
           this%fire_decomp_closs_vr_col(c,j)
    enddo
  enddo
  end subroutine summary
end module BeTR_carbonfluxRecvType
