module BeTR_nitrogenstateRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_nitrogenstate_recv_type
    real(r8), pointer :: cwdn_col(:) => null()
    real(r8), pointer :: totlitn_col(:) => null()
    real(r8), pointer :: totsomn_col(:) => null()
    real(r8), pointer :: sminn_col(:) => null()
    real(r8), pointer :: sminn_nh4_col(:)=> null()
    real(r8), pointer :: sminn_no3_col(:)=> null()
    real(r8), pointer :: totlitn_1m_col(:) => null()
    real(r8), pointer :: totsomn_1m_col(:) => null()
    real(r8), pointer :: cwdn_vr_col(:,:) => null()
    real(r8), pointer :: totlitn_vr_col(:,:) => null()
    real(r8), pointer :: totsomn_vr_col(:,:) => null()
    real(r8), pointer :: sminn_vr_col(:,:) => null()
    real(r8), pointer :: sminn_nh4_vr_col(:,:) => null()
    real(r8), pointer :: sminn_no3_vr_col(:,:) => null()

  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
    procedure, public  :: summary
  end type betr_nitrogenstate_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_nitrogenstate_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
  implicit none
  class(betr_nitrogenstate_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%cwdn_col(begc:endc))
  allocate(this%totlitn_col(begc:endc))
  allocate(this%totsomn_col(begc:endc))
  allocate(this%sminn_col(begc:endc))
  allocate(this%sminn_nh4_col(begc:endc))
  allocate(this%sminn_no3_col(begc:endc))
  allocate(this%totlitn_1m_col(begc:endc))
  allocate(this%totsomn_1m_col(begc:endc))

  allocate(this%cwdn_vr_col(begc:endc,lbj:ubj)); this%cwdn_vr_col(:,:) = nan
  allocate(this%totlitn_vr_col(begc:endc,lbj:ubj)); this%totlitn_vr_col(:,:)=nan
  allocate(this%totsomn_vr_col(begc:endc,lbj:ubj)); this%totsomn_vr_col(:,:) =nan
  allocate(this%sminn_vr_col(begc:endc,lbj:ubj)); this%sminn_vr_col(:,:) =nan
  allocate(this%sminn_nh4_vr_col(begc:endc,lbj:ubj)); this%sminn_nh4_vr_col(:,:) =nan
  allocate(this%sminn_no3_vr_col(begc:endc,lbj:ubj)); this%sminn_no3_vr_col(:,:) =nan
  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_nitrogenstate_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%cwdn_vr_col(:,:) = value_column
  this%totlitn_vr_col(:,:) = value_column
  this%totsomn_vr_col(:,:) = value_column
  this%sminn_vr_col(:,:) = value_column
  this%sminn_nh4_vr_col(:,:) = value_column
  this%sminn_no3_vr_col(:,:) = value_column

  end subroutine reset

  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz, zs)
  use tracer_varcon, only : natomw
  implicit none
  class(betr_nitrogenstate_recv_type),intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)
  real(r8), intent(in) :: zs(bounds%begc:bounds%endc,lbj:ubj)


  integer :: c, j

  this%cwdn_col(:) = 0._r8
  this%totlitn_col(:) = 0._r8
  this%totsomn_col(:) = 0._r8
  this%sminn_col(:) = 0._r8
  this%sminn_nh4_col(:) = 0._r8
  this%sminn_no3_col(:) = 0._r8
  this%totlitn_1m_col(:) = 0._r8
  this%totsomn_1m_col(:) = 0._r8
  do j = lbj, ubj
    do c = bounds%begc, bounds%endc
      this%cwdn_col(c) = this%cwdn_col(c) + dz(c,j) * this%cwdn_vr_col(c,j)
      this%totlitn_col(c) = this%totlitn_col(c) + dz(c,j)*this%totlitn_vr_col(c,j)
      this%totsomn_col(c) = this%totsomn_col(c) + dz(c,j)*this%totsomn_vr_col(c,j)
      this%sminn_col(c) = this%sminn_col(c) + dz(c,j)*this%sminn_vr_col(c,j)
      this%sminn_nh4_col(c) = this%sminn_nh4_col(c) + dz(c,j)*this%sminn_nh4_vr_col(c,j)
      this%sminn_no3_col(c) = this%sminn_no3_col(c) + dz(c,j)*this%sminn_no3_vr_col(c,j)
      if(zs(c,j)<1._r8)then
        if(zs(c,j+1)>1._r8)then
          this%totlitn_1m_col(c) = this%totlitn_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totlitn_vr_col(c,j)
          this%totsomn_1m_col(c) = this%totsomn_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totsomn_vr_col(c,j)
        else
          this%totlitn_1m_col(c) = this%totlitn_1m_col(c) + dz(c,j)*this%totlitn_vr_col(c,j)
          this%totsomn_1m_col(c) = this%totsomn_1m_col(c) + dz(c,j)*this%totsomn_vr_col(c,j)
        endif
      endif
    enddo
  enddo
  end subroutine summary
end module BeTR_nitrogenstateRecvType
