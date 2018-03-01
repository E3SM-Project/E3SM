module BeTR_carbonstateRecvType
  use bshr_kind_mod  , only : r8 => shr_kind_r8
  use betr_decompMod , only : betr_bounds_type
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type, public :: betr_carbonstate_recv_type
     real(r8), pointer :: cwdc_col(:) => null()
     real(r8), pointer :: totlitc_col(:)  => null()
     real(r8), pointer :: totsomc_col(:) => null()
     real(r8), pointer :: totlitc_1m_col(:) => null()
     real(r8), pointer :: totsomc_1m_col(:) => null()
     real(r8), pointer :: cwdc_vr_col(:,:) => null()
     real(r8), pointer :: totlitc_vr_col(:,:) => null()
     real(r8), pointer :: totsomc_vr_col(:,:) => null()
 contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
    procedure, public  :: reset
    procedure, public  :: summary
  end type betr_carbonstate_recv_type

 contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

  implicit none
  class(betr_carbonstate_recv_type), intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
  use bshr_infnan_mod     , only : nan => shr_infnan_nan, assignment(=)
  implicit none
  class(betr_carbonstate_recv_type), intent(inout) :: this
  type(betr_bounds_type), intent(in) :: bounds
  !temporary variables
  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  allocate(this%cwdc_col(begc:endc))
  allocate(this%totlitc_col(begc:endc))
  allocate(this%totsomc_col(begc:endc))
  allocate(this%totlitc_1m_col(begc:endc))
  allocate(this%totsomc_1m_col(begc:endc))

  allocate(this%cwdc_vr_col(begc:endc,lbj:ubj)); this%cwdc_vr_col(:,:)   = nan
  allocate(this%totlitc_vr_col(begc:endc,lbj:ubj)); this%totlitc_vr_col(:,:) = nan
  allocate(this%totsomc_vr_col(begc:endc,lbj:ubj)); this%totsomc_vr_col(:,:) = nan

  end subroutine InitAllocate



  !------------------------------------------------------------------------
  subroutine reset(this, value_column)
  implicit none
  class(betr_carbonstate_recv_type)  :: this
  real(r8), intent(in) :: value_column

  this%cwdc_vr_col(:,:) = value_column
  this%totlitc_vr_col(:,:) = value_column
  this%totsomc_vr_col(:,:) = value_column
  
  
  end subroutine reset


  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz, zs)

  implicit none
  class(betr_carbonstate_recv_type),intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)
  real(r8), intent(in) :: zs(bounds%begc:bounds%endc,lbj:ubj)

  integer :: c, j

  this%cwdc_col(:) = 0._r8
  this%totlitc_col(:) = 0._r8
  this%totsomc_col(:) = 0._r8
  this%totlitc_1m_col(:) = 0._r8
  this%totsomc_1m_col(:) = 0._r8
  do j = lbj, ubj 
    do c = bounds%begc, bounds%endc
      this%cwdc_col(c) = this%cwdc_col(c) + dz(c,j) * this%cwdc_vr_col(c,j)
      this%totlitc_col(c) = this%totlitc_col(c) + dz(c,j)*this%totlitc_vr_col(c,j)
      this%totsomc_col(c) = this%totsomc_col(c) + dz(c,j)*this%totsomc_vr_col(c,j)

      if(zs(c,j)<1._r8)then
        if(zs(c,j+1)>1._r8)then
          this%totlitc_1m_col(c) = this%totlitc_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totlitc_vr_col(c,j)
          this%totsomc_1m_col(c) = this%totsomc_1m_col(c) + (dz(c,j)-(zs(c,j)-1._r8))*this%totsomc_vr_col(c,j)
        else
          this%totlitc_1m_col(c) = this%totlitc_1m_col(c) + dz(c,j)*this%totlitc_vr_col(c,j)
          this%totsomc_1m_col(c) = this%totsomc_1m_col(c) + dz(c,j)*this%totsomc_vr_col(c,j)
        endif
      endif
    enddo
  enddo
  end subroutine summary
end module BeTR_carbonstateRecvType
