module BeTR_biogeoStateType

  use bshr_kind_mod   , only : r8 => shr_kind_r8
  use betr_decompMod  , only : betr_bounds_type
  use BeTR_carbonstateRecvType, only : betr_carbonstate_recv_type
  use BeTR_nitrogenstateRecvType, only : betr_nitrogenstate_recv_type
  use BeTR_phosphorusstateRecvType, only : betr_phosphorusstate_recv_type
  use tracer_varcon, only : use_c13_betr, use_c14_betr
implicit none

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  type betr_biogeo_state_type
    real(r8), pointer :: zwts_col           (:)   => null() ! the shallower between zwt_perch and zwt

    type(betr_carbonstate_recv_type) :: c12state_vars
    type(betr_carbonstate_recv_type) :: c13state_vars
    type(betr_carbonstate_recv_type) :: c14state_vars
    type(betr_nitrogenstate_recv_type):: n14state_vars
    type(betr_phosphorusstate_recv_type) :: p31state_vars
 contains
      procedure, public  :: Init
      procedure, private :: InitAllocate
      procedure, public  :: summary
      procedure, public  :: reset
  end type betr_biogeo_state_type

  public :: create_betr_biogeo_state

contains

  function create_betr_biogeo_state()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(betr_biogeo_state_type), pointer :: create_betr_biogeo_state
    class(betr_biogeo_state_type), pointer :: biogeo_state

    allocate(biogeo_state)
    create_betr_biogeo_state => biogeo_state

  end function create_betr_biogeo_state

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, active_soibgc)

  class(betr_biogeo_state_type)      :: this
  type(betr_bounds_type), intent(in) :: bounds
  logical, intent(in) :: active_soibgc

  if(active_soibgc)then
    call this%c12state_vars%Init(bounds)
    if(use_c13_betr)then
      call this%c13state_vars%Init(bounds)
    endif
    if(use_c14_betr)then
      call this%c14state_vars%Init(bounds)
    endif
    call this%n14state_vars%Init(bounds)
    call this%p31state_vars%Init(bounds)
  endif
  call this%InitAllocate(bounds)
  end subroutine Init
  !------------------------------------------------------------------------
  subroutine reset(this, value_column, active_soibgc)
  implicit none
  class(betr_biogeo_state_type)  :: this
  real(r8), intent(in) :: value_column
  logical, intent(in) :: active_soibgc

  if(active_soibgc)then
    call this%c12state_vars%reset(value_column)
    if(use_c13_betr)then
      call this%c13state_vars%reset(value_column)
    endif
    if(use_c14_betr)then
      call this%c14state_vars%reset(value_column)
    endif
    call this%n14state_vars%reset(value_column)
    call this%p31state_vars%reset(value_column)
  endif
  end subroutine reset
  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  class(betr_biogeo_state_type)      :: this
  type(betr_bounds_type), intent(in) :: bounds

  integer :: begp, endp
  integer :: begc, endc
  integer :: lbj, ubj
  begp = bounds%begp ; endp=bounds%endp
  begc = bounds%begc ; endc=bounds%endc
  lbj = bounds%lbj   ; ubj=bounds%ubj

  !soilhydrology
  allocate(this%zwts_col           (begc:endc) ) ! the shallower between zwt_perch and zwt

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine summary(this, bounds, lbj, ubj, dz, zs, active_soibgc)

  implicit none
  class(betr_biogeo_state_type),intent(inout)  :: this
  type(betr_bounds_type), intent(in) :: bounds
  integer , intent(in) :: lbj, ubj
  real(r8), intent(in) :: dz(bounds%begc:bounds%endc,lbj:ubj)
  real(r8), intent(in) :: zs(bounds%begc:bounds%endc,lbj:ubj)
  logical,  intent(in) :: active_soibgc

  if(active_soibgc)then
    call this%c12state_vars%summary(bounds, lbj, ubj, dz(bounds%begc:bounds%endc,lbj:ubj), zs(bounds%begc:bounds%endc,lbj:ubj))
    if(use_c13_betr)then
      call this%c13state_vars%summary(bounds, lbj, ubj, dz(bounds%begc:bounds%endc,lbj:ubj), zs(bounds%begc:bounds%endc,lbj:ubj))
    endif
    if(use_c14_betr)then
      call this%c14state_vars%summary(bounds, lbj, ubj, dz(bounds%begc:bounds%endc,lbj:ubj),zs(bounds%begc:bounds%endc,lbj:ubj))
    endif

    call this%n14state_vars%summary(bounds, lbj, ubj, dz(bounds%begc:bounds%endc,lbj:ubj), zs(bounds%begc:bounds%endc,lbj:ubj))
    call this%p31state_vars%summary(bounds, lbj, ubj, dz(bounds%begc:bounds%endc,lbj:ubj), zs(bounds%begc:bounds%endc,lbj:ubj))
  endif
  end subroutine summary
end module BeTR_biogeoStateType
