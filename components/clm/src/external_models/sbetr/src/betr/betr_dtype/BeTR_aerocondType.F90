module BeTR_aerocondType
  !DESCRIPTION
  !module for aerodynamic calculation
  use bshr_kind_mod   , only : r8 => shr_kind_r8
  use bshr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
  use BeTR_decompMod  , only : bounds_type  => betr_bounds_type
  use betr_varcon     , only : spval => bspval, ispval => bispval

  implicit none

  private

  character(len=*), private, parameter :: mod_filename = &
       __FILE__
!
! !PUBLIC DATA:
!
  type, public :: betr_aerecond_type
   real(r8), pointer :: tempavg_agnpp_patch  (:) => null() !  temporary average above-ground NPP (gC/m2/s)
   real(r8), pointer :: annavg_agnpp_patch   (:) => null() !  annual average above-ground NPP (gC/m2/s)
   real(r8), pointer :: tempavg_bgnpp_patch  (:)=> null() !  temporary average below-ground NPP (gC/m2/s)
   real(r8), pointer :: annavg_bgnpp_patch   (:) => null() !  annual average below-ground NPP (gC/m2/s)
   real(r8), pointer :: plant_frootsc_patch  (:) => null() !
  contains
     procedure , public  :: Init
     procedure , private :: InitAllocate
  end type betr_aerecond_type

 contains

  subroutine Init(this, bounds)

  implicit none
  class(betr_aerecond_type) :: this
  type(bounds_type), intent(in) :: bounds

  call this%InitAllocate(bounds)
  end subroutine Init
  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)

  !Description
  !allocate memory
  implicit none
  class(betr_aerecond_type) :: this
  type(bounds_type), intent(in) :: bounds
  integer :: begp, endp

  begp = bounds%begp; endp=bounds%endp

  allocate(this%plant_frootsc_patch (begp:endp)); this%plant_frootsc_patch (:)   = nan
  allocate(this%annavg_agnpp_patch  (begp:endp)); this%annavg_agnpp_patch  (:) = spval ! To detect first year
  allocate(this%annavg_bgnpp_patch  (begp:endp)); this%annavg_bgnpp_patch  (:) = spval ! To detect first year
  allocate(this%tempavg_agnpp_patch (begp:endp)); this%tempavg_agnpp_patch (:) = spval
  allocate(this%tempavg_bgnpp_patch (begp:endp)); this%tempavg_bgnpp_patch (:) = spval

  end subroutine InitAllocate
end module BeTR_aerocondType
