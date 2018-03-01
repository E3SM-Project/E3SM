module TemperatureType

  !------------------------------------------------------------------------------
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use shr_log_mod     , only : errMsg => shr_log_errMsg
  use shr_infnan_mod  , only : nan => shr_infnan_nan, assignment(=)
  use decompMod       , only : bounds_type
  use abortutils      , only : endrun

  implicit none
  save
  private

!----------------------------------------------------
! column energy state variables structure
!----------------------------------------------------
  type, public :: temperature_type
    real(r8), pointer :: t_soisno_col(:,:)   => null()      !soil temperature (Kelvin)  (-nlevsno+1:nlevgrnd)
    real(r8), pointer :: t_soi10cm_col(:)     => null()    !soil temperature in top 10cm of soil (Kelvin)
    real(r8), pointer :: t_veg_patch              (:) => null()  ! patch vegetation temperature (Kelvin)
  contains
     procedure, public  :: Init
     procedure, private :: InitAllocate

  end type temperature_type

  contains


  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !ARGUMENTS:
    class(temperature_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg
    integer :: lbj,  ubj
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    lbj  = bounds%lbj;  ubj = bounds%ubj
    begp = bounds%begp; endp = bounds%endp
    allocate(this%t_soisno_col(begc:endc, lbj:ubj));  this%t_soisno_col(:,:) = nan
    allocate(this%t_soi10cm_col(begc:endc))           ;  this%t_soi10cm_col(:) = nan
    allocate(this%t_veg_patch  (begp:endp))        ; this%t_veg_patch              (:)   = nan
  end subroutine InitAllocate
end module TemperatureType
