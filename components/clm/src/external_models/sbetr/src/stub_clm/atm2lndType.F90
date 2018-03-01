module atm2lndType

! DESCRIPTION
! derived data for atmospheric/land exchange variables

  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type

implicit none

  type, public :: atm2lnd_type

  real(r8), pointer :: forc_pbot_downscaled_col      (:)   => null() ! downscaled atm pressure (Pa)
  real(r8), pointer :: forc_t_downscaled_col         (:)   => null() ! downscaled atm temperature (Kelvin)
  contains
    procedure, public  :: Init
    procedure, private :: InitAllocate
  end type atm2lnd_type
  !----------------------------------------------------

  contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds

    call this%InitAllocate(bounds)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    ! !ARGUMENTS:
    class(atm2lnd_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc
    integer :: begg, endg

    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc


    allocate(this%forc_pbot_downscaled_col(begc:endc));    this%forc_pbot_downscaled_col(:) = nan
    allocate(this%forc_t_downscaled_col(begc:endc));    this%forc_t_downscaled_col(:) = nan

  end subroutine InitAllocate
end module atm2lndType
