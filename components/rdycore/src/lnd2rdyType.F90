module lnd2rdyType

  use shr_kind_mod  , only : r8 => shr_kind_r8
  use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod   , only : errMsg => shr_log_errMsg
  use rdydecompMod  , only : rdy_bounds_type

  implicit none
  private
  save

  type, public :: lnd2rdy_type

     real(r8), pointer :: forc_qsur(:) => null() ! liquid surface runoff    [m/s]
     real(r8), pointer :: forc_qsub(:) => null() ! liquid subsurface runoff [m/s]

   contains

     procedure, public :: Init
     procedure, public :: Destroy

  end type lnd2rdy_type

contains

  !------------------------------------------------------------------------
  subroutine Init(this, rdy_bounds)
    !
    ! !DESCRIPTION:
    ! Allocates memory
    !
    implicit none
    !
    ! !ARGUMENTS:
    class(lnd2rdy_type)   :: this
    type(rdy_bounds_type) :: rdy_bounds
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival = 0._r8 ! initial value
    integer  :: begg, endg   ! indices

    begg = rdy_bounds%begg
    endg = rdy_bounds%endg

    allocate(this%forc_qsur(begg:endg)); this%forc_qsur(:) = ival
    allocate(this%forc_qsub(begg:endg)); this%forc_qsub(:) = ival

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine Destroy(this)
    !
    ! !DESCRIPTION:
    ! Frees up memory
    !
    implicit none
    !
    class(lnd2rdy_type) :: this
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ival = 0._r8 ! initial value

    deallocate(this%forc_qsur)
    deallocate(this%forc_qsub)

  end subroutine Destroy

end module lnd2rdyType
