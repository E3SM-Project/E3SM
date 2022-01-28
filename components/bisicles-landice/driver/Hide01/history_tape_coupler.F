module history_tape_coupler

  ! Defines a class for controlling history frequency based on the coupler's history
  ! frequency.

  use history_tape_base, only : history_tape_base_type

  implicit none
  private
  save

  public :: history_tape_coupler_type
  type, extends(history_tape_base_type) :: history_tape_coupler_type
     private
   contains
     ! Logical function saying whether it's time to write a history file
     procedure :: is_time_to_write_hist

     ! Function returning a string describing the history frequency
     procedure :: history_frequency_string
  end type history_tape_coupler_type

  interface history_tape_coupler_type
     module procedure constructor
  end interface history_tape_coupler_type

contains

  !-----------------------------------------------------------------------
  function constructor(history_vars)
    !
    ! !DESCRIPTION:
    ! Creates a history_tape_coupler_type object
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(history_tape_coupler_type) :: constructor  ! function result

    ! List of variables to write to file
    character(len=*), intent(in) :: history_vars

    !-----------------------------------------------------------------------
  
    call constructor%set_history_vars(history_vars)
  end function constructor

  !-----------------------------------------------------------------------
  logical function is_time_to_write_hist(this, EClock)
    !
    ! !DESCRIPTION:
    ! Returns true if it is time to write the history tape associated with this controller.
    !
    ! !USES:
    use esmf, only: ESMF_Clock
    use seq_timemgr_mod, only : seq_timemgr_HistoryAlarmIsOn
    !
    ! !ARGUMENTS:
    class(history_tape_coupler_type), intent(in) :: this
    type(ESMF_Clock), intent(in) :: EClock

    !-----------------------------------------------------------------------

    is_time_to_write_hist = seq_timemgr_HistoryAlarmIsOn(EClock)

  end function is_time_to_write_hist

  !-----------------------------------------------------------------------
  function history_frequency_string(this)
    !
    ! !DESCRIPTION:
    ! Returns a string representation of this history frequency
    !
    ! TODO(wjs, 2015-02-17) This needs to be implemented. It is currently difficult (or
    ! impossible) to extract the frequency information from the coupler. Hopefully this
    ! will become easier once the coupler implements the necessary functionality for
    ! metadata on its own history files.
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: history_frequency_string  ! function result
    class(history_tape_coupler_type), intent(in) :: this

    !-----------------------------------------------------------------------

    history_frequency_string = '(matches coupler history frequency)'

  end function history_frequency_string
  
end module history_tape_coupler
