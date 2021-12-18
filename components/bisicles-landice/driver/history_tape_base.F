module history_tape_base

  ! This module defines an abstract base class to implement a single history tape.

  use shr_kind_mod, only : CXX => SHR_KIND_CXX

  implicit none
  private
  save

  ! max character length
  integer, parameter, public :: len_history_vars = CXX

  public :: history_tape_base_type
  type, abstract :: history_tape_base_type
     private

     ! Names of CISM variables to be output in cesm history files
     !
     ! COMPILER_BUG(wjs, 2015-02-19, pgi 14.7) Ideally, this would be an allocatable
     ! character variable. But when written that way, it gets filled with garbage by pgi
     ! 14.7. So for now, I'm using a declared maximum length together with a check that
     ! it's not being set to greater than its length.
     character(len=len_history_vars) :: history_vars

   contains
     ! ------------------------------------------------------------------------
     ! Public methods
     ! ------------------------------------------------------------------------
     procedure :: write_history    ! write history, if it's time to do so
     procedure :: set_history_vars ! set the list of history variables

     ! ------------------------------------------------------------------------
     ! The following are public simply because they need to be overridden by derived
     ! classes. They should not be called directly by clients.
     ! ------------------------------------------------------------------------
     ! Logical function saying whether it's time to write a history file
     procedure(is_time_to_write_hist_interface), deferred :: is_time_to_write_hist

     ! Function returning a string describing the history frequency
     procedure(history_frequency_string_interface), deferred :: history_frequency_string
  end type history_tape_base_type

  abstract interface
     
     logical function is_time_to_write_hist_interface(this, EClock)
       use esmf, only : ESMF_Clock
       import :: history_tape_base_type

       class(history_tape_base_type), intent(in) :: this
       type(ESMF_Clock), intent(in) :: EClock
     end function is_time_to_write_hist_interface

     function history_frequency_string_interface(this)
       import :: history_tape_base_type

       character(len=:), allocatable :: history_frequency_string_interface
       class(history_tape_base_type), intent(in) :: this
     end function history_frequency_string_interface

  end interface

contains

  !-----------------------------------------------------------------------
  subroutine write_history(this, instance, EClock, force_write)
    !
    ! !DESCRIPTION:
    ! Write a CISM history file, if it's time to do so.
    !
    ! This routine should be called every time step. It will return without doing
    ! anything if it isn't yet time to write a history file.
    !
    ! If force_write is present and true, then a history file is written regardless of
    ! the check for whether it's time to do so.
    !
    ! !USES:
    use glc_io, only : glc_io_write_history
    use glad_type, only : glad_instance
    use esmf, only: ESMF_Clock
    !
    ! !ARGUMENTS:
    class(history_tape_base_type), intent(in) :: this
    type(glad_instance), intent(inout) :: instance
    type(ESMF_Clock),     intent(in)    :: EClock
    logical, intent(in), optional :: force_write
    !
    ! !LOCAL VARIABLES:
    logical :: l_force_write   ! local version of force_write

    character(len=*), parameter :: subname = 'write_history'
    !-----------------------------------------------------------------------

    l_force_write = .false.
    if (present(force_write)) then
       l_force_write = force_write
    end if

    if (l_force_write .or. this%is_time_to_write_hist(EClock)) then
       call glc_io_write_history(instance, EClock, &
            this%history_vars, this%history_frequency_string())
    end if
    
  end subroutine write_history

  !-----------------------------------------------------------------------
  subroutine set_history_vars(this, history_vars)
    !
    ! !DESCRIPTION:
    ! Set the list of history variables
    !
    ! !USES:
    use glc_exit_mod, only : exit_glc, sigAbort
    use glc_constants, only : stdout
    !
    ! !ARGUMENTS:
    class(history_tape_base_type), intent(inout) :: this
    character(len=*), intent(in) :: history_vars
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'set_history_vars'
    !-----------------------------------------------------------------------

    if (len_trim(history_vars) > len(this%history_vars)) then
       write(stdout,*) subname//' ERROR: too-long history vars: <', trim(history_vars), '>'
       call exit_glc(sigAbort, subname//' ERROR: too-long history vars')
    end if
       
    this%history_vars = trim(history_vars)
    
  end subroutine set_history_vars

end module history_tape_base
