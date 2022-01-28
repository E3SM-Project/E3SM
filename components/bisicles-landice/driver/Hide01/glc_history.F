!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

module glc_history

  !BOP
  ! !MODULE: glc_history

  ! !DESCRIPTION:
  ! Contains routines for handling history output.
  !
  ! Usage:
  !
  ! - In initialization, call glc_history_init
  !
  ! - Every time through the run loop, call glc_history_write
  !
  ! !USES:
  use glc_kinds_mod
  use history_tape_base , only : history_tape_base_type, len_history_vars
  use shr_kind_mod      , only : CL=>SHR_KIND_CL, CXX=>SHR_KIND_CXX
  use glc_exit_mod      , only : exit_glc, sigAbort
  use glc_constants     , only : nml_in, stdout, blank_fmt, ndelim_fmt
  
  implicit none
  private
  save

  ! !PUBLIC ROUTINES:
  public :: glc_history_init  ! initialize the history_tape instance
  public :: glc_history_write ! write to history file, if it's time to do so
  
  ! !PRIVATE ROUTINES:
  private :: read_namelist
  
  ! !PRIVATE MODULE VARIABLES:

  ! TODO(wjs, 2015-02-18) Eventually, we may want to allow for multiple history tapes. In
  ! that case, we should replace this scalar variable with an array. We would also need
  ! to modify the code in this module to read namelist options for all history tapes, and
  ! then have a loop that creates all history tape objects. Note that the history tape
  ! index should become a field in the history tape class; this is needed to create
  ! unique time flags for each history tape (and possibly other things).
  class(history_tape_base_type), allocatable :: history_tape

  ! max character lengths
  integer, parameter :: len_history_option = CL
contains

  !------------------------------------------------------------------------
  ! PUBLIC ROUTINES
  !------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine glc_history_init
    !
    ! !DESCRIPTION:
    ! Initialize the history_tape instance
    !
    ! Should be called once, in model initialization
    !
    ! !USES:
    use glc_time_management, only : freq_opt_nyear
    use history_tape_standard, only : history_tape_standard_type
    use history_tape_coupler, only : history_tape_coupler_type
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    character(len=len_history_vars) :: cesm_history_vars
    character(len=len_history_option) :: history_option
    integer(int_kind) :: history_frequency
    
    character(len=*), parameter :: subname = 'glc_history_init'
    !-----------------------------------------------------------------------

    call read_namelist(cesm_history_vars, history_option, history_frequency)
    
    select case (history_option)
    case ('nyears')
       allocate(history_tape, source = history_tape_standard_type( &
            history_vars = cesm_history_vars, freq_opt = freq_opt_nyear, &
            freq = history_frequency))
    case ('coupler')
       allocate(history_tape, source = history_tape_coupler_type( &
            history_vars = cesm_history_vars))
    case default
       write(stdout,*) subname//' ERROR: Unhandled history_option: ', trim(history_option)
       call exit_glc(sigAbort, subname//' ERROR: Unhandled history_option')
    end select
       
  end subroutine glc_history_init

  !-----------------------------------------------------------------------
  subroutine glc_history_write(instance, EClock, force_write)
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
    use glad_type, only : glad_instance
    use esmf, only: ESMF_Clock
    !
    ! !ARGUMENTS:
    type(glad_instance), intent(inout) :: instance
    type(ESMF_Clock),     intent(in)    :: EClock
    logical, intent(in), optional :: force_write
    !-----------------------------------------------------------------------

    call history_tape%write_history(instance, EClock, force_write)
    
  end subroutine glc_history_write


  !------------------------------------------------------------------------
  ! PRIVATE ROUTINES
  !------------------------------------------------------------------------

  !-----------------------------------------------------------------------
  subroutine read_namelist(cesm_history_vars, history_option, history_frequency)
    !
    ! !DESCRIPTION:
    ! Reads the namelist containing history options
    !
    ! !USES:
    use glc_communicate , only: my_task, master_task
    use glc_files       , only: nml_filename
    use glc_broadcast   , only: broadcast_scalar
    !
    ! !ARGUMENTS:
    character(len=len_history_vars), intent(out) :: cesm_history_vars
    character(len=len_history_option), intent(out) :: history_option
    integer(int_kind), intent(out) :: history_frequency
    !
    ! !LOCAL VARIABLES:
    
    integer :: nml_error
    
    character(len=*), parameter :: subname = 'read_namelist'
    !-----------------------------------------------------------------------

    namelist /cism_history/ cesm_history_vars, history_option, history_frequency

    ! Set default values
    cesm_history_vars = ' '
    history_option = ' '
    history_frequency = 1
    
    if (my_task == master_task) then
       open(nml_in, file=nml_filename, status='old', iostat=nml_error)
       if (nml_error /= 0) then
          nml_error = -1
       else
          nml_error =  1
       end if
       do while (nml_error > 0)
          read(nml_in, nml=cism_history, iostat=nml_error)
       end do
       if (nml_error == 0) then
          close(nml_in)
       end if
    end if

    call broadcast_scalar(nml_error, master_task)
    if (nml_error /= 0) then
       call exit_glc(sigAbort,'ERROR reading cism_history namelist')
    end if

    ! Write namelist settings
    if (my_task == master_task) then
       write(stdout,blank_fmt)
       write(stdout,ndelim_fmt)
       write(stdout,blank_fmt)
       write(stdout,*) ' cism_history namelist settings:'
       write(stdout,blank_fmt)
       write(stdout, cism_history)
    end if

    ! Send namelist settings to all procs
    call broadcast_scalar(cesm_history_vars, master_task)
    call broadcast_scalar(history_option, master_task)
    call broadcast_scalar(history_frequency, master_task)

    if ((len_trim(cesm_history_vars)+3) >= len(cesm_history_vars)) then
       ! Assume that if we get within 3 spaces of the variable legth (excluding spaces)
       ! then we may be truncating the intended value
       call exit_glc(sigAbort, subname// &
            ' ERROR: The value of cesm_history_vars is too long for the variable')
    end if
    
  end subroutine read_namelist

  
end module glc_history
