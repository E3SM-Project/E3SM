module history_tape_standard
  
  ! Defines a class for a standard history tape. This class writes history files at a
  ! specified frequency (e.g., every N years).

  use history_tape_base, only : history_tape_base_type
  use glc_kinds_mod
  use glc_time_management, only : init_time_flag, check_time_flag
  
  implicit none
  private
  save

  public :: history_tape_standard_type
  type, extends(history_tape_base_type) :: history_tape_standard_type
     private
     integer(int_kind) :: freq_opt  ! frequency option (as defined in glc_time_management)
     integer(int_kind) :: freq      ! frequency (e.g., number of years for freq_opt=freq_opt_nyear)
     integer(int_kind) :: time_flag ! reference to a time flag in glc_time_management
   contains
     ! Logical function saying whether it's time to write a history file
     procedure :: is_time_to_write_hist

     ! Function returning a string describing the history frequency
     procedure :: history_frequency_string
  end type history_tape_standard_type

  interface history_tape_standard_type
     module procedure constructor
  end interface history_tape_standard_type

contains

  !-----------------------------------------------------------------------
  function constructor(history_vars, freq_opt, freq)
    !
    ! !DESCRIPTION:
    ! Creates a history_tape_standard_type object
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(history_tape_standard_type) :: constructor  ! function result

    ! List of variables to write to file
    character(len=*), intent(in) :: history_vars
    
    ! Frequency option; should be one of the options defined in glc_time_management
    ! (e.g., freq_opt_nyear)
    integer(int_kind), intent(in) :: freq_opt

    ! Frequency; e.g., for freq_opt = freq_opt_nyear, history files will be written every
    ! freq years
    integer(int_kind), intent(in) :: freq

    !-----------------------------------------------------------------------

    call constructor%set_history_vars(history_vars)
    constructor%freq_opt = freq_opt
    constructor%freq = freq

    ! TODO(wjs, 2015-02-18) If we allow multiple history tapes, we should construct a time
    ! flag name that includes the history tape index, so that we have unique time flags
    ! for each history tape. In that case, the history tape index should be passed into
    ! the constructor, and stored as a component of the class.
    constructor%time_flag = init_time_flag('do_hist', freq_opt = freq_opt, freq = freq)
    
  end function constructor

  !-----------------------------------------------------------------------
  logical function is_time_to_write_hist(this, EClock)
    !
    ! !DESCRIPTION:
    ! Returns true if it is time to write the history tape associated with this
    ! controller.
    !
    ! Note that EClock is unused in this implementation; it is simply included for
    ! compatibility with the generic interface.
    !
    ! !USES:
    use esmf, only: ESMF_Clock
    !
    ! !ARGUMENTS:
    class(history_tape_standard_type), intent(in) :: this
    type(ESMF_Clock), intent(in) :: EClock

    !-----------------------------------------------------------------------

    is_time_to_write_hist = check_time_flag(this%time_flag)
    
  end function is_time_to_write_hist

  !-----------------------------------------------------------------------
  function history_frequency_string(this)
    !
    ! !DESCRIPTION:
    ! Returns a string representation of this history frequency
    !
    ! This string representation is based on the convention for the time_period_freq
    ! metadata defined here:
    ! http://www.cesm.ucar.edu/models/cesm2.0/filename_conventions_cesm.html
    !
    ! NOTE(wjs, 2015-02-17) Design note: Really, this functionality should be delegated to
    ! a method on a time frequency class; this class would contain the various time
    ! frequency methods defined in glc_time_management, as well as this to_string
    ! method. However, I don't want to go to the effort of making a new class for that
    ! right now, so for now I'm putting this behavior here.
    !
    ! !USES:
    use glc_exit_mod, only : exit_glc, sigAbort
    use glc_constants, only : stdout
    use glc_time_management, only : freq_opt_nyear, freq_opt_nmonth, freq_opt_nday, &
         freq_opt_nhour, freq_opt_nsecond
    !
    ! !ARGUMENTS:
    character(len=:), allocatable :: history_frequency_string  ! function result
    class(history_tape_standard_type), intent(in) :: this
    !
    ! !LOCAL VARIABLES:
    character(len=:), allocatable :: freq_opt_string
    integer, parameter :: max_digits = 20

    character(len=*), parameter :: subname = 'history_frequency_string'
    !-----------------------------------------------------------------------

    select case (this%freq_opt)
    case (freq_opt_nyear)
       freq_opt_string = 'year_'
    case (freq_opt_nmonth)
       freq_opt_string = 'month_'
    case (freq_opt_nday)
       freq_opt_string = 'day_'
    case (freq_opt_nhour)
       freq_opt_string = 'hour_'
    case (freq_opt_nsecond)
       freq_opt_string = 'second_'
    case default
       write(stdout,*) subname//' ERROR: Unhandled freq_opt: ', this%freq_opt
       call exit_glc(sigAbort, subname//' ERROR: Unhandled freq_opt')
    end select

    allocate(character(len = len(freq_opt_string) + max_digits) :: &
         history_frequency_string)

    write(history_frequency_string, '(a, i0)') freq_opt_string, this%freq
    
  end function history_frequency_string

end module history_tape_standard
