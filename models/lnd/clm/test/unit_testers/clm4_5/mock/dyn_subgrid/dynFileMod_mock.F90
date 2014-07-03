module dynFileMod

  ! This is a mock replacement for dynFileMod. It bypasses all of the netcdf-related
  ! stuff, instead allowing direct specification of the possible set of years and the
  ! current year. Thus, it is essentially just a wrapper to a dyn_time_info variable.

  use dynTimeInfoMod, only : time_info_type
  use ncdio_pio, only : file_desc_t
  implicit none
  save
  private

  public :: dyn_file_type

  ! Note that this is intended to be used with the mock form of file_desc_t, defined in
  ! ncdio_pio_mock.F90
  type, extends(file_desc_t) :: dyn_file_type
     private
     type(time_info_type) :: time_info
   contains
     procedure :: update_time_info     ! should be called every time step to update time information

     ! The following are pass-through methods to time_info:
     procedure :: get_nt1              ! get lower bound index of current interval
     procedure :: get_nt2              ! get upper bound index of current interval
     procedure :: get_year             ! get the year associated with a given time index
     procedure :: is_within_bounds     ! return true if we are currently within the bounds of this file
  end type dyn_file_type

  interface dyn_file_type
     module procedure constructor  ! initialize a new dyn_file_type object
  end interface dyn_file_type

contains
  
  ! ======================================================================
  ! Constructors
  ! ======================================================================

  type(dyn_file_type) function constructor(my_years, cur_year)
    ! Note that this should be used with the mock form of file_desc_t, defined in
    ! ncdio_pio_mock.F90

    integer, intent(in) :: my_years(:)  ! all years desired for the time_info variable
    integer, intent(in) :: cur_year     ! current model year

    ! The following only works if we're using the mock form of file_desc_t, defined in
    ! ncdio_pio_mock.F90
    constructor%file_desc_t = file_desc_t()

    constructor%time_info = time_info_type(my_years, cur_year)
  end function constructor

  ! ======================================================================
  ! Public methods
  ! ======================================================================

  subroutine update_time_info(this, cur_year)
    class(dyn_file_type), intent(inout) :: this
    integer, intent(in) :: cur_year  ! current model year

    call this%time_info%update_time_info(cur_year)
  end subroutine update_time_info

  pure integer function get_nt1(this)
    class(dyn_file_type), intent(in) :: this

    get_nt1 = this%time_info%get_nt1()
  end function get_nt1

  pure integer function get_nt2(this)
    class(dyn_file_type), intent(in) :: this

    get_nt2 = this%time_info%get_nt2()
  end function get_nt2

  pure integer function get_year(this, nt)
    class(dyn_file_type), intent(in) :: this
    integer, intent(in) :: nt  ! time index

    get_year = this%time_info%get_year(nt)
  end function get_year

  pure logical function is_within_bounds(this)
    class(dyn_file_type), intent(in) :: this

    is_within_bounds = this%time_info%is_within_bounds()
  end function is_within_bounds

end module dynFileMod
