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
     type(time_info_type) :: time_info
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

end module dynFileMod
