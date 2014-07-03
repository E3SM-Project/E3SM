module dynFileMod
  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains a derived type and associated methods for storing and working with
  ! information for a single dynamic landuse file. This largely consists of the time
  ! information, which is handled by the time_info object contained here. But there is
  ! also other file information, such as the ncid. Note that the time_info class assumes
  ! that there is a single time sample per year.
  !
  ! !USES:
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use dynTimeInfoMod , only : time_info_type
  use ncdio_pio      , only : file_desc_t, ncd_pio_openfile, ncd_inqdid, ncd_inqdlen, ncd_io
  use clm_varctl     , only : iulog
  use abortutils     , only : endrun
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: dyn_file_type

  type, extends(file_desc_t) :: dyn_file_type
     private
     character(len=256)   :: locfn     ! local file name
     type(time_info_type) :: time_info ! time information for this file
   contains
     procedure :: update_time_info     ! should be called every time step to update time information

     ! The following are pass-through methods to time_info:
     procedure :: get_nt1               ! get lower bound index of current interval
     procedure :: get_nt2               ! get upper bound index of current interval
     procedure :: get_year              ! get the year associated with a given time index
     procedure :: is_within_bounds      ! return true if we are currently within the bounds of this file
     procedure :: is_before_time_series ! returns true if we are currently prior to the bounds of this file
     procedure :: is_after_time_series  ! returns true if we are currently after the bounds of this file (if the last year of the file is (e.g.) 2005, then this is TRUE if the current year is 2005)
     
  end type dyn_file_type

  interface dyn_file_type
     module procedure constructor  ! initialize a new dyn_file_type object
  end interface dyn_file_type

contains
  
  ! ======================================================================
  ! Constructors
  ! ======================================================================

  !-----------------------------------------------------------------------
  type(dyn_file_type) function constructor(filename)
    !
    ! !DESCRIPTION:
    ! Initialize a dyn_file_type object
    !
    ! Opens the file associated with filename for reading, reads the 'YEAR' variable from
    ! this file (assumed to have dimension 'time'), and initializes a dyn_time_info object
    ! based on this YEAR variable and the current model year.
    !
    ! !USES:
    use fileutils        , only : getfil
    use clm_time_manager , only : get_curr_date
    !
    ! !ARGUMENTS:
    character(len=*), intent(in) :: filename
    !
    ! !LOCAL VARIABLES:
    integer :: ier                   ! error code
    integer :: ntimes                ! number of time samples
    integer :: varid                 ! netcdf variable ID
    integer :: cur_year              ! year (0, ...) for nstep+1
    integer :: mon                   ! month (1, ..., 12) for nstep+1
    integer :: day                   ! day of month (1, ..., 31) for nstep+1
    integer :: sec                   ! seconds into current date for nstep+1
    integer, allocatable :: years(:) ! years in the file

    character(len=*), parameter :: subname = 'dyn_file_type constructor'
    !-----------------------------------------------------------------------

    ! Obtain file

    call getfil(filename, constructor%locfn, 0)
    call ncd_pio_openfile(constructor, constructor%locfn, 0)
    
    ! Obtain years

    call ncd_inqdid(constructor, 'time', varid)
    call ncd_inqdlen(constructor, varid, ntimes)
    allocate(years(ntimes), stat=ier)
    if (ier /= 0) then
       call endrun(msg=' allocation error for years'//errMsg(__FILE__, __LINE__))
    end if
    call ncd_io(ncid=constructor, varname='YEAR', flag='read', data=years)
    
    ! Initialize object containing time information for the file

    call get_curr_date(cur_year, mon, day, sec)

    constructor%time_info = time_info_type(years, cur_year)

    deallocate(years)

  end function constructor


  ! ======================================================================
  ! Public methods
  ! ======================================================================

  !-----------------------------------------------------------------------
  subroutine update_time_info(this)
    !
    ! !DESCRIPTION:
    ! Update time information for this file, based on the current model time
    ! 
    ! Should be called every time step
    !
    ! !USES:
    use clm_time_manager , only : get_curr_date
    !
    ! !ARGUMENTS:
    class(dyn_file_type), intent(inout) :: this ! this object
    !
    ! !LOCAL VARIABLES:
    integer  :: year             ! year (0, ...) for nstep+1
    integer  :: mon              ! month (1, ..., 12) for nstep+1
    integer  :: day              ! day of month (1, ..., 31) for nstep+1
    integer  :: sec              ! seconds into current date for nstep+1
    !-----------------------------------------------------------------------
    
    call get_curr_date(year, mon, day, sec)
    call this%time_info%update_time_info(year)
  end subroutine update_time_info

  ! ----------------------------------------------------------------------
  ! The following are pass-through methods to the time_info object
  ! ----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  pure integer function get_nt1(this)
    ! !DESCRIPTION: Get lower bound index of current interval
    !
    ! !ARGUMENTS:
    class(dyn_file_type), intent(in) :: this
    !-----------------------------------------------------------------------

    get_nt1 = this%time_info%get_nt1()
  end function get_nt1

  !-----------------------------------------------------------------------
  pure integer function get_nt2(this)
    ! !DESCRIPTION: Get upper bound index of current interval
    !
    ! !ARGUMENTS:
    class(dyn_file_type), intent(in) :: this
    !-----------------------------------------------------------------------

    get_nt2 = this%time_info%get_nt2()
  end function get_nt2

  !-----------------------------------------------------------------------
  integer function get_year(this, nt)
    ! !DESCRIPTION: Get the year associated with time index nt
    !
    ! !ARGUMENTS:
    class(dyn_file_type), intent(in) :: this
    integer             , intent(in) :: nt    ! time index
    !-----------------------------------------------------------------------
    
    get_year = this%time_info%get_year(nt)
  end function get_year

  !-----------------------------------------------------------------------
  pure logical function is_within_bounds(this)
    ! !DESCRIPTION: Returns true if we are currently within the bounds of this file
    !
    ! !ARGUMENTS:
    class(dyn_file_type), intent(in) :: this
    !-----------------------------------------------------------------------
    
    is_within_bounds = this%time_info%is_within_bounds()
  end function is_within_bounds

  !-----------------------------------------------------------------------
  pure logical function is_before_time_series(this)
    ! !DESCRIPTION: Returns true if we are currently prior to the bounds of this file
    !
    ! !ARGUMENTS:
    class(dyn_file_type), intent(in) :: this
    !-----------------------------------------------------------------------
    
    is_before_time_series = this%time_info%is_before_time_series()
  end function is_before_time_series

  !-----------------------------------------------------------------------
  pure logical function is_after_time_series(this)
    ! !DESCRIPTION: Returns true if we are currently after the bounds of this file
    !
    ! If the last year of the file is (e.g.) 2005, then this is TRUE if the current year
    ! is 2005
    !
    ! !ARGUMENTS:
    class(dyn_file_type), intent(in) :: this
    !-----------------------------------------------------------------------
    
    is_after_time_series = this%time_info%is_after_time_series()
  end function is_after_time_series

end module dynFileMod
