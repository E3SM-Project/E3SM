module dynTimeInfoMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Contains a derived type and associated methods for storing and working with time
  ! information for a single dynamic landuse file. The assumption is that there is a
  ! single time sample per year.
  !
  ! !USES:
  use clm_varctl     , only : iulog
  use shr_log_mod    , only : errMsg => shr_log_errMsg
  use abortutils     , only : endrun

  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: time_info_type

  type time_info_type
     private
     ! Static information about the file:
     integer :: nyears                ! number of years in the file
     integer, allocatable :: years(:) ! all years in this file
     
     ! Information that potentially changes each time step:
     integer :: time_index_lower           ! lower bound index of the current interval
     integer :: time_index_upper           ! upper bound index of the current interval

   contains
     generic :: set_current_year => &   ! should be called every time step to update info with the current model year
          set_current_year_get_year, &
          set_current_year_from_year
     procedure :: set_current_year_get_year  ! version of set_current_year that obtains current model year
     procedure :: set_current_year_from_year ! version of set_current_year that assumes you already have obtained the current model year
     procedure :: get_time_index_lower               ! get lower bound index of current interval
     procedure :: get_time_index_upper               ! get upper bound index of current interval
     procedure :: get_year              ! get the year associated with a given time index
     procedure :: is_within_bounds      ! return true if we are currently within the bounds of this file
     procedure :: is_before_time_series ! returns true if we are currently prior to the bounds of this file
     procedure :: is_after_time_series  ! returns true if we are currently after the bounds of this file (if the last year of the file is (e.g.) 2005, then this is TRUE if the current year is 2005)

     procedure, private :: year_in_current_interval ! returns true if the current year is in the current interval
  end type time_info_type

  interface time_info_type
     module procedure constructor                 ! initialize a time_info_type object
  end interface time_info_type

contains

  ! ======================================================================
  ! Constructors
  ! ======================================================================

  !-----------------------------------------------------------------------
  type(time_info_type) function constructor(my_years, cur_year)
    !
    ! !DESCRIPTION:
    ! Initialize a time_info_type object
    !
    ! !ARGUMENTS:
    integer, intent(in) :: my_years(:) ! all years in this file
    integer, intent(in) :: cur_year    ! current model year
    !-----------------------------------------------------------------------

    constructor%nyears = size(my_years)
    allocate(constructor%years(constructor%nyears))
    constructor%years = my_years

    ! Set time_index_lower and time_index_upper arbitrarily; they'll get set correctly by set_current_year
    constructor%time_index_lower = 1
    constructor%time_index_upper = 1

    ! Set time_index_lower and time_index_upper to their correct values
    call constructor%set_current_year(cur_year)

  end function constructor


  ! ======================================================================
  ! Public methods
  ! ======================================================================

  !-----------------------------------------------------------------------
  subroutine set_current_year_get_year(this)
    !
    ! !DESCRIPTION:
    ! Update time information (time_index_lower and time_index_upper)
    !
    ! Should be called every time step
    !
    ! This version of set_current_year obtains the current model year for you.
    !
    ! !USES:
    use clm_time_manager , only : get_curr_date
    !
    ! !ARGUMENTS:
    class(time_info_type), intent(inout) :: this      ! this object
    !
    ! !LOCAL VARIABLES:
    integer  :: year             ! year (0, ...) for nstep+1
    integer  :: mon              ! month (1, ..., 12) for nstep+1
    integer  :: day              ! day of month (1, ..., 31) for nstep+1
    integer  :: sec              ! seconds into current date for nstep+1
    
    character(len=*), parameter :: subname = 'set_current_year_get_year'
    !-----------------------------------------------------------------------
    
    call get_curr_date(year, mon, day, sec)
    call this%set_current_year(year)

  end subroutine set_current_year_get_year


  !-----------------------------------------------------------------------
  subroutine set_current_year_from_year(this, cur_year)
    !
    ! !DESCRIPTION: 
    ! Update time information (time_index_lower and time_index_upper)
    !
    ! Should be called every time step
    !
    ! This version of set_current_year assumes you have already obtained the current
    ! model year (passed via the cur_year argument).
    !
    ! !ARGUMENTS:
    class(time_info_type), intent(inout) :: this      ! this object
    integer, intent(in)                  :: cur_year  ! current model year
    !
    ! !LOCAL VARIABLES:
    logical :: found   ! has the correct interval been found?
    integer :: n       ! interval index

    character(len=*), parameter :: subname = 'set_current_year_from_year'
    !-----------------------------------------------------------------------

    ! Determine if current date spans the years
    !
    ! If current year is less than first timeseries year, then use the first year from
    ! dynamic land cover file for both time_index_lower and time_index_upper, forcing constant weights until the
    ! model year enters the dynamic land cover dataset timeseries range.
    !
    ! If current year is equal to or greater than the last timeseries year, then use the
    ! last year for both time_index_lower and time_index_upper, forcing constant weights for the remainder of the
    ! simulation.
    !
    ! This mechanism permits the introduction of a dynamic pft period in the middle of a
    ! simulation, with constant weights before and after the dynamic period.

    associate( &
         nyears => this%nyears, &
         years  => this%years, &
         time_index_lower    => this%time_index_lower, &
         time_index_upper    => this%time_index_upper)

    if (year_in_current_interval(this, cur_year)) then
       ! DO NOTHING - NT1 AND NT2 ARE ALREADY CORRECT
    else
       if (cur_year < years(1)) then
          ! prior to the first interval
          time_index_lower = 1
          time_index_upper = 1
       else if (cur_year >= years(nyears)) then
          ! past the last interval
          time_index_lower = nyears
          time_index_upper = nyears
       else
          ! within the time bounds of the file
          found = .false.
          do n = 1, nyears-1
             if (cur_year == years(n)) then
                time_index_lower = n
                time_index_upper = n + 1
                found = .true.
                exit
             end if
          end do
          if (.not. found) then
             write(iulog,*) subname//' ERROR: model year not found in pftdyn timeseries'
             write(iulog,*)'model year = ',cur_year
             call endrun(msg=errMsg(__FILE__, __LINE__))
          end if
       end if
    end if

    SHR_ASSERT(time_index_upper <= nyears, subname // ': time_index_upper should not be greater than nyears')
          
    end associate

  end subroutine set_current_year_from_year


  ! ----------------------------------------------------------------------
  ! Various getter routines
  ! ----------------------------------------------------------------------

  !-----------------------------------------------------------------------
  pure integer function get_time_index_lower(this)
    ! !DESCRIPTION: Get lower bound index of current interval
    !
    ! !ARGUMENTS:
    class(time_info_type), intent(in) :: this
    !-----------------------------------------------------------------------

    get_time_index_lower = this%time_index_lower
  end function get_time_index_lower

  !-----------------------------------------------------------------------
  pure integer function get_time_index_upper(this)
    ! !DESCRIPTION: Get upper bound index of current interval
    !
    ! !ARGUMENTS:
    class(time_info_type), intent(in) :: this
    !-----------------------------------------------------------------------

    get_time_index_upper = this%time_index_upper
  end function get_time_index_upper

  !-----------------------------------------------------------------------
  integer function get_year(this, nt)
    ! !DESCRIPTION: Get the year associated with time index nt
    !
    ! Note this can't be a pure function because of the call to shr_assert
    !
    ! !ARGUMENTS:
    class(time_info_type), intent(in) :: this
    integer              , intent(in) :: nt    ! time index

    character(len=*), parameter :: subname = 'get_year'
    !-----------------------------------------------------------------------

    SHR_ASSERT(1 <= nt .and. nt <= this%nyears, subname // ': nt out of bounds')
    get_year = this%years(nt)
  end function get_year

  !-----------------------------------------------------------------------
  pure logical function is_within_bounds(this)
    ! !DESCRIPTION: Returns true if we are currently within the bounds of this file
    !
    ! !ARGUMENTS:
    class(time_info_type), intent(in) :: this
    !-----------------------------------------------------------------------
    
    is_within_bounds = ((.not. this%is_before_time_series()) .and. &
         (.not. this%is_after_time_series()))

  end function is_within_bounds

  ! ======================================================================
  ! Private methods
  ! ======================================================================

  !-----------------------------------------------------------------------
  pure logical function year_in_current_interval(this, cur_year)
    ! !DESCRIPTION:
    ! Returns true if the current year is in the current interval
    !
    ! !ARGUMENTS:
    class(time_info_type), intent(in) :: this      ! this object
    integer, intent(in)               :: cur_year  ! current model year
    !-----------------------------------------------------------------------

    if (this%years(this%time_index_lower) == cur_year .and. this%years(this%time_index_upper) == (cur_year + 1)) then
       ! Normal case: we're within the time series, in the same interval as before
       year_in_current_interval = .true.
    else if (this%is_before_time_series() .and. cur_year < this%years(1)) then
       ! We were and still are before the time series
       year_in_current_interval = .true.
    else if (this%is_after_time_series() .and. cur_year >= this%years(this%nyears)) then
       ! We were and still are after the time series
       year_in_current_interval = .true.
    else
       year_in_current_interval = .false.
    end if

  end function year_in_current_interval

  !-----------------------------------------------------------------------
  pure logical function is_before_time_series(this)
    ! !DESCRIPTION: Returns true if we are currently prior to the bounds of this file
    !
    ! !ARGUMENTS:
    class(time_info_type), intent(in) :: this
    !-----------------------------------------------------------------------

    if (this%time_index_upper == 1) then
       is_before_time_series = .true.
    else
       is_before_time_series = .false.
    end if
  end function is_before_time_series

  !-----------------------------------------------------------------------
  pure logical function is_after_time_series(this)
    ! !DESCRIPTION: Returns true if we are currently after the bounds of this file
    !
    ! If the last year of the file is (e.g.) 2005, then this is TRUE if the current year
    ! is 2005
    !
    ! !ARGUMENTS:
    class(time_info_type), intent(in) :: this
    !-----------------------------------------------------------------------

    if (this%time_index_lower == this%nyears) then
       is_after_time_series = .true.
    else
       is_after_time_series = .false.
    end if
  end function is_after_time_series

end module dynTimeInfoMod
