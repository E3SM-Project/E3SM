!$Id: calendar.F90 5867 2012-07-03 21:06:44Z dschanen@uwm.edu $
module calendar

  implicit none

  public :: gregorian2julian_date, julian2gregorian_date,  & 
            leap_year, compute_current_date, & 
            gregorian2julian_day

  private ! Default Scope

  ! Constant Parameters

  ! 3 Letter Month Abbreviations
  character(len=3), dimension(12), public, parameter :: & 
    month_names = (/'JAN','FEB','MAR','APR','MAY','JUN', & 
                    'JUL','AUG','SEP','OCT','NOV','DEC'/)

  ! Number of days per month (Jan..Dec) for a non leap year
  integer, public, dimension(12), parameter :: & 
    days_per_month = (/31, 28, 31, 30, 31, 30, & 
                       31, 31, 30, 31, 30, 31/)

  contains
!-----------------------------------------------------------------------
  integer function gregorian2julian_date( day, month, year )
!
! Description:
!   Computes the Julian Date (gregorian2julian), or the number of days since
!   1 January 4713 BC, given a Gregorian Calender date (day, month, year).
!
! Reference:
!   Fliegel, H. F. and van Flandern, T. C.,
!   Communications of the ACM, Vol. 11, No. 10 (October, 1968)
!----------------------------------------------------------------------

    implicit none

    ! Input Variables
    integer, intent(in) ::  & 
      day,        & ! Gregorian Calendar Day for given Month        [dd]
      month,      & ! Gregorian Calendar Month for given Year       [mm]
      year          ! Gregorian Calendar Year                       [yyyy]

    ! Local Variables
    integer :: I,J,K

    I = year
    J = month
    K = day

    gregorian2julian_date = K-32075+1461*(I+4800+(J-14)/12)/4+367* & 
           (J-2-(J-14)/12*12)/12-3*((I+4900+(J-14)/12)/100)/4

    return
  end function gregorian2julian_date

!------------------------------------------------------------------
  subroutine julian2gregorian_date & 
               ( julian_date, day, month, year )
!
! Description:
!   Computes the Gregorina Calendar date (day, month, year)
!   given the Julian date (julian_date).
!
! Reference:
!   Fliegel, H. F. and van Flandern, T. C.,
!   Communications of the ACM, Vol. 11, No. 10 (October, 1968)
!   http://portal.acm.org/citation.cfm?id=364097
!------------------------------------------------------------------
    implicit none

    ! Input Variable(s)
    integer, intent(in) :: julian_date ! Julian date being converted from

    ! Output Variable(s)
    integer, intent(out)::  & 
      day,     & ! Gregorian calender day for given Month       [dd]
      month,   & ! Gregorian calender month for given Year      [mm]
      year       ! Gregorian calender year                      [yyyy]

    ! Local Variables
    integer :: i, j, k, n, l

    ! ---- Begin Code ----

    L = julian_date+68569 ! Known magic number
    N = 4*L/146097 ! Known magic number
    L = L-(146097*N+3)/4 ! Known magic number
    I = 4000*(L+1)/1461001 ! Known magic number
    L = L-1461*I/4+31 ! Known magic number
    J = 80*L/2447 ! Known magic number
    K = L-2447*J/80 ! Known magic number
    L = J/11 ! Known magic number
    J = J+2-12*L ! Known magic number
    I = 100*(N-49)+I+L ! Known magic number

    year = I
    month = J
    day = K

    return

  end subroutine julian2gregorian_date

!-----------------------------------------------------------------------------
  logical function leap_year( year )
!
! Description:
!   Determines if the given year is a leap year.
!
! References:
!   None
!-----------------------------------------------------------------------------
    implicit none

    ! External
    intrinsic :: mod

    ! Input Variable(s)
    integer, intent(in) :: year ! Gregorian Calendar Year [yyyy]

    ! ---- Begin Code ----

    leap_year = ( (mod( year, 4 ) == 0) .and. & 
         (.not.(  mod( year, 100 ) == 0 .and. mod( year, 400 ) /= 0 ) ) )

    return
  end function leap_year

!----------------------------------------------------------------------------
  subroutine compute_current_date( previous_day, previous_month, & 
                                   previous_year,  & 
                                   seconds_since_previous_date, & 
                                   current_day, current_month, & 
                                   current_year, & 
                                   seconds_since_current_date )
!
! Description: 
!   Computes the current Gregorian date from a previous date and
!   the seconds that have transpired since that date.
!
! References:
!   None
!----------------------------------------------------------------------------
    use clubb_precision, only: & 
      time_precision  ! Variable(s)

    use constants_clubb, only: & 
      sec_per_day     ! Variable(s)

    implicit none

    ! Input Variable(s)

    ! Previous date
    integer, intent(in) :: & 
      previous_day,    & ! Day of the month      [dd]
      previous_month,  & ! Month of the year     [mm]
      previous_year      ! Year                  [yyyy]

    real(kind=time_precision), intent(in) :: & 
      seconds_since_previous_date ! [s]

    ! Output Variable(s)

    ! Current date
    integer, intent(out) :: & 
      current_day,     & ! Day of the month      [dd]
      current_month,   & ! Month of the year     [mm]
      current_year       ! Year                  [yyyy]

    real(kind=time_precision), intent(out) :: & 
      seconds_since_current_date

    integer :: & 
      days_since_1jan4713bc, & 
      days_since_start

    ! ---- Begin Code ----

    ! Using Julian dates we are able to add the days that the model
    ! has been running

    ! Determine the Julian Date of the starting date,
    !    written in Gregorian (day, month, year) form
    days_since_1jan4713bc = gregorian2julian_date( previous_day,  & 
                                     previous_month, previous_year )

    ! Determine the amount of days that have passed since start date
    days_since_start =  & 
          floor( seconds_since_previous_date / sec_per_day )

    ! Set days_since_1jan4713 to the present Julian date
    days_since_1jan4713bc = days_since_1jan4713bc + days_since_start

    ! Set Present time to be seconds since the Julian date
    seconds_since_current_date = seconds_since_previous_date &
      - ( real( days_since_start, kind=time_precision ) * sec_per_day )

    call julian2gregorian_date & 
           ( days_since_1jan4713bc, & 
             current_day, current_month, current_year )

    return
  end subroutine compute_current_date

!-------------------------------------------------------------------------------------
  integer function gregorian2julian_day( day, month, year )
!
! Description: 
!   This subroutine determines the Julian day (1-366)
!   for a given Gregorian calendar date(e.g. July 1, 2008).
!
! References:
!   None
!-------------------------------------------------------------------------------------

    implicit none

    ! External
    intrinsic :: sum

    ! Input Variable(s)
    integer, intent(in) :: & 
     day,             & ! Day of the Month      [dd]
     month,           & ! Month of the Year     [mm]
     year               ! Year                  [yyyy]

    ! ---- Begin Code ----

    ! Add the days from the previous months
    gregorian2julian_day = day + sum( days_per_month(1:month-1) )

    ! Kluge for a leap year
    ! If the date were 29 Feb 2000 this would not increment julian_day
    ! However 01 March 2000 would need the 1 day bump
    if ( leap_year( year ) .and. month > 2 ) then
      gregorian2julian_day = gregorian2julian_day + 1
    end if

    if ( ( leap_year( year ) .and. gregorian2julian_day > 366 ) .or. & 
         ( .not. leap_year( year ) .and. gregorian2julian_day > 365 ) ) then
      stop "Problem with Julian day conversion in gregorian2julian_day."
    end if

    return
  end function gregorian2julian_day

end module calendar
