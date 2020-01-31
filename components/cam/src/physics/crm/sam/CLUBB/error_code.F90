!-------------------------------------------------------------------------------
! $Id: error_code.F90 5906 2012-08-10 23:20:05Z dschanen@uwm.edu $
!-------------------------------------------------------------------------------

module error_code

! Description:
!   Since f90/95 lacks enumeration, we're stuck numbering each
!   error code by hand like this.

!   We are "enumerating" error codes to be used with CLUBB. Adding
!   additional codes is as simple adding an additional integer
!   parameter. The error codes are ranked by severity, the higher
!   number being more servere. When two errors occur, assign the
!   most servere to the output.

!   This code also handles subroutines related to debug_level. See
!   the 'set_clubb_debug_level' description for more detail.

! References:
!   None
!-------------------------------------------------------------------------------

  implicit none

  private ! Default Scope

  public :: & 
    reportError,  & 
    fatal_error, & 
    lapack_error,     & 
    clubb_at_least_debug_level,  & 
    set_clubb_debug_level, & 
    clubb_debug

  private :: clubb_debug_level

  ! Model-Wide Debug Level
  integer, save :: clubb_debug_level = 0

!$omp threadprivate(clubb_debug_level)

  ! Error Code Values
  integer, parameter, public :: & 
    clubb_no_error                 =  0, & 
    clubb_var_less_than_zero       =  1, & 
    clubb_var_equals_NaN           =  2, & 
    clubb_singular_matrix          =  3, & 
    clubb_bad_lapack_arg           =  4, & 
    clubb_rtm_level_not_found      =  5, & 
    clubb_var_out_of_bounds        =  6, &
    clubb_var_out_of_range         =  7

  contains

!-------------------------------------------------------------------------------
  subroutine reportError( err_code )
!
! Description: 
!   Reports meaning of error code to console.
!
!-------------------------------------------------------------------------------

    use constants_clubb, only: & 
        fstderr ! Variable(s)

    implicit none

    ! Input Variable
    integer, intent(in) :: err_code ! Error Code being examined

    ! ---- Begin Code ----

    select case ( err_code )

    case ( clubb_no_error )
      write(fstderr,*) "No errors reported."

    case ( clubb_var_less_than_zero )
      write(fstderr,*) "Variable in CLUBB is less than zero."

    case ( clubb_singular_matrix )
      write(fstderr,*) "Singular Matrix in CLUBB."

    case ( clubb_var_equals_NaN )
      write(fstderr,*) "Variable in CLUBB is NaN."

    case ( clubb_bad_lapack_arg )
      write(fstderr,*) "Argument passed to a LAPACK procedure is invalid."

    case ( clubb_rtm_level_not_found )
      write(fstderr,*) "rtm level not found"

    case ( clubb_var_out_of_bounds )
      write(fstderr,*) "Input variable is out of bounds."

    case ( clubb_var_out_of_range )
      write(fstderr,*) "A CLUBB variable had a value outside the valid range."

    case default
      write(fstderr,*) "Unknown error: ", err_code

    end select

    return
  end subroutine reportError
!-------------------------------------------------------------------------------
  elemental function lapack_error( err_code )
!
! Description: 
!   Checks to see if the err_code is equal to one
!   caused by an error encountered using LAPACK.
! Reference:
!   None
!-------------------------------------------------------------------------------
    implicit none

    ! Input variable
    integer,intent(in) :: err_code ! Error Code being examined

    ! Output variable
    logical :: lapack_error

    ! ---- Begin Code ----

    lapack_error = (err_code == clubb_singular_matrix .or. & 
        err_code == clubb_bad_lapack_arg )

    return
  end function lapack_error

!-------------------------------------------------------------------------------
  elemental function fatal_error( err_code )
!
! Description: Checks to see if the err_code is one that usually
!   causes an exit in other parts of CLUBB.
! References:
!   None
!-------------------------------------------------------------------------------
    implicit none

    ! Input Variable
    integer, intent(in) :: err_code ! Error Code being examined

    ! Output variable
    logical :: fatal_error

    ! ---- Begin Code ----

    fatal_error = err_code /= clubb_no_error .and. & 
                  err_code /= clubb_var_less_than_zero
    return
  end function fatal_error

!------------------------------------------------------------------	
  logical function clubb_at_least_debug_level( level )
!
! Description:
!   Checks to see if clubb has been set to a specified debug level
!------------------------------------------------------------------
    implicit none

    ! Input variable
    integer, intent(in) :: level   ! The debug level being checked against the current setting

    ! ---- Begin Code ----

    clubb_at_least_debug_level = ( level <= clubb_debug_level )

    return
  end function clubb_at_least_debug_level

!-------------------------------------------------------------------------------
  subroutine set_clubb_debug_level( level )
!
!  Description:
!    Accessor for clubb_debug_level
!
!   0 => Print no debug messages to the screen
!   1 => Print lightweight debug messages, e.g. print statements
!   2 => Print debug messages that require extra testing,
!        e.g. checks for NaNs and spurious negative values.
!  References:
!    None
!-------------------------------------------------------------------------------
    implicit none

    ! Input variable
    integer, intent(in) :: level ! The debug level being checked against the current setting

    ! ---- Begin Code ----

    clubb_debug_level = level

    return
  end subroutine set_clubb_debug_level

!-------------------------------------------------------------------------------
  subroutine clubb_debug( level, str )
!
! Description:
!   Prints a message to file unit fstderr if the level is greater
!   than or equal to the current debug level.
!-------------------------------------------------------------------------------
    use constants_clubb, only: & 
        fstderr ! Variable(s)

    implicit none

    ! Input Variable(s)

    character(len=*), intent(in) :: str ! The message being reported

    ! The debug level being checked against the current setting
    integer, intent(in) :: level

    ! ---- Begin Code ----

    if ( level <= clubb_debug_level ) then
      write(fstderr,*) str
    end if

    return
  end subroutine clubb_debug

end module error_code
!-------------------------------------------------------------------------------
