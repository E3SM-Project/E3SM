!BOP ===========================================================================
!
! !MODULE: shr_log_mod -- variables and methods for logging
!
! !DESCRIPTION:
!    Low-level shared variables for logging.
!
!    Also, routines for generating log file messages.
!
! !INTERFACE: ------------------------------------------------------------------

module shr_log_mod

! !USES:

  use shr_kind_mod
  use shr_strconvert_mod, only: toString

  use, intrinsic :: iso_fortran_env, only: output_unit

  implicit none
  private

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

  public :: shr_log_errMsg
  public :: shr_log_OOBMsg

! !PUBLIC DATA MEMBERS:

  public :: shr_log_Level
  public :: shr_log_Unit

!EOP

  ! low-level shared variables for logging, these may not be parameters
  integer(SHR_KIND_IN) :: shr_log_Level = 1
  integer(SHR_KIND_IN) :: shr_log_Unit  = output_unit

contains

!===============================================================================
!BOP ===========================================================================
!
! !IROUTINE: shr_log_errMsg -- Return an error message containing file & line info
!
! !DESCRIPTION:
!     Return an error message containing file & line info
!     \newline
!     errMsg = shr\_log\_errMsg(__FILE__, __LINE__)
!
! This is meant to be used when a routine expects a string argument for some message,
! but you want to provide file and line information.
!
! However: Note that the performance of this function can be very bad. It is currently
! maintained because it is used by old code, but you should probably avoid using this
! in new code if possible.
!
! !REVISION HISTORY:
!     2013-July-23 - Bill Sacks
!
! !INTERFACE: ------------------------------------------------------------------

pure function shr_log_errMsg(file, line)

! !INPUT/OUTPUT PARAMETERS:

  character(len=SHR_KIND_CX)   :: shr_log_errMsg
  character(len=*), intent(in) :: file
  integer         , intent(in) :: line

!EOP

  shr_log_errMsg = 'ERROR in '//trim(file)//' at line '//toString(line)

end function shr_log_errMsg

! Create a message for an out of bounds error.
pure function shr_log_OOBMsg(operation, bounds, idx) result(OOBMsg)

  ! A name for the operation being attempted when the bounds error
  ! occurred. A string containing the subroutine name is ideal, but more
  ! generic descriptions such as "read", "modify", or "insert" could be used.
  character(len=*), intent(in) :: operation

  ! Upper and lower bounds allowed for the operation.
  integer, intent(in) :: bounds(2)

  ! Index at which access was attempted.
  integer, intent(in) :: idx

  ! Output message
  character(len=:), allocatable :: OOBMsg

  allocate(OOBMsg, source=(operation//": "//toString(idx)//" not in range ["//&
       toString(bounds(1))//", "//toString(bounds(2))//"]."))

end function shr_log_OOBMsg

end module shr_log_mod
