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

  implicit none
  private

! !PUBLIC TYPES:

   ! no public types

! !PUBLIC MEMBER FUNCTIONS:

  public :: shr_log_errMsg

! !PUBLIC DATA MEMBERS:

  public :: shr_log_Level
  public :: shr_log_Unit

!EOP

  ! low-level shared variables for logging, these may not be parameters
  integer(SHR_KIND_IN) :: shr_log_Level = 1
  integer(SHR_KIND_IN) :: shr_log_Unit  = 6

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
! !REVISION HISTORY:
!     2013-July-23 - Bill Sacks
!
! !INTERFACE: ------------------------------------------------------------------

function shr_log_errMsg(file, line)
  
  implicit none

! !INPUT/OUTPUT PARAMETERS:

  character(len=SHR_KIND_CX)   :: shr_log_errMsg
  character(len=*), intent(in) :: file
  integer         , intent(in) :: line

!EOP

  write(shr_log_errMsg, '(a, a, a, i0)') 'ERROR in ', trim(file), ' at line ', line

end function shr_log_errMsg

end module shr_log_mod
