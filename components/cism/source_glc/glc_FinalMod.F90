!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_FinalMod

!BOP
! !MODULE: glc_FinalMod
! !DESCRIPTION:
!  This module contains the glc finalization method that shuts down glc
!  gracefully (we hope).  It exits the message environment and checks 
!  for successful execution.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id: POP_FinalMod.F90 808 2006-04-28 17:06:38Z njn01 $
!  Adapted by William Lipscomb from POP_FinalMod.F90
!
! !USES:

   use glc_kinds_mod
   use glc_ErrorMod
   use glc_communicate, only: exit_message_environment
   use glc_global_fields, only: ice_sheet
   use glint_main, only: end_glint
   use glc_constants, only: stdout

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: glc_final

!EOP
!BOC
!-----------------------------------------------------------------------
!
!     module variables
!
!-----------------------------------------------------------------------

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: glc_final
! !INTERFACE:

 subroutine glc_final(ErrorCode)

! !DESCRIPTION:
!  This routine shuts down glc by exiting all relevent environments.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (i4), intent(inout) :: &
      ErrorCode              ! On input, error code from Init,Run method
                             ! On output, status of this routine

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
!
!  exit glint gracefully
!
!-----------------------------------------------------------------------

   call end_glint(ice_sheet, close_logfile=.false.)

!-----------------------------------------------------------------------
!
!  call Error Logging to print any error messages.
!
!-----------------------------------------------------------------------

   call glc_ErrorPrint(ErrorCode)

!-----------------------------------------------------------------------
!
!  write final message to glc output log
!
!-----------------------------------------------------------------------
    write(stdout,*) '==================='
    write(stdout,*) 'completed glc_final'
    write(stdout,*) '==================='

!-----------------------------------------------------------------------
!
!  exit the communication environment
!
!-----------------------------------------------------------------------

   !call glc_CommExitEnvironment(ErrorCode)
   call exit_message_environment(ErrorCode)

!-----------------------------------------------------------------------
!EOC

 end subroutine glc_final

!***********************************************************************

 end module glc_FinalMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
