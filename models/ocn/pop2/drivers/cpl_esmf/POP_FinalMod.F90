!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_FinalMod

!BOP
! !MODULE: POP_FinalMod
! !DESCRIPTION:
!  This module contains the POP finalization method that shuts down POP
!  gracefully (we hope).  It exits the message environment and checks 
!  for successful execution.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  SVN:$Id$
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_IOUnitsMod, only: POP_stdout
   use communicate
   use io_types
   use timers, only: timer_print_all

   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_Final

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
! !IROUTINE: POP_Final
! !INTERFACE:

 subroutine POP_Final(ErrorCode)

! !DESCRIPTION:
!  This routine shuts down POP by exiting all relevent environments.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !INPUT/OUTPUT PARAMETERS:

   integer (POP_i4), intent(inout) :: &
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
!  call Error Logging to print any error messages.
!
!-----------------------------------------------------------------------

   call POP_ErrorPrint(errorCode, printTask=master_task)

!-----------------------------------------------------------------------
!
!  clear any open displays and print all timers with statistics
!
!-----------------------------------------------------------------------

   call timer_print_all(stats=.true.)

!-----------------------------------------------------------------------
!
!  write final message to pop output log
!
!-----------------------------------------------------------------------
    if (my_task == master_task) then
      write(stdout,*) '==================='
      write(stdout,*) 'completed POP_Final'
      write(stdout,*) '==================='
    endif

!-----------------------------------------------------------------------
!
!  exit the communication environment
!
!-----------------------------------------------------------------------

   !call POP_CommExitEnvironment(ErrorCode)
   call exit_message_environment(ErrorCode)

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_Final

!***********************************************************************

 end module POP_FinalMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
