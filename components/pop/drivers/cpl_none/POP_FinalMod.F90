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
!  SVN:$Id: POP_FinalMod.F90 8528 2008-01-15 01:49:19Z dennis $
!
! !USES:

   use POP_KindsMod
   use POP_ErrorMod
   use POP_IOUnitsMod, only: POP_stdout
   use communicate
   use timers, only: timer_print_all
   !use POP_CommMod
   !use esmf_mod

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

   call POP_ErrorPrint(ErrorCode)

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
      write(POP_stdout,*) '==================='
      write(POP_stdout,*) 'completed POP_Final'
      write(POP_stdout,*) '==================='
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
