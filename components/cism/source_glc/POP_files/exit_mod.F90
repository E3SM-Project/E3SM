!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module exit_mod

!BOP
! !MODULE: exit_mod
!
! !DESCRIPTION:
!  This module provides a means for a graceful exit from POP when
!  encountering an error.  it contains only the routines exit\_POP
!  and flushm
!
! !REVISION HISTORY:
!  SVN:$Id: exit_mod.F90 808 2006-04-28 17:06:38Z njn01 $

! !USES:

   use kinds_mod
   use communicate
   use constants
   use shr_sys_mod


   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: exit_POP, flushm

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public :: &
      sigExit  =  0,    &! signal for normal exit
      sigAbort = -1      ! signal for aborting (exit due to error)

!EOP
!BOC
!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: exit_POP
! !INTERFACE:

 subroutine exit_POP(exit_mode, exit_message)

! !DESCRIPTION:
!  This routine prints a message, exits any message environment
!  and cleans up before stopping

! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     exit_mode    ! method for exiting (normal exit or abort)

   character (*), intent(in) :: &
     exit_message ! message to print before stopping

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) :: ierr  ! error flag

!-----------------------------------------------------------------------
!
!  print message - must use unit 6 in place of stdout to
!  prevent circular dependence with io module
!
!-----------------------------------------------------------------------

   if (my_task == master_task) then
      write (6,delim_fmt)
      write (6,blank_fmt)
      call shr_sys_flush(6)

      select case(exit_mode)
      case(sigExit)
         write (6,'(a14)') 'POP exiting...'
      case(sigAbort)
         write (6,'(a15)') 'POP aborting...'
      case default
         write (6,'(a37)') 'POP exiting with unknown exit mode...'
      end select

      write (6,*) exit_message
      write (6,blank_fmt)
      write (6,delim_fmt)
      call shr_sys_flush(6)
   endif

!-----------------------------------------------------------------------
!
!  exit or abort the message-passing environment if required
!
!-----------------------------------------------------------------------

   select case(exit_mode)
   case(sigExit)
      call exit_message_environment(ierr)
   case(sigAbort)
      call abort_message_environment(ierr)
   case default
   end select

!-----------------------------------------------------------------------
!
!  now we can stop
!
!-----------------------------------------------------------------------

   stop

!-----------------------------------------------------------------------
!EOC

 end subroutine exit_POP

!***********************************************************************
!BOP
! !IROUTINE: flushm (iunit)
! !INTERFACE:

 subroutine flushm (iunit)

! !DESCRIPTION:
!  This routine flushes the stdout buffer for the master_task only
!
! !REVISION HISTORY:
!  same as module
 integer (int_kind), intent(in) :: iunit
  
 if (my_task == master_task) then
   call shr_sys_flush (iunit)
 endif

 end subroutine flushm

!***********************************************************************

 end module exit_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
