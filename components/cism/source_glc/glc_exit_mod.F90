!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module glc_exit_mod

!BOP
! !MODULE: glc_exit_mod
!
! !DESCRIPTION:
!  This module provides a means for a graceful exit from glc when
!  encountering an error.  it contains only the routines exit\_glc
!  and flushm
!
! !REVISION HISTORY:
!  SVN:$Id: exit_mod.F90 808 2006-04-28 17:06:38Z njn01 $
!  Adapted by William Lipscomb from exit_mod.F90 in POP

! !USES:

   use glc_kinds_mod
   use glc_communicate
   use glc_constants
   use shr_sys_mod


   implicit none
   private
   save

! !PUBLIC MEMBER FUNCTIONS:

   public :: exit_glc, flushm

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
! !IROUTINE: exit_glc
! !INTERFACE:

 subroutine exit_glc(exit_mode, exit_message)

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

   if (my_task == master_task) then
      write (stdout,delim_fmt)
      write (stdout,blank_fmt)
      call shr_sys_flush(stdout)

      select case(exit_mode)
      case(sigExit)
         write (stdout,'(a14)') 'glc exiting...'
      case(sigAbort)
         write (stdout,'(a15)') 'glc aborting...'
      case default
         write (stdout,'(a37)') 'glc exiting with unknown exit mode...'
      end select

      write (stdout,*) exit_message
      write (stdout,blank_fmt)
      write (stdout,delim_fmt)
      call shr_sys_flush(stdout)
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

 end subroutine exit_glc

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

 end module glc_exit_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
