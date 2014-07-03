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
!  SVN:$Id: exit_mod.F90 35325 2012-03-09 00:48:12Z njn01 $

! !USES:

   use kinds_mod
   use communicate
   use constants
   use POP_IOUnitsMod

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

 subroutine exit_POP(exit_mode, exit_message, out_unit)

! !DESCRIPTION:
!  This routine prints a message, exits any message environment
!  and cleans up before stopping

! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (int_kind), intent(in) :: &
     exit_mode    ! method for exiting (normal exit or abort)

   integer (int_kind), optional, intent(in) :: &
     out_unit    ! optional output unit specifier

   character (*), intent(in) :: &
     exit_message ! message to print before stopping

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (int_kind) ::  &
      ierr,               &! error flag
      local_unit

!-----------------------------------------------------------------------
!
!  print message - must use unit 6 in place of stdout to
!  prevent circular dependence with io module
!
!-----------------------------------------------------------------------

   if (present(out_unit)) then
     local_unit = out_unit
   else
     local_unit = 6
   endif

#ifndef CCSMCOUPLED
   if (my_task == master_task) then
#endif
      write (local_unit,delim_fmt)
      write (local_unit,blank_fmt)
      call POP_IOUnitsFlush(local_unit)

      select case(exit_mode)
      case(sigExit)
         write (local_unit,'(a14)') 'POP exiting...'
      case(sigAbort)
         write (local_unit,'(a15)') 'POP aborting...'
      case default
         write (local_unit,'(a37)') 'POP exiting with unknown exit mode...'
      end select

      write (local_unit,*) exit_message
      write (local_unit,blank_fmt)
      write (local_unit,delim_fmt)
#ifndef CCSMCOUPLED
   endif
#endif
   call POP_IOUnitsFlush(local_unit)
   call POP_IOUnitsFlush(6)

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
   call POP_IOUnitsFlush(iunit)
 endif

 end subroutine flushm

!***********************************************************************

 end module exit_mod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
