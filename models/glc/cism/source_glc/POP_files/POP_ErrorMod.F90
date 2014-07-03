!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_ErrorMod

!BOP
! !MODULE: POP_ErrorMod
! !DESCRIPTION:
!  This module contains POP error flags and facilities for logging and
!  printing error messages.  Note that error flags are local to a 
!  process and there is no synchronization of error flags across 
!  processes.  As routines trap error flags, they may add a message
!  to the error log to aid in tracking the call sequence.
!
! !USERDOC:
!  Users should not need to change any values in this module.
!
! !REFDOC:
!  All routines in POP which encounter an error should return to
!  the calling routine with the POP\_Fail error code set and a message
!  added to the error log using the POP\_ErrorSet function.  Also,
!  routines in POP should check error codes returned by called routines
!  and add a message to the error log to help users track the calling
!  sequence that generated the error.  This process
!  enables the error code to be propagated to the highest level or
!  to a coupler for a proper call to the POP finalize method.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_ErrorMod.F90 808 2006-04-28 17:06:38Z njn01 $
!
! !USES:

   use POP_KindsMod
   !use POP_CommMod
   use communicate
   use constants
   use POP_IOUnitsMod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer (POP_i4), parameter, public :: &
      POP_Success =  0,           & ! standard POP error flags
      POP_Fail    = -1

! !PUBLIC MEMBER FUNCTIONS:

   public :: POP_ErrorSet, &
             POP_ErrorPrint

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  module variables
!
!-----------------------------------------------------------------------

   integer (POP_i4), parameter :: &
      POP_ErrorLogDepth = 20   ! Max depth of call tree to properly
                               ! size the error log array

   integer (POP_i4) ::         &
      POP_ErrorMsgCount =  0   ! tracks current number of log messages

   character (POP_CharLength), dimension(POP_ErrorLogDepth) :: &
      POP_ErrorLog             ! list of error messages to be output

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_ErrorSet -- sets error code and logs error message
! !INTERFACE:

 subroutine POP_ErrorSet(ErrorCode, ErrorMsg)

! !DESCRIPTION:
!  This routine sets an error code to POP\_Fail and adds a message to 
!  the error log for later printing.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      ErrorCode              ! Error code to set to fail

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      ErrorMsg               ! message to add to error log for printing

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  Set error code to fail
!
!-----------------------------------------------------------------------

   ErrorCode = POP_Fail

!-----------------------------------------------------------------------
!
!  Add error message to error log
!
!-----------------------------------------------------------------------

   POP_ErrorMsgCount = POP_ErrorMsgCount + 1

   if (POP_ErrorMsgCount <= POP_ErrorLogDepth) then
      POP_ErrorLog(POP_ErrorMsgCount) = ErrorMsg
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ErrorSet

!***********************************************************************
!BOP
! !IROUTINE: POP_ErrorPrint -- prints the error log
! !INTERFACE:

 subroutine POP_ErrorPrint(ErrorCode, PrintTask)

! !DESCRIPTION:
!  This routine prints all messages in the error log.  If a PrintTask
!  is specified, only the log on that task will be printed.
!
! !USERDOC:
!
! !REFDOC:
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      ErrorCode              ! input error code to check success/fail

   integer (POP_i4), intent(in), optional :: &
      PrintTask              ! Task from which to print error log

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: n

!-----------------------------------------------------------------------
!
!  Print all error messages to stdout
!
!-----------------------------------------------------------------------

   if (present(PrintTask)) then

      if (my_Task == PrintTask) then
      !if (POP_myTask == PrintTask) then

         write(POP_stdout,blank_fmt)
         write(POP_stdout,'(a34)') '----------------------------------'

         if (POP_ErrorMsgCount == 0) then ! no errors

            write(POP_stdout,'(a34)') &
                                'Successful completion of POP model'

         else

            write(POP_stdout,'(a14)') 'POP Exiting...'
            do n=1,min(POP_ErrorMsgCount,POP_ErrorLogDepth)
               write(POP_stderr,'(a)') trim(POP_ErrorLog(n))
               if (POP_stdout /= POP_stderr) then
                  write(POP_stdout,'(a)') trim(POP_ErrorLog(n))
               endif
            end do
            if (POP_ErrorMsgCount > POP_ErrorLogDepth) then
               write(POP_stderr,'(a)') 'Too many error messages'
               if (POP_stdout /= POP_stderr) then
                  write(POP_stdout,'(a)') 'Too many error messages'
               endif
            endif

         endif

         write(POP_stdout,'(a34)') '----------------------------------'

      endif

   else

      write(POP_stdout,'(a34)') '----------------------------------'

      if (POP_ErrorMsgCount == 0) then ! no errors

         write(POP_stdout,'(a34)') 'Successful completion of POP model'

      else

         write(POP_stdout,'(a14)') 'POP Exiting...'
         do n=1,min(POP_ErrorMsgCount,POP_ErrorLogDepth)
            write(POP_stderr,'(a)') trim(POP_ErrorLog(n))
            if (POP_stdout /= POP_stderr) then
               write(POP_stdout,'(a)') trim(POP_ErrorLog(n))
            endif
         end do
         if (POP_ErrorMsgCount > POP_ErrorLogDepth) then
            write(POP_stderr,'(a)') 'Too many error messages'
            if (POP_stdout /= POP_stderr) then
               write(POP_stdout,'(a)') 'Too many error messages'
            endif
         endif

      endif

      write(POP_stdout,'(a34)') '----------------------------------'

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ErrorPrint

!***********************************************************************

 end module POP_ErrorMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
