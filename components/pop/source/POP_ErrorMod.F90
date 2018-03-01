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
!  Users should not need to change any values in this module.
!
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
!  SVN:$Id: POP_ErrorMod.F90 16641 2009-06-12 17:00:31Z njn01 $
!  2006-07-10: Phil Jones
!     Added new error module for logging and printing error messages.
!
! !USES:

   use POP_KindsMod
   use POP_CommMod
   use POP_IOUnitsMod
#ifdef CCSMCOUPLED
   use io_types
#endif

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
      POP_errorLogDepth = 20   ! Max depth of call tree to properly
                               ! size the error log array

   integer (POP_i4) ::         &
      POP_errorMsgCount =  0   ! tracks current number of log messages

   character (POP_CharLength), dimension(POP_ErrorLogDepth) :: &
      POP_errorLog             ! list of error messages to be output

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: POP_ErrorSet -- sets error code and logs error message
! !INTERFACE:

 subroutine POP_ErrorSet(errorCode, errorMsg)

! !DESCRIPTION:
!  This routine sets an error code to POP\_Fail and adds a message to
!  the error log for later printing.
!
! !REVISION HISTORY:
!  same as module

! !OUTPUT PARAMETERS:

   integer (POP_i4), intent(out) :: &
      errorCode              ! Error code to set to fail

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      errorMsg               ! message to add to error log for printing

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  Set error code to fail
!
!-----------------------------------------------------------------------

   errorCode = POP_Fail

!-----------------------------------------------------------------------
!
!  Add error message to error log
!
!-----------------------------------------------------------------------

   POP_errorMsgCount = POP_errorMsgCount + 1

   if (POP_errorMsgCount <= POP_errorLogDepth) then
      POP_errorLog(POP_errorMsgCount) = errorMsg
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ErrorSet

!***********************************************************************
!BOP
! !IROUTINE: POP_ErrorPrint -- prints the error log
! !INTERFACE:

 subroutine POP_ErrorPrint(errorCode, printTask)

! !DESCRIPTION:
!  This routine prints all messages in the error log.  If a printTask
!  is specified, only the message log on that task will be printed.
!
! !REVISION HISTORY:
!  same as module

! !INPUT PARAMETERS:

   integer (POP_i4), intent(in) :: &
      errorCode              ! input error code to check success/fail

   integer (POP_i4), intent(in), optional :: &
      printTask              ! Task from which to print error log

!EOP
!BOC
!-----------------------------------------------------------------------
!
!  local variables
!
!-----------------------------------------------------------------------

   integer (POP_i4) :: n, local_unit

!-----------------------------------------------------------------------
!
!  Print all error messages to stdout
!
!-----------------------------------------------------------------------

#ifdef CCSMCOUPLED
    ! unnecessary after fully converting to new pop2 infrastructure
    local_unit = stdout
#else
    local_unit = POP_stdout
#endif
   if (present(printTask)) then

      if (POP_myTask == printTask) then

         write(local_unit,POP_blankFormat)
         write(local_unit,POP_delimFormat)
         write(local_unit,POP_blankFormat)

         if (POP_errorMsgCount == 0) then ! no errors

            write(local_unit,'(a34)') &
                                'Successful completion of POP model'

         else

            write(local_unit,'(a14)') 'POP Exiting...'

            do n=1,min(POP_errorMsgCount,POP_errorLogDepth)
               write(POP_stderr,'(a)') trim(POP_errorLog(n))
               if (local_unit /= POP_stderr) then
                  write(local_unit,'(a)') trim(POP_errorLog(n))
               endif
            end do

            if (POP_errorMsgCount > POP_errorLogDepth) then
               write(POP_stderr,'(a23)') 'Too many error messages'
               if (local_unit /= POP_stderr) then
                  write(local_unit,'(a23)') 'Too many error messages'
               endif
            endif

         endif

         write(local_unit,POP_blankFormat)
         write(local_unit,POP_delimFormat)
         write(local_unit,POP_blankFormat)

      endif

   else

      write(local_unit,POP_blankFormat)
      write(local_unit,POP_delimFormat)
      write(local_unit,POP_blankFormat)

      if (POP_errorMsgCount == 0) then ! no errors

         write(local_unit,'(a34)') 'Successful completion of POP model'

      else

         write(local_unit,'(a14)') 'POP Exiting...'

         do n=1,min(POP_errorMsgCount,POP_errorLogDepth)
            write(POP_stderr,'(a)') trim(POP_errorLog(n))
            if (local_unit /= POP_stderr) then
               write(local_unit,'(a)') trim(POP_errorLog(n))
            endif
         end do

         if (POP_errorMsgCount > POP_errorLogDepth) then
            write(POP_stderr,'(a23)') 'Too many error messages'
            if (local_unit /= POP_stderr) then
               write(local_unit,'(a23)') 'Too many error messages'
            endif
         endif

      endif

      write(local_unit,POP_blankFormat)
      write(local_unit,POP_delimFormat)
      write(local_unit,POP_blankFormat)

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine POP_ErrorPrint

!***********************************************************************

 end module POP_ErrorMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
