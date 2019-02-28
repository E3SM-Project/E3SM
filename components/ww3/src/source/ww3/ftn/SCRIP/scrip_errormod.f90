!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module SCRIP_ErrorMod

!BOP
! !MODULE: SCRIP_ErrorMod
! !DESCRIPTION:
! This module contains SCRIP error flags and facilities for logging and
! printing error messages. Note that error flags are local to a
! process and there is no synchronization of error flags across
! processes. As routines trap error flags, they may add a message
! to the error log to aid in tracking the call sequence.
!
! Users should not need to change any values in this module.
!
! All routines in SCRIP which encounter an error should return to
! the calling routine with the SCRIP\_Fail error code set and a message
! added to the error log using the SCRIP\_ErrorCheck or
! SCRIP\_ErrorSet function. Routines in SCRIP should also check
! error codes returned by called routines and add a message to the
! error log to help users track the calling sequence that generated
! the error. This process enables the error code to be propagated
! to the highest level or calling routine to enable a graceful
! exit. At that level, the SCRIP_ErrorPrint call can be used to
! print the entire error trace or error log.
!
! !REVISION HISTORY:
! SVN:$Id: SCRIP_ErrorMod.F90 14 2006-08-17 17:07:05Z $
!
! !USES:

   use SCRIP_KindsMod
   !use SCRIP_CommMod
   use SCRIP_IOUnitsMod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer (SCRIP_i4), parameter, public :: &
      SCRIP_Success = 0, & ! standard SCRIP error flags
      SCRIP_Fail = -1

! !PUBLIC MEMBER FUNCTIONS:

   public :: SCRIP_ErrorSet, &
             SCRIP_ErrorCheck, &
             SCRIP_ErrorPrint

!EOP
!BOC
!-----------------------------------------------------------------------
!
! module variables
!
!-----------------------------------------------------------------------

   integer (SCRIP_i4), parameter :: &
      SCRIP_errorLogDepth = 20 ! Max depth of call tree to properly
                               ! size the error log array

   integer (SCRIP_i4) :: &
      SCRIP_errorMsgCount = 0 ! tracks current number of log messages

   character (SCRIP_CharLength), dimension(SCRIP_ErrorLogDepth) :: &
      SCRIP_errorLog ! list of error messages to be output

!EOC
!***********************************************************************

 contains

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_ErrorSet -- sets error code and logs error message
! !INTERFACE:

 subroutine SCRIP_ErrorSet(errorCode, rtnName, errorMsg)

! !DESCRIPTION:
! This routine sets an error code to SCRIP\_Fail and adds a message to
! the error log for later printing.
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

   integer (SCRIP_i4), intent(out) :: &
      errorCode ! Error code to set to fail

! !INPUT PARAMETERS:

   character (*), intent(in) :: &
      rtnName, &! name of calling routine
      errorMsg ! message to add to error log for printing

!EOP
!BOC
!-----------------------------------------------------------------------
!
! Local variables
!
!-----------------------------------------------------------------------

   character(SCRIP_charLength) :: &
      logErrorMsg ! constructed error message with routine name

!-----------------------------------------------------------------------
!
! Set error code to fail
!
!-----------------------------------------------------------------------

   errorCode = SCRIP_Fail

!-----------------------------------------------------------------------
!
! Add error message to error log
!
!-----------------------------------------------------------------------

   SCRIP_errorMsgCount = SCRIP_errorMsgCount + 1

   if (SCRIP_errorMsgCount <= SCRIP_errorLogDepth) then
      write(logErrorMsg,'(a,a2,a)') rtnName,': ',errorMsg
      SCRIP_errorLog(SCRIP_errorMsgCount) = logErrorMsg
   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine SCRIP_ErrorSet

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_ErrorCheck -- checks error code and logs error message
! !INTERFACE:

 function SCRIP_ErrorCheck(errorCode, rtnName, errorMsg)

! !DESCRIPTION:
! This function checks an error code and adds a message to the error
! log for later printing. It is a more compact form of the ErrorSet
! routine that is especially useful for checking an error code after
! returning from a routine or function. If the errorCode is the
! failure code SCRIP\_Fail, it returns a logical true value so that
! it can be used in a typical call like:
! \begin{verbatim}
! if (SCRIP_ErrorCheck(errorCode, rtnName, errorMsg)) return
! \end{verbatim}
!
! !REVISION HISTORY:
! same as module

! !OUTPUT PARAMETERS:

   logical (SCRIP_logical) :: &
      SCRIP_ErrorCheck

! !INPUT PARAMETERS:

   integer (SCRIP_i4), intent(in) :: &
      errorCode ! Error code to check

   character (*), intent(in) :: &
      rtnName, &! name of calling routine
      errorMsg ! message to add to error log for printing

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   character (SCRIP_charLength) :: &
      logErrorMsg ! constructed error message with routine name

!-----------------------------------------------------------------------
!
! If the error code is success, set the return value to false.
!
!-----------------------------------------------------------------------

   if (errorCode == SCRIP_Success) then
      SCRIP_ErrorCheck = .false.

!-----------------------------------------------------------------------
!
! If the error code is a fail, set the return value to true and
! add the error message to the log.
!
!-----------------------------------------------------------------------

   else
      SCRIP_ErrorCheck = .true.

      SCRIP_errorMsgCount = SCRIP_errorMsgCount + 1

      if (SCRIP_errorMsgCount <= SCRIP_errorLogDepth) then
         write(logErrorMsg,'(a,a2,a)') rtnName,': ',errorMsg
         SCRIP_errorLog(SCRIP_errorMsgCount) = logErrorMsg
      endif
   endif

!-----------------------------------------------------------------------
!EOC

 end function SCRIP_ErrorCheck

!***********************************************************************
!BOP
! !IROUTINE: SCRIP_ErrorPrint -- prints the error log
! !INTERFACE:

 subroutine SCRIP_ErrorPrint(printTask)

! !DESCRIPTION:
! This routine prints all messages in the error log. If a printTask
! is specified, only the message log on that task will be printed.
!
! !REVISION HISTORY:
! same as module

! !INPUT PARAMETERS:

!   integer (SCRIP_i4), intent(in) :: &
!      errorCode ! input error code to check success/fail

   !*** currently this has no meaning, but will be used in parallel
   !*** SCRIP version
   integer (SCRIP_i4), intent(in), optional :: &
      printTask ! Task from which to print error log

!EOP
!BOC
!-----------------------------------------------------------------------
!
! local variables
!
!-----------------------------------------------------------------------

   integer (SCRIP_i4) :: n

!-----------------------------------------------------------------------
!
! Print all error messages to stdout
!
!-----------------------------------------------------------------------

   if (present(printTask)) then

      !*** parallel SCRIP not yet supported
      !if (SCRIP_myTask == printTask) then

         write(SCRIP_stdout,SCRIP_blankFormat)
         write(SCRIP_stdout,SCRIP_delimFormat)
         write(SCRIP_stdout,SCRIP_blankFormat)

         if (SCRIP_errorMsgCount == 0) then ! no errors

            write(SCRIP_stdout,'(a34)') &
                                'Successful completion of SCRIP model'

         else

            write(SCRIP_stdout,'(a14)') 'SCRIP Exiting...'

            do n=1,min(SCRIP_errorMsgCount,SCRIP_errorLogDepth)
               write(SCRIP_stderr,'(a)') trim(SCRIP_errorLog(n))
               if (SCRIP_stdout /= SCRIP_stderr) then
                  write(SCRIP_stdout,'(a)') trim(SCRIP_errorLog(n))
               endif
            end do

            if (SCRIP_errorMsgCount > SCRIP_errorLogDepth) then
               write(SCRIP_stderr,'(a23)') 'Too many error messages'
               if (SCRIP_stdout /= SCRIP_stderr) then
                  write(SCRIP_stdout,'(a23)') 'Too many error messages'
               endif
            endif

         endif

         write(SCRIP_stdout,SCRIP_blankFormat)
         write(SCRIP_stdout,SCRIP_delimFormat)
         write(SCRIP_stdout,SCRIP_blankFormat)

      !endif

   else

      write(SCRIP_stdout,SCRIP_blankFormat)
      write(SCRIP_stdout,SCRIP_delimFormat)
      write(SCRIP_stdout,SCRIP_blankFormat)

      if (SCRIP_errorMsgCount == 0) then ! no errors

         write(SCRIP_stdout,'(a34)') 'Successful completion of SCRIP'

      else

         write(SCRIP_stdout,'(a14)') 'SCRIP Exiting...'

         do n=1,min(SCRIP_errorMsgCount,SCRIP_errorLogDepth)
            write(SCRIP_stderr,'(a)') trim(SCRIP_errorLog(n))
            if (SCRIP_stdout /= SCRIP_stderr) then
               write(SCRIP_stdout,'(a)') trim(SCRIP_errorLog(n))
            endif
         end do

         if (SCRIP_errorMsgCount > SCRIP_errorLogDepth) then
            write(SCRIP_stderr,'(a23)') 'Too many error messages'
            if (SCRIP_stdout /= SCRIP_stderr) then
               write(SCRIP_stdout,'(a23)') 'Too many error messages'
            endif
         endif

      endif

      write(SCRIP_stdout,SCRIP_blankFormat)
      write(SCRIP_stdout,SCRIP_delimFormat)
      write(SCRIP_stdout,SCRIP_blankFormat)

   endif

!-----------------------------------------------------------------------
!EOC

 end subroutine SCRIP_ErrorPrint

!***********************************************************************

 end module SCRIP_ErrorMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
