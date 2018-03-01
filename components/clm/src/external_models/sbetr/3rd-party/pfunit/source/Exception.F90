!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: PrivateException
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune, NASA/GSFC 
!!
!! @date
!! 07 Nov 2013
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 07 Nov 2013 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------
module PrivateException_mod
   use SourceLocation_mod
   implicit none
   private

   public :: Exception
   public :: newException
   public :: ExceptionList
   public :: newExceptionList

   public :: MAXLEN_MESSAGE
   public :: MAXLEN_FILE_NAME
   public :: NULL_MESSAGE
   public :: UNKNOWN_LINE_NUMBER
   public :: UNKNOWN_FILE_NAME

   integer, parameter :: MAXLEN_MESSAGE = 80*15
   integer, parameter :: MAXLEN_FILE_NAME = 80
   character(len=*), parameter :: NULL_MESSAGE = ''

   type Exception
      character(len=MAXLEN_MESSAGE) :: message = NULL_MESSAGE
      type (SourceLocation) :: location = UNKNOWN_SOURCE_LOCATION
      logical :: nullFlag = .true.
   contains
      procedure :: getMessage
      procedure :: getLineNumber
      procedure :: getFileName
      procedure :: isNull
   end type Exception

   type (Exception), parameter :: NULL_EXCEPTION = Exception('NULL EXCEPTION', UNKNOWN_SOURCE_LOCATION, .true.)

   type ExceptionList
      type (Exception), allocatable :: exceptions(:)
   contains
      
      procedure :: getNumExceptions

      procedure :: catch_any
      procedure :: catchNext
      procedure :: gather
      procedure :: catch_message
      generic :: catch => catch_any
      generic :: catch => catch_message
      procedure :: getExceptions
      procedure :: noExceptions
      procedure :: anyExceptions
      procedure :: clearAll
      procedure, private :: deleteIthException

      generic :: throw => throwMessage
      generic :: throw => throwException

      procedure :: throwMessage
      procedure :: throwException
!TODO - NAG does not yet support FINAL keyword
!!$$      final :: delete
   end type ExceptionList

   interface newException
      module procedure Exception_
   end interface

contains

   type(Exception) function Exception_(message, location)
      character(len=*), optional, intent(in) :: message
      type (SourceLocation), optional, intent(in) :: location

      if (present(message)) then
         Exception_%message = trim(message)
      else
         Exception_%message = NULL_MESSAGE
      end if

      if (present(location)) then
         Exception_%location = location
      else
         Exception_%location = UNKNOWN_SOURCE_LOCATION
      end if

      Exception_%nullFlag = .false.

    end function Exception_

   function getMessage(this) result(message)
      class (Exception), intent(in) :: this
      character(len=len_trim(this%message)) :: message
      message = trim(this%message)
   end function getMessage

   integer function getLineNumber(this) 
      class (Exception), intent(in) :: this
      getLineNumber = this%location%lineNumber
   end function getLineNumber

   character(len=MAXLEN_FILE_NAME) function getFileName(this) 
      class (Exception), intent(in) :: this
      getFileName = trim(this%location%fileName)
   end function getFileName

   logical function isNull(this) 
      class (Exception), intent(in) :: this
      isNull = this%nullFlag
   end function isNull

   function newExceptionList() result(list)
      type (ExceptionList) :: list
      allocate(list%exceptions(0))
   end function newExceptionList

   integer function getNumExceptions(this)
      class (ExceptionList), intent(in) :: this
      getNumExceptions = size(this%exceptions)
   end function getNumExceptions

   subroutine throwMessage(this, message, location)
      class (ExceptionList), intent(inOut) :: this
      character(len=*), intent(in) :: message
      type (SourceLocation), optional, intent(in) :: location

      call this%throw(newException(message, location))

   end subroutine throwMessage

   subroutine throwException(this, anException)
      class (ExceptionList), intent(inOut) :: this
      type (Exception), intent(in) :: anException

      type (Exception), allocatable :: tmp(:)
      integer :: n

      n = size(this%exceptions)
      allocate(tmp(n+1))
      tmp(1:n) = this%exceptions
      tmp(n+1) = anException
      deallocate(this%exceptions)
      allocate(this%exceptions(n+1))
      this%exceptions = tmp
      deallocate(tmp)

   end subroutine throwException

   function catchNext(this, preserve) result(anException)
      class (ExceptionList), intent(inOut) :: this
      logical, optional, intent(in) :: preserve
      type (Exception) :: anException
      if (size(this%exceptions) > 0) then
         anException = this%exceptions(1)
         call this%deleteIthException(1, preserve)
      else
         anException = NULL_EXCEPTION
      end if

   end function catchNext

   subroutine gather(this, context)
      use ParallelContext_mod
      class (ExceptionList), intent(inOut) :: this
      class (ParallelContext), intent(in) :: context

      type (ExceptionList) :: list
      integer :: globalExceptionCount
!      character(len=MAXLEN_MESSAGE) :: msg
      integer :: i

      globalExceptionCount = context%sum(size(this%exceptions))

      if (globalExceptionCount > 0) then

         allocate(list%exceptions(globalExceptionCount))
         
         do i = 1, this%getNumExceptions()
            call context%labelProcess(this%exceptions(i)%message)
         end do

         call context%gather(this%exceptions(:)%nullFlag, list%exceptions(:)%nullFlag)
         call context%gather(this%exceptions(:)%location%fileName, list%exceptions(:)%location%fileName)
         call context%gather(this%exceptions(:)%location%lineNumber, list%exceptions(:)%location%lineNumber)
         call context%gather(this%exceptions(:)%message, list%exceptions(:)%message)

         call clearAll(this)
      
         if (context%isRootProcess()) then
            deallocate(this%exceptions)
            allocate(this%exceptions(globalExceptionCount))
            this%exceptions(:) = list%exceptions
         end if

         call clearAll(list)

      end if

   end subroutine gather

   logical function noExceptions(this)
      class (ExceptionList), intent(inOut) :: this

      noExceptions = .not. this%anyExceptions()

   end function noExceptions

   logical function anyExceptions(this)
      class (ExceptionList), intent(inOut) :: this

      anyExceptions = (this%getNumExceptions() > 0)

   end function anyExceptions

   ! Fortran does not require "short-circuit" so be careful with 
   ! evaluation of optional arguments.
   logical function preserveMessage(preserve)
      logical, optional, intent(in) :: preserve

      preserveMessage = .false. ! default
      if (present(preserve)) preserveMessage = preserve

   end function preserveMessage

   subroutine deleteIthException(this, i, preserve)
      class (ExceptionList), intent(inOut) :: this
      integer, intent(in) :: i
      logical, optional, intent(in) :: preserve

      type (Exception), allocatable :: tmp(:)
      integer :: n

      if (preserveMessage(preserve)) return

      n = this%getNumExceptions()
      if (n == 0) return ! cannot throw exceptions here, alas
      allocate(tmp(n-1))
      tmp(1:i-1) = this%exceptions(1:i-1)
      tmp(i:n-1) = this%exceptions(i+1:n)
      deallocate(this%exceptions)
      allocate(this%exceptions(n-1))
      this%exceptions = tmp
      deallocate(tmp)

   end subroutine deleteIthException

   logical function catch_any(this, preserve) 
      class (ExceptionList), intent(inOut) :: this
      logical, optional, intent(in) :: preserve

      integer :: n
      logical :: preserve_ ! for default value

      n = this%getNumExceptions()

      if (n >= 1) then
         catch_any =.true.
         preserve_ = .false.
         if (present(preserve)) preserve_ = preserve
         if (.not. preserve_) call this%deleteIthException(n, preserve)
         return
      end if

      catch_any =.false.

   end function catch_any

   logical function catch_message(this, message, preserve) 
      class (ExceptionList), intent(inOut) :: this
      character(len=*), intent(in) :: message
      logical, optional, intent(in) :: preserve

      integer :: i, n
      logical :: preserve_ ! for default value

      n = this%getNumExceptions()

      do i = 1, n
         if (trim(message) == this%exceptions(i)%getMessage()) then
            catch_message =.true.
            preserve_ = .false.
            if (present(preserve)) preserve_ = preserve
            if (.not. preserve_) call this%deleteIthException(i, preserve)
            return
         end if
      end do
      catch_message =.false.

   end function catch_message

   function getExceptions(this) result(exceptions)
      type (Exception), allocatable :: exceptions(:)
      class (ExceptionList), intent(inOut) :: this

      call move_alloc(from=this%exceptions, to=exceptions)
      allocate(this%exceptions(0))
      
   end function getExceptions

   subroutine clearAll(this)
      class (ExceptionList), intent(inOut) :: this
      deallocate(this%exceptions)
      allocate(this%exceptions(0))
   end subroutine clearAll

   subroutine delete(this)
      type (ExceptionList), intent(inOut) :: this
      if (allocated(this%exceptions)) deallocate(this%exceptions)
   end subroutine delete

end module PrivateException_mod

module Exception_mod
   use SourceLocation_mod
   use PrivateException_mod
   implicit none
   private

   public :: Exception
   public :: newException
   public :: ExceptionList
   public :: newExceptionList

   public :: MAXLEN_MESSAGE
   public :: NULL_MESSAGE
   public :: UNKNOWN_LINE_NUMBER
   public :: UNKNOWN_FILE_NAME

   public :: getNumExceptions
   public :: throw
   public :: gatherExceptions
   public :: catchNext
   public :: catch
   public :: getExceptions
   public :: noExceptions
   public :: anyExceptions
   public :: anyErrors
   public :: clearAll

   public :: initializeGlobalExceptionList

   type (ExceptionList) :: globalExceptionList ! private

  interface throw
    module procedure throw_message
  end interface

  interface catch
     module procedure catch_any
     module procedure catch_message
  end interface catch

  interface anyExceptions
     module procedure anyExceptions_local
  end interface anyExceptions

  interface getNumExceptions
     module procedure getNumExceptions_local
  end interface getNumExceptions

contains

   subroutine initializeGlobalExceptionList()
      globalExceptionList = newExceptionList()
   end subroutine initializeGlobalExceptionList


   integer function getNumExceptions_local() result(numExceptions)

      if (.not. allocated(globalExceptionList%exceptions)) then
         call initializeGlobalExceptionList()
      end if

      numExceptions = globalExceptionList%getNumExceptions()
   end function getNumExceptions_local

   subroutine throw_message(message, location)
      character(len=*), intent(in) :: message
      type (SourceLocation), optional, intent(in) :: location

      !$omp critical
      if (.not. allocated(globalExceptionList%exceptions)) then
         call initializeGlobalExceptionList()
      end if

      call globalExceptionList%throw(message, location)
      !$omp end critical

   end subroutine throw_message

   function catchNext(preserve) result(anException)
      logical, optional, intent(in) :: preserve
      type (Exception) :: anException

      if (.not. allocated(globalExceptionList%exceptions)) then
         call initializeGlobalExceptionList()
      end if

      anException = globalExceptionList%catchNext(preserve)
   end function catchNext

   logical function catch_any(preserve)
      logical, optional, intent(in) :: preserve

      if (.not. allocated(globalExceptionList%exceptions)) then
         call initializeGlobalExceptionList()
      end if

      catch_any = globalExceptionList%catch(preserve)
   end function catch_any

   logical function catch_message(message, preserve)
      character(len=*), intent(in) :: message
      logical, optional, intent(in) :: preserve

      if (.not. allocated(globalExceptionList%exceptions)) then
         call initializeGlobalExceptionList()
      end if

      catch_message = globalExceptionList%catch(message, preserve)
   end function catch_message

   function getExceptions() result(exceptions)
      type (Exception), allocatable :: exceptions(:)

      if (.not. allocated(globalExceptionList%exceptions)) then
         call initializeGlobalExceptionList()
      end if

      exceptions = globalExceptionList%getExceptions()
   end function getExceptions

   logical function noExceptions()

      if (.not. allocated(globalExceptionList%exceptions)) then
         call initializeGlobalExceptionList()
      end if

      noExceptions = globalExceptionList%noExceptions()
   end function noExceptions

   logical function anyExceptions_local() result(anyExceptions)

      if (.not. allocated(globalExceptionList%exceptions)) then
         call initializeGlobalExceptionList()
      end if

      anyExceptions = globalExceptionList%anyExceptions()
   end function anyExceptions_local

   logical function anyErrors()
      integer :: i
      integer :: n

      do i = 1, globalExceptionList%getNumExceptions()
         n = min(14,len(globalExceptionList%exceptions(i)%message))
         if (globalExceptionList%exceptions(i)%message(1:n) == 'RUNTIME-ERROR:') then
            anyErrors = .true.
            return
         end if
      end do
      anyErrors = .false.
   end function anyErrors

   subroutine gatherExceptions(context)
      use ParallelContext_mod
      class (ParallelContext), intent(in) :: context
      call globalExceptionList%gather(context)
   end subroutine gatherExceptions

   subroutine clearAll()
      call globalExceptionList%clearAll()
   end subroutine clearAll

end module Exception_mod
