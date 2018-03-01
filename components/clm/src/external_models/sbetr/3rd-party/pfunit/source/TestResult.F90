!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: TestResult
!
!> @brief
!! <BriefDescription>
!! Note: A possible extension point for user-specialized TestResults.
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
module TestResult_mod
   use AbstractTestResult_mod
   use SurrogateTestCase_mod
   use TestListener_mod
   use TestFailure_mod

   implicit none
   private

   public :: TestResult
   public :: newTestResult

#ifndef DEFERRED_LENGTH_CHARACTER
   integer, parameter :: MAX_LENGTH_NAME = 64
#endif

   type, extends(AbstractTestResult) :: TestResult
      private
      integer :: numFailed = 0
      integer :: numErrors = 0
      integer :: numRun = 0
      integer :: numSuccesses
      real    :: runTime
      type (ListenerPointer), allocatable :: listeners(:)
      type (TestFailure), allocatable :: failures(:)
      type (TestFailure), allocatable :: errors(:)
      type (TestFailure), allocatable :: successes(:)
#ifdef DEFERRED_LENGTH_CHARACTER
      character(:), allocatable :: name
#else
      character(len=MAX_LENGTH_NAME) :: name
#endif
   contains
      procedure :: addFailure
      procedure :: addError
      procedure :: addSuccess
      procedure :: zeroRunTime
      procedure :: addRunTime
      procedure :: getRunTime
      procedure :: failureCount
      procedure :: errorCount
      procedure :: startTest
      procedure :: endTest
      procedure :: runCount
      procedure :: run
      procedure :: addListener
      procedure :: wasSuccessful
      procedure :: getSuccesses
      procedure :: getErrors
      procedure :: getFailures
      procedure :: getIthFailure
      procedure :: getName
      procedure :: setName
   end type TestResult

contains

   function newTestResult(name)
      type (TestResult) :: newTestResult
      character(len=*), intent(in), optional :: name
      allocate(newTestResult%listeners(0))
      allocate(newTestResult%failures(0))
      allocate(newTestResult%errors(0))
      allocate(newTestResult%successes(0))
      newTestResult%numFailed = 0
      newTestResult%numErrors = 0
      newTestResult%numRun = 0
      newTestResult%numSuccesses = 0
      newTestResult%runTime = 0
      if(present(name)) then
         newTestResult%name = name
      else
         newTestResult%name = 'default_suite_name'
      end if
   end function newTestResult

   subroutine addFailure(this, aTest, exceptions)
      use Exception_mod, only: Exception
      use TestFailure_mod
      class (TestResult), intent(inout) :: this
      class (SurrogateTestCase), intent(in) :: aTest
      type (Exception), intent(in) :: exceptions(:)

      integer :: i, n
      type (TestFailure), allocatable :: tmp(:)

      n = this%numFailed
      allocate(tmp(n))
      tmp(1:n) = this%failures(1:n)
      deallocate(this%failures)
      allocate(this%failures(n+1))
      this%failures(1:n) = tmp
      deallocate(tmp)
      this%failures(n+1) = TestFailure(aTest%getName(), exceptions)

      this%numFailed = n + 1
      do i = 1, size(this%listeners)
         call this%listeners(i)%pListener%addFailure(aTest%getName(), exceptions)
      end do

   end subroutine addFailure

   subroutine addError(this, aTest, exceptions)
      use Exception_mod, only: Exception
      use TestFailure_mod
      class (TestResult), intent(inout) :: this
      class (SurrogateTestCase), intent(in) :: aTest
      type (Exception), intent(in) :: exceptions(:)

      integer :: i, n
      type (TestFailure), allocatable :: tmp(:)

      n = this%numErrors
      allocate(tmp(n))
      tmp(1:n) = this%errors(1:n)
      deallocate(this%errors)
      allocate(this%errors(n+1))
      this%errors(1:n) = tmp
      deallocate(tmp)
      this%errors(n+1) = TestFailure(aTest%getName(), exceptions)

      this%numErrors = n + 1
      do i = 1, size(this%listeners)
         call this%listeners(i)%pListener%addError(aTest%getName(), exceptions)
      end do

   end subroutine addError

   subroutine addSuccess(this, aTest)
      use Exception_mod, only: Exception
      use TestFailure_mod
      class (TestResult), intent(inout) :: this
      class (SurrogateTestCase), intent(in) :: aTest

!      integer :: i, n
      integer :: n
      type (TestFailure), allocatable :: tmp(:)
      type (Exception), allocatable :: noExceptions(:)
!      type (Exception) :: noExceptions(0)
      
      allocate(noExceptions(0))

      n = this%numSuccesses
      allocate(tmp(n))
      tmp(1:n) = this%successes(1:n)
      deallocate(this%successes)
      allocate(this%successes(n+1))
      this%successes(1:n) = tmp
      deallocate(tmp)
      this%successes(n+1) = TestFailure(aTest%getName(), noExceptions)

      this%numSuccesses = n + 1

   end subroutine addSuccess

   integer function failureCount(this)
      class (TestResult), intent(in) :: this
      failureCount = this%numFailed
   end function failureCount

   integer function errorCount(this)
      class (TestResult), intent(in) :: this
      errorCount = this%numErrors
   end function errorCount

   subroutine startTest(this, aTest)
     use StringConversionUtilities_mod, only : toString
      class (TestResult), intent(inout) :: this
      class (SurrogateTestCase), intent(in) :: aTest

      integer :: i


      this%numRun = this%numRun + 1

      !print *,'1000 starting: '//toString(this%numRun)//' '//trim(aTest%getName())

      do i = 1, size(this%listeners)
         call this%listeners(i)%pListener%startTest(aTest%getName())
      end do

   end subroutine startTest

   subroutine endTest(this, aTest)
      class (TestResult), intent(inout) :: this
      class (SurrogateTestCase), intent(in) :: aTest

      integer :: i

      do i = 1, size(this%listeners)
         call this%listeners(i)%pListener%endTest(aTest%getName())
      end do
   end subroutine endTest

   integer function runCount(this)
      class (TestResult), intent(in) :: this
      runCount = this%numRun
   end function runCount

! only invoked for a "real" test, not suites etc.
   recursive subroutine run(this, test, context)
      use Exception_mod
      use ParallelContext_mod
      class (TestResult), intent(inout) :: this
      class (SurrogateTestCase) :: test 
      class (ParallelContext), intent(in) :: context

!      type (Exception) :: anException

      if (context%isRootProcess()) call this%startTest(test)

      call test%runBare()

      if (context%isRootProcess()) then
         if (anyErrors()) then
            call this%addError(test, getExceptions())
         elseif (anyExceptions()) then
            call this%addFailure(test, getExceptions())
         else
            call this%addSuccess(test)
         end if
      end if

      if (context%isRootProcess()) call this%endTest(test)

   end subroutine run

   subroutine addListener(this, listener)
      use TestListener_mod, only: TestListener
      class (TestResult), intent(inOut) :: this
      class (TestListener), target, intent(in) :: listener

      integer :: n

      call extend(this%listeners)
      n = size(this%listeners)
      this%listeners(n)%pListener => listener

   contains

      subroutine extend(listeners)
         type (ListenerPointer), allocatable, intent(inout) :: listeners(:)
         type (ListenerPointer), allocatable :: temp(:)
         integer :: n

         n = size(listeners)
         temp = listeners
         deallocate(listeners)

         allocate(listeners(n+1))
         listeners(:n) = temp
         deallocate(temp)

      end subroutine extend

   end subroutine addListener

   function getIthFailure(this, i) result(failure)
      class (TestResult), intent(in) :: this
      integer, intent(in) :: i
      type (TestFailure) :: failure

      failure = this%failures(i)

   end function getIthFailure

   logical function wasSuccessful(this)
      class (TestResult), intent(in) :: this
      wasSuccessful = (this%failureCount() ==  0) .and. (this%errorCount() == 0)
   end function wasSuccessful

   function getSuccesses(this) result(successes)
     class (TestResult), intent(in) :: this
     type (TestFailure), allocatable :: successes(:)
     successes = this%successes
   end function getSuccesses

   function getErrors(this) result(errors)
      class (TestResult), intent(in) :: this
      type (TestFailure), allocatable :: errors(:)
      errors = this%errors
   end function getErrors

   function getFailures(this) result(failures)
      class (TestResult), intent(in) :: this
      type (TestFailure), allocatable :: failures(:)
      failures = this%failures
   end function getFailures

   subroutine zeroRunTime(this)
     class (TestResult), intent(inout) :: this
     this%runTime = 0
   end subroutine zeroRunTime

   subroutine addRunTime(this, time)
     class (TestResult), intent(inout) :: this
     real, intent(in) :: time
     this%runTime = this%runTime + time
   end subroutine addRunTime

   function getRunTime(this) result(duration)
     class (TestResult), intent(in) :: this
     real :: duration
     duration = this%runTime
   end function getRunTime

   function getName(this) result(name)
      class (TestResult), intent(in) :: this
      character(:), allocatable :: name
      name = this%name
   end function getName

   subroutine setName(this, name)
      class (TestResult), intent(inout) :: this
      character(len=*),intent(in) :: name

      this%name = trim(name)
   end subroutine setName

end module TestResult_mod
