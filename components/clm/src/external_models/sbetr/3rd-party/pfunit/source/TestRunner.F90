!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: TestRunner
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

module TestRunner_mod
   use Test_mod
   use BaseTestRunner_mod
   use TestListener_mod
   implicit none
   private

   public :: TestRunner
   public :: newTestRunner

   type, extends(BaseTestRunner) :: TestRunner
      type (ListenerPointer), allocatable :: extListeners(:)
   contains
      procedure :: run
      procedure :: createTestResult
      procedure :: addFailure
      procedure :: startTest
      procedure :: endTest
      procedure :: endRun
   end type TestRunner

   interface newTestRunner
      module procedure newTestRunner_default
      module procedure newTestRunner_unit
   end interface

contains

   function newTestRunner_default() result(runner)
!mlr-      use iso_fortran_env, only: OUTPUT_UNIT
      type (TestRunner) :: runner
      allocate(runner%extListeners(0))
   end function newTestRunner_default

   function newTestRunner_unit(extListeners) result(runner)
      type(ListenerPointer), intent(in) :: extListeners(:)
      type (TestRunner) :: runner
      allocate(runner%extListeners(size(extListeners)), source=extListeners)
   end function newTestRunner_unit

   function createTestResult(this) result(tstResult)
      use TestResult_mod
      class (TestRunner), intent(inout) :: this
      type (TestResult) :: tstResult

      tstResult = newTestResult()

    end function createTestResult

    function run(this, aTest, context) result(result)
      use Test_mod
      use TestSuite_mod
      use TestCase_mod
      use TestResult_mod
      use ParallelContext_mod

      type (TestResult) :: result
      class (TestRunner), intent(inout) :: this
      class (Test), intent(inout) :: aTest
      class (ParallelContext), intent(in) :: context
      
      integer :: clockStart
      integer :: clockStop
      integer :: clockRate
      integer :: i


      call system_clock(clockStart)

      result = this%createTestResult()
      call result%setName(aTest%getName())
! Add the extListeners to the listeners list.

      do i=1,size(this%extListeners)
         call result%addListener(this%extListeners(i)%pListener)
      end do
      call aTest%run(result, context)
      call system_clock(clockStop, clockRate)

      call result%addRunTime(real(clockStop - clockStart) / clockRate)

! Post run printing. Q: Should we do this for listeners too?
! E.g. and end-run method & move this up to basetestrunner...

! e.g. call result%endRun()...
      if (context%isRootProcess())  then
         do i=1,size(this%extListeners)
            call this%extListeners(i)%pListener%endRun(result)
         end do
      end if
!tc: 2+1 lists -- extListeners, listeners and in testresult too...
   end function run

! Recall, runner is also a listener and these will be called from
! TestResult, adding the ability to put in functionality here. In
! addition to other listeners added above.
!
    subroutine startTest(this, testName)
       class (TestRunner), intent(inout) :: this
       character(len=*), intent(in) :: testName
    end subroutine startTest

    subroutine endTest(this, testName)
       class (TestRunner), intent(inout) :: this
       character(len=*), intent(in) :: testName
    end subroutine endTest

    subroutine endRun(this, result)
      use AbstractTestResult_mod, only : AbstractTestResult
      class (TestRunner), intent(inout) :: this
      class (AbstractTestResult), intent(in) :: result
    end subroutine endRun

    subroutine addFailure(this, testName, exceptions)
       use Exception_mod
       class (TestRunner), intent(inout) :: this
       character(len=*), intent(in) :: testName
       type (Exception), intent(in) :: exceptions(:)
    end subroutine addFailure

end module TestRunner_mod
