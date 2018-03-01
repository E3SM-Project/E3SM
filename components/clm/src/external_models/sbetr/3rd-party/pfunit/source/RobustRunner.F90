
!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: RobustRunner
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
module RobustRunner_mod
   use Test_mod
   use TestCase_mod
   use BaseTestRunner_mod
   use TestListener_mod
   use UnixProcess_mod
   implicit none
   private

   public :: RobustRunner
#ifndef DEFERRED_LENGTH_CHARACTER
   integer, parameter :: MAX_LENGTH_COMMAND=80
#endif

   type, extends(BaseTestRunner) :: RobustRunner
      private
#ifdef DEFERRED_LENGTH_CHARACTER
      character(len=:), allocatable :: remoteRunCommand
#else
      character(len=MAX_LENGTH_COMMAND) :: remoteRunCommand
#endif
      integer :: numSkip
      type (ListenerPointer), allocatable :: extListeners(:)
      type (UnixProcess) :: remoteProcess
      real :: maxLaunchDuration
      real :: maxTimeoutDuration
   contains
      procedure :: run
      procedure :: runWithResult
      procedure :: startTest
      procedure :: endTest
      procedure :: endRun
      procedure :: addFailure
      procedure :: addError
      procedure :: launchRemoteRunner
      procedure :: createTestResult
   end type RobustRunner

   interface RobustRunner
!      module procedure newRobustRunner
      module procedure newRobustRunner_extListeners
   end interface RobustRunner

   type, extends(TestCase) :: TestCaseMonitor
      private
      type (UnixProcess), pointer :: process
   contains
      procedure :: runMethod
   end type TestCaseMonitor

!!! Inject dependency through constructor...
   real, parameter :: MAX_TIME_LAUNCH = 5.00 ! in seconds
   real, parameter :: MAX_TIME_TEST   = 0.11 ! in seconds

contains

!   function newRobustRunner(remoteRunCommand,maxLaunchDuration) result(runner)
!      type (RobustRunner) :: runner
!      character(len=*), intent(in) :: remoteRunCommand
!      real, optional, intent(in) :: maxLaunchDuration
!
!      if(.not.present(maxLaunchDuration))then
!         runner%maxLaunchDuration = MAX_TIME_LAUNCH
!      else
!         runner%maxLaunchDuration = maxLaunchDuration
!      end if
!      
!      runner%remoteRunCommand = trim(remoteRunCommand)
!      allocate(runner%extListeners(0))
!   end function newRobustRunner

   function newRobustRunner_extListeners( &
        & remoteRunCommand   &
        & ,extListeners      &
        & ,maxLaunchDuration &
        & ,maxTimeoutDuration &
        & ) result(runner)
      type (RobustRunner) :: runner
      character(len=*), intent(in) :: remoteRunCommand
      type(ListenerPointer), optional, intent(in) :: extListeners(:)
      
      real, optional, intent(in) :: maxLaunchDuration
      real, optional, intent(in) :: maxTimeoutDuration

      if(.not.present(maxLaunchDuration))then
         runner%maxLaunchDuration = MAX_TIME_LAUNCH
      else
         runner%maxLaunchDuration = maxLaunchDuration
      end if

      if(.not.present(maxTimeoutDuration))then
         runner%maxTimeoutDuration = MAX_TIME_TEST
      else
         runner%maxTimeoutDuration = maxTimeoutDuration
      end if

      if(present(extListeners))then
         allocate(runner%extListeners(size(extListeners)), source=extListeners)
      end if
      
      runner%remoteRunCommand = trim(remoteRunCommand)
      runner%numSkip = 0
      
   end function newRobustRunner_extListeners

   subroutine runMethod(this)
      class (TestCaseMonitor), intent(inout) :: this
   end subroutine runMethod

   function run(this, aTest, context) result(result)
      use Test_mod
      use TestSuite_mod
      use TestResult_mod
      use ParallelContext_mod

      type (TestResult) :: result
      class (RobustRunner), intent(inout) :: this
      class (Test), intent(inout) :: aTest
      class (ParallelContext), intent(in) :: context

      result = this%createTestResult()
      call result%setName(aTest%getName())
      call this%runWithResult(aTest, context, result)

   end function run

   subroutine runWithResult(this, aTest, context, result)
      use Test_mod
      use ParallelContext_mod
      use TestResult_mod
      use RemoteProxyTestCase_mod
      use TestSuite_mod
      use Exception_mod
      class (RobustRunner), intent(inout) :: this
      class (Test), intent(inout) :: aTest
      class (ParallelContext), intent(in) :: context
      type (TestResult), intent(inout) :: result

      type (TestCaseReference), allocatable :: testCases(:)
      type (RemoteProxyTestCase) :: proxy
      integer :: i
      integer :: clockStart, clockStop, clockRate

      call system_clock(clockStart)

      do i=1,size(this%extListeners)
         call result%addListener(this%extListeners(i)%pListener)
      end do
      call result%addListener( this ) ! - monitoring

      select type (aTest)
      class is (TestSuite)
!!! if defined(PGI) || (defined(__INTEL_COMPILER) && (INTEL_13))
#if (defined(__INTEL_COMPILER) && (INTEL_13))         
         testCases = aTest%getTestCases()
#else
         call aTest%getTestCases(testCases)
#endif
      class is (TestCase)
         allocate(testCases(1))
         allocate(testCases(1)%test, source= aTest)
      class default
         stop
      end select

! mlr q: set up named pipes or units to handle comm between remote processes
      ! mlr q: and the root... being done at ukmet?
      do i = 1, size(testCases)
         if (.not. this%remoteProcess%isActive()) then
            call this%launchRemoteRunner(numSkip=i-1)
         end if
         proxy = RemoteProxyTestCase( &
              &     testCases(i)%test &
              &     ,this%remoteProcess &
              &     ,maxTimeoutDuration=this%maxTimeoutDuration &
              &  )
         call proxy%run(result, context)
      end do

      call system_clock(clockStop, clockRate)

      call result%addRunTime(real(clockStop - clockStart) / clockRate)

      ! Maybe push this call up into parent, i.e. loop over all of the listeners there...
      if (context%isRootProcess())  then
         do i=1,size(this%extListeners)
            call this%extListeners(i)%pListener%endRun(result)
         end do
      end if

   end subroutine runWithResult

   subroutine launchRemoteRunner(this, numSkip)
      use UnixProcess_mod
      use Exception_mod
      class (RobustRunner), intent(inout) :: this
      integer, intent(in) :: numSkip

      character(len=:), allocatable :: command

      integer, parameter :: MAX_LEN=8
      character(len=MAX_LEN) :: suffix

      character(len=80) :: timeCommand
      type (UnixProcess) :: timerProcess
      character(len=:), allocatable :: line
      character(len=100) :: throwMessage

      
      write(suffix,'(i0)') numSkip
      command = trim(this%remoteRunCommand) // ' -skip ' // suffix


      this%remoteProcess = UnixProcess(command, runInBackground=.true.)

      ! Check for successful launch - prevents MPI launch time from counting against
      ! first test's time limit.
      write(timeCommand,'(a, f10.3,a,i0,a)') &
           & "(sleep ",this%maxLaunchDuration," && kill -9 ", &
           & this%remoteProcess%getPid(), &
           & ") > /dev/null 2>&1"
      timerProcess = UnixProcess(trim(timeCommand), runInBackground=.true.)

      do
         line = this%remoteProcess%getLine()
         if (len(line) == 0) then
            if (.not. this%remoteProcess%isActive()) then
               write(throwMessage,'(a,f0.3,a)') &
                    & ' (max launch duration = ',this%maxLaunchDuration,')'
               call throw('RUNTIME-ERROR: terminated before starting'//trim(throwMessage))
               call timerProcess%terminate()
               return
            else
!!$               call timerProcess%terminate()
!!$               timerProcess = UnixProcess(trim(timeCommand), runInBackground=.true.)
               cycle ! might just not be ready yet
            end if
         else
            if ('*LAUNCHED*' /= line) then
               call throw(&
	       &    'Failure to launch in RobustRunner. ' &
	       &    //"Expected: '*LAUNCHED*' Found: '"//line//"'" )
               return
            else
               ! successfully launched
               call timerProcess%terminate()
               exit
            end if
         end if
      end do

   end subroutine launchRemoteRunner

   ! No matter what, we don't want to rerun this test, so
   ! we need to increment numSkip here.
   subroutine startTest(this, testName)
      class (RobustRunner), intent(inout) :: this
      character(len=*), intent(in) :: testName
      
      this%numSkip = this%numSkip + 1

   end subroutine startTest

   subroutine endTest(this, testName)
      class (RobustRunner), intent(inout) :: this
      character(len=*), intent(in) :: testName
   end subroutine endTest

   subroutine endRun(this, result)
     use AbstractTestResult_mod
     class (RobustRunner), intent(inout) :: this
     class (AbstractTestResult), intent(in) :: result
   end subroutine endRun

   subroutine addFailure(this, testName, exceptions)
      use Exception_mod
      class (RobustRunner), intent(inout) :: this
      character(len=*), intent(in) :: testName
      type (Exception), intent(in) :: exceptions(:)

   end subroutine addFailure

   subroutine addError(this, testName, exceptions)
      use Exception_mod
      class (RobustRunner), intent(inout) :: this
      character(len=*), intent(in) :: testName
      type (Exception), intent(in) :: exceptions(:)

   end subroutine addError

   function createTestResult(this) result(tstResult)
      use TestResult_mod
      class (RobustRunner), intent(inout) :: this
      type (TestResult) :: tstResult

      tstResult = newTestResult()
    end function createTestResult

end module RobustRunner_mod
