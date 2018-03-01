!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: ResultPrinter
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
module ResultPrinter_mod
   use Exception_mod
   use TestListener_mod, only : TestListener
   implicit none
   private

   public :: ResultPrinter
   public :: newResultPrinter

   type, extends(TestListener) :: ResultPrinter
      integer :: unit
      integer :: column
   contains
      procedure :: addFailure
      procedure :: addError
      procedure :: startTest
      procedure :: endTest
      procedure :: endRun
      procedure :: print
      procedure :: printHeader
      procedure :: printFailures
      procedure :: printFooter
      procedure :: incrementColumn
   end type ResultPrinter

   integer, parameter :: MAX_COLUMN = 80
   logical, parameter :: DEBUG = .false.
!!$   logical, parameter :: DEBUG = .true.

contains

  function newResultPrinter(unit)
     type (ResultPrinter) :: newResultPrinter
     integer, intent(in) :: unit

     newResultPrinter%unit = unit
     newResultPrinter%column = 0

  end function newResultPrinter

  subroutine addFailure(this, testName, exceptions)
     use Exception_mod
     class (ResultPrinter), intent(inOut) :: this
     character(len=*), intent(in) :: testName
     type (Exception), intent(in) :: exceptions(:)

     write(this%unit,'("F")', advance='no')
     call this%incrementColumn()

  end subroutine addFailure

  subroutine addError(this, testName, exceptions)
     use Exception_mod
     class (ResultPrinter), intent(inOut) :: this
     character(len=*), intent(in) :: testName
     type (Exception), intent(in) :: exceptions(:)

     write(this%unit,'("E")', advance='no')
     call this%incrementColumn()

  end subroutine addError

  subroutine startTest(this, testName)
     class (ResultPrinter), intent(inOut) :: this
     character(len=*), intent(in) :: testName

     write(this%unit,'(".")', advance='no')
     call this%incrementColumn()

     if (DEBUG) then
        write(this%unit,*)trim(testName)
        call flush(this%unit)
     end if

   end subroutine startTest

  subroutine endTest(this, testName)
     class (ResultPrinter), intent(inOut) :: this
     character(len=*), intent(in) :: testName

     if (DEBUG) then
        write(this%unit,*)trim(testName)
        call flush(this%unit)
     end if

   end subroutine endTest

   subroutine endRun(this, result)
      use AbstractTestResult_mod, only : AbstractTestResult
      class (ResultPrinter), intent(inout) :: this
      class (AbstractTestResult), intent(in) :: result

      call this%print(result)

    end subroutine endRun

   subroutine print(this, result)
      use AbstractTestResult_mod, only : AbstractTestResult
      class (ResultPrinter), intent(in) :: this
      class (AbstractTestResult), intent(in) :: result

      call this%printHeader(result%getRunTime())
      call this%printFailures('Error', result%getErrors())
      call this%printFailures('Failure', result%getFailures())
      call this%printFooter(result)

   end subroutine print

   subroutine printHeader(this, runTime)
      class (ResultPrinter), intent(in) :: this
      real, intent(in) :: runTime

      write(this%unit,*)
      write(this%unit,'(a,1x,f12.3,1x,a)') 'Time: ', runTime, 'seconds'
      write(this%unit,*)" "

   end subroutine printHeader

   subroutine printFailures(this, label, failures)
!?      u TestResult_mod
      use TestFailure_mod
      use SourceLocation_mod
      class (ResultPrinter), intent(in) :: this
      character(len=*), intent(in) :: label
      type (TestFailure), intent(in) :: failures(:)

      type (TestFailure) :: aFailedTest
      integer :: i, j
      character(len=80) :: locationString

      do i = 1, size(failures)
         aFailedTest = failures(i)

         do j= 1, size(aFailedTest%exceptions)
            locationString = aFailedTest%exceptions(j)%location%toString()

            write(this%unit,*) label,' in: ', trim(aFailedTest%testName)
            write(this%unit,*) '  Location: ', trim(locationString)
            write(this%unit,'(a,1x,a)') aFailedTest%exceptions(j)%getMessage()
            write(this%unit,*)' '
         end do
      end do

   end subroutine printFailures

   subroutine printFooter(this, result)
      use AbstractTestResult_mod
      class (ResultPrinter), intent(in) :: this
      class (AbstractTestResult), intent(in) :: result

      if (result%wasSuccessful()) then
         write(this%unit,*)"OK"
         write(this%unit,'(a,i0,a)',advance='no')" (", result%runCount(), " test"
         if (result%runCount() > 1) then
            write(this%unit,'(a)')"s)"
         else
            write(this%unit,'(a)')")"
         end if
      else
         write(this%unit,*)"FAILURES!!!"
         write(this%unit,'(a,i0,a,i0,a,i0)')"Tests run: ", result%runCount(), &
              & ", Failures: ",result%failureCount(), &
              & ", Errors: ",result%errorCount()

      end if

   end subroutine printFooter

     subroutine incrementColumn(this)
        class (ResultPrinter), intent(inout) :: this

        this%column = this%column + 1

        if (this%column >= MAX_COLUMN) then
           write(this%unit,*) ! newline
           this%column = 0
        end if

     end subroutine incrementColumn

end module ResultPrinter_mod
