!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: DebugListener
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
module DebugListener_mod
   use TestListener_mod
   implicit none
   private

   public :: DebugListener

   type, extends(TestListener) :: DebugListener
      integer :: unit
   contains
      procedure :: addFailure
      procedure :: startTest
      procedure :: endTest
      procedure :: endRun
   end type DebugListener

   interface DebugListener
      module procedure newDebugListener_unit
      module procedure newDebugListener_default
   end interface DebugListener

contains

   function newDebugListener_unit(unit) result(listener)
      type (DebugListener) :: listener
      integer, intent(in) :: unit
      call listener%setDebug()
      listener%unit = unit
   end function newDebugListener_unit

   function newDebugListener_default() result(listener)
      use iso_fortran_env, only: OUTPUT_UNIT
      type (DebugListener) :: listener
      call listener%setDebug()
      listener = DebugListener(OUTPUT_UNIT)
   end function newDebugListener_default

   subroutine addFailure(this, testName, exceptions)
     use Exception_mod
     class (DebugListener), intent(inOut) :: this
     character(len=*), intent(in) :: testName
     type (Exception), intent(in) :: exceptions(:)

     write(this%unit,*)'Failure in <',trim(testName),'>'
     flush(this%unit)

  end subroutine addFailure

  subroutine startTest(this, testName)
     class (DebugListener), intent(inOut) :: this
     character(len=*), intent(in) :: testName

     write(this%unit,*)new_line('A')
     write(this%unit,*)'Start: <',trim(testName),'>'
     flush(this%unit)
   end subroutine startTest

  subroutine endTest(this, testName)
     class (DebugListener), intent(inOut) :: this
     character(len=*), intent(in) :: testName

     write(this%unit,*)'  end: <',trim(testName),'>'
     flush(this%unit)

   end subroutine endTest

   subroutine endRun(this, result)
     use AbstractTestResult_mod, only : AbstractTestResult
     class (DebugListener), intent(inout) :: this
     class (AbstractTestResult), intent(in) :: result
   end subroutine endRun

end module DebugListener_mod
