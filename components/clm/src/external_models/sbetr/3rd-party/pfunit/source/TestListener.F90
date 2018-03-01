!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: TestListener
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
module TestListener_mod
   implicit none
   private

   public :: TestListener
   public :: ListenerPointer

   type, abstract :: TestListener
      private
      logical :: useDebug = .false.
   contains
     procedure(addFailure), deferred :: addFailure
     procedure(startTest), deferred :: startTest
     procedure(endTest), deferred :: endTest
!     procedure(startRun), deferred :: startRun  ! make deferred when ready
     procedure(endRun), deferred :: endRun    ! make deferred when ready
     procedure :: addError
     procedure :: setDebug
     procedure :: debug
   end type TestListener

   type ListenerPointer
     class (TestListener), pointer :: pListener
   end type ListenerPointer

   abstract interface
      subroutine addFailure(this, testName, exceptions)
         use Exception_mod
         import TestListener
         class (TestListener), intent(inout) :: this
         character(len=*), intent(in) :: testName
         type (Exception), intent(in) :: exceptions(:)
      end subroutine addFailure

      subroutine startTest(this, testName)
         import TestListener
         class (TestListener), intent(inout) :: this
         character(len=*), intent(in) :: testName
      end subroutine startTest
    
      subroutine endTest(this, testName)
         import TestListener
         class (TestListener), intent(inout) :: this
         character(len=*), intent(in) :: testName
      end subroutine endTest

!      ! Stub for future implementation.
!      subroutine startRun(this)
!         import TestListener
!         class (TestListener), intent(inout) :: this
!      end subroutine startRun
!
      ! Stub for future implementation.
      subroutine endRun(this, result)
         use AbstractTestResult_mod, only : AbstractTestResult
         import TestListener
         class (TestListener), intent(inout) :: this
         class (AbstractTestResult), intent(in) :: result
      end subroutine endRun

   end interface

contains

   ! Most scenarios in Fortran cannot diagnose true errors, so
   ! an empty stub is provided here for convenience.
   subroutine addError(this, testName, exceptions)
      use Exception_mod, only: Exception
      class (TestListener), intent(inout) :: this
      character(len=*), intent(in) :: testName
      type (Exception), intent(in) :: exceptions(:)
   end subroutine addError

   ! Promoted from BaseTestRunner.F90. Every listener can have debug
   ! behaviors.
    subroutine setDebug(this)
       class (TestListener), intent(inout) :: this
       this%useDebug = .true.
    end subroutine setDebug

    logical function debug(this)
       class (TestListener), intent(inout) :: this
       debug = this%useDebug
    end function debug

 end module TestListener_mod
