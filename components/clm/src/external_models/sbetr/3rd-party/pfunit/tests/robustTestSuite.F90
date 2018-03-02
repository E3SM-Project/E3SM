!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: robustTestSuite_mod
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune,  NASA/GSFC
!!
!! @date
!! 21 Mar 2015
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 21 Mar 2015 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------
module robustTestSuite_mod
   use pFUnit_mod
   implicit none
   private

   public :: suite

contains

   function suite()
      use TestSuite_mod, only: TestSuite, newTestSuite
      use TestMethod_mod, only: newTestMethod
      type (TestSuite) :: suite

      suite = newTestSuite('StringConversionUtilities')
!#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))

      call suite%addTest( &
           &   newTestMethod('testRunSucceeds', &
           &                  testRunSucceeds))
      call suite%addTest( &
           &   newTestMethod('testRunMultipleExceptions', &
           &                  testRunMultipleExceptions))
      call suite%addTest( &
           &   newTestMethod('testRunAssertFailure', &
           &                  testRunAssertFailure))
      call suite%addTest( &
           &   newTestMethod('testRunStops', &
           &                  testRunStops))
      call suite%addTest( &
           &   newTestMethod('testRunHangs', &
           &                  testRunHangs))

   end function suite

   subroutine testRunSucceeds()
     ! do nothing
   end subroutine testRunSucceeds

   subroutine testRunMultipleExceptions()
     use Assert_mod
     ! do nothing
     call assertTrue(1 == 2)
     call assertTrue(1 == 3)
     call assertTrue(1 == 4)
   end subroutine testRunMultipleExceptions

   subroutine testRunAssertFailure()
      use Assert_mod
      ! do nothing
      call assertTrue(1 == 2)
   end subroutine testRunAssertFailure

   subroutine testRunStops()
      call runStops()
   end subroutine testRunStops

   ! This will stop the framework cold.  The robust runner
   ! should detect this - reporting it as an Error.
   subroutine runStops()
      stop
   end subroutine runStops

   ! This test will hang.  The robust runner
   ! should detect this - reporting it as a hung process.
   subroutine testRunHangs()
      do
      end do
   end subroutine testRunHangs

end module robustTestSuite_mod
