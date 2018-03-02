!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_Assert_mod
!
!> @brief
!! <BriefDescription>
!!
!! @author
!! Tom Clune,  NASA/GSFC
!!
!! @date
!! 20 Mar 2015
!! 
!! @note <A note here.>
!! <Or starting here...>
!
! REVISION HISTORY:
!
! 20 Mar 2015 - Added the prologue for the compliance with Doxygen. 
!
!-------------------------------------------------------------------------------
module Test_Assert_mod
   use TestSuite_mod
   use Assert_mod
   use Exception_mod, only: NULL_MESSAGE
   use Exception_mod, only: catch
   use Exception_mod, only: getNumExceptions
   implicit none
   private

   public :: suite 

contains

   function suite() result(aSuite)
      use Test_mod
      use TestMethod_mod
      use TestSuite_mod
      type (TestSuite) :: aSuite

      aSuite = newTestSuite('Assert')

!#define ADD(method) call aSuite%addTest(newTestMethod(REFLECT(method)))

      call aSuite%addTest( &
           &   newTestMethod('testAssertEqualStringDiffer1st', &
           &                  testAssertEqualStringDiffer1st))

      call aSuite%addTest( &
           &   newTestMethod('testAssertWithLocation', &
           &                  testAssertWithLocation))
   end function suite


   subroutine testAssertEqualStringDiffer1st()
      call assertEqual(expected="a string A", found="string B")
      call assertTrue(catch('String assertion failed:' // new_line('A') // &
           & '    expected: <"a string A">' // new_line('A') // &
           & '   but found: <"string B">' // new_line('A') // &
           & '  first diff:   ^'))
   end subroutine testAssertEqualStringDiffer1st

   subroutine testAssertWithLocation
      use SourceLocation_mod
      call assertTrue(.false., 'intentional fail', &
           & SourceLocation(fileName='nowhere', lineNumber=5))
      call assertTrue(catch('intentional fail'))
   end subroutine testAssertWithLocation

end module Test_Assert_mod
