!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_AssertInteger_mod
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
module Test_AssertInteger_mod
   use AssertBasic_mod
   use Assert_mod, only: assertEqual
   use Assert_mod, only: assertLessThan
   use Assert_mod, only: assertGreaterThan
   use Assert_mod, only: assertLessThanOrEqual
   use Assert_mod, only: assertGreaterThanOrEqual
!   use AssertInteger_mod
   use TestSuite_mod, only: TestSuite, newTestSuite
   use Params_mod, only: i32, i64
   implicit none
   private

   public :: suite

contains

   function suite()
      use TestSuite_mod, only: TestSuite, newTestSuite
      use TestMethod_mod, only: newTestMethod
      use Test_mod

      type (TestSuite) :: suite

      suite = newTestSuite('AssertIntegerTests')
!#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))

      call suite%addTest( &
           &   newTestMethod('testAssertEqual_equal', &
           &                  testAssertEqual_equal))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual_unequal', &
           &                  testAssertEqual_unequal))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual_unequalWithMessage', &
           &                  testAssertEqual_unequalWithMessage))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual1D1D_equal', &
           &                  testAssertEqual1D1D_equal))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual1D1D_nonconforming', &
           &                  testAssertEqual1D1D_nonconforming))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual1D1D_conforming', &
           &                  testAssertEqual1D1D_conforming))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual1D1D_unequalA', &
           &                  testAssertEqual1D1D_unequalA))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual1D1D_unequalB', &
           &                  testAssertEqual1D1D_unequalB))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual2D2D_equal', &
           &                  testAssertEqual2D2D_equal))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual2D2D_nonconforming', &
           &                  testAssertEqual2D2D_nonconforming))
      call suite%addTest( &
           &   newTestMethod('testAssertEqual2D2D_unequal', &
           &                  testAssertEqual2D2D_unequal))

      call suite%addTest( &
           &   newTestMethod('testAssertLessThan_falseA', &
           &                  testAssertLessThan_falseA))
      call suite%addTest( &
           &   newTestMethod('testAssertLessThan_falseB', &
           &                  testAssertLessThan_falseB))
      call suite%addTest( &
           &   newTestMethod('testAssertLessThan_true', &
           &                  testAssertLessThan_true))

      call suite%addTest( &
           &   newTestMethod('testAssertLessThanOrEqual_false', &
           &                  testAssertLessThanOrEqual_false))
      call suite%addTest( &
           &   newTestMethod('testAssertLessThanOrEqual_trueA', &
           &                  testAssertLessThanOrEqual_trueA))
      call suite%addTest( &
           &   newTestMethod('testAssertLessThanOrEqual_trueB', &
           &                  testAssertLessThanOrEqual_trueB))

      call suite%addTest( &
           &   newTestMethod('testAssertGreaterThan_falseA', &
           &                  testAssertGreaterThan_falseA))
      call suite%addTest( &
           &   newTestMethod('testAssertGreaterThan_falseB', &
           &                  testAssertGreaterThan_falseB))
      call suite%addTest( &
           &   newTestMethod('testAssertGreaterThan_true', &
           &                  testAssertGreaterThan_true))

      call suite%addTest( &
           &   newTestMethod('testAssertGreaterThanOrEqual_false', &
           &                  testAssertGreaterThanOrEqual_false))
      call suite%addTest( &
           &   newTestMethod('testAssertGreaterThanOrEqual_trueA', &
           &                  testAssertGreaterThanOrEqual_trueA))
      call suite%addTest( &
           &   newTestMethod('testAssertGreaterThanOrEqual_trueB', &
           &                  testAssertGreaterThanOrEqual_trueB))

   end function suite

   subroutine testAssertEqual_equal()
      call assertEqual(2,2)
   end subroutine testAssertEqual_equal

   subroutine testAssertEqual_unequal()
      call assertEqual(2,3)
      call assertExceptionRaised('expected 2 but found: 3;  difference: |1|.')
   end subroutine testAssertEqual_unequal

   subroutine testAssertEqual_unequalWithMessage()
      call assertEqual(2,3,'what?')
      call assertExceptionRaised('what? expected 2 but found: 3;  difference: |1|.')
   end subroutine testAssertEqual_unequalWithMessage

   subroutine testAssertEqual1D1D_equal()
      call assertEqual([1,2],[1,2])
   end subroutine testAssertEqual1D1D_equal

   subroutine testAssertEqual1D1D_nonconforming()
      call assertEqual([1,2],[1,2,3])
      call assertExceptionRaised('nonconforming arrays - expected shape: [2] but found shape: [3]')
   end subroutine testAssertEqual1D1D_nonconforming

   subroutine testAssertEqual1D1D_conforming()
      call assertEqual(1,[1,1,1])
   end subroutine testAssertEqual1D1D_conforming

   subroutine testAssertEqual1D1D_unequalA()
      call assertEqual([1,2,3],[1,3,3])
      call assertExceptionRaised('expected 2 but found: 3;  difference: |1|;  first difference at element [2].')
   end subroutine testAssertEqual1D1D_unequalA

   subroutine testAssertEqual1D1D_unequalB()
      call assertEqual(1, [1,2,1])
      call assertExceptionRaised('expected 1 but found: 2;  difference: |1|;  first difference at element [2].')
   end subroutine testAssertEqual1D1D_unequalB

   subroutine testAssertEqual2D2D_equal()
      integer :: array(2,3)
      array = reshape([1,2,3,4,5,6],[2,3])
      call assertEqual(array, array)
   end subroutine testAssertEqual2D2D_equal

   subroutine testAssertEqual2D2D_nonconforming()
      integer :: expected(2,3)
      integer :: found(3,5)

      expected = 1
      found = 1
      call assertEqual(expected, found)
      call assertExceptionRaised('nonconforming arrays - expected shape: [2,3] but found shape: [3,5]')

   end subroutine testAssertEqual2D2D_nonconforming

   subroutine testAssertEqual2D2D_unequal()
      integer(kind=i32) :: expected(2,3)
      integer(kind=i64) :: found(2,3)

      expected = 1
      found = 1
      found(1,2) = -1

      call assertEqual(expected, found)
      call assertExceptionRaised('expected 1 but found: -1;  difference: |2|;  first difference at element [1, 2].')

      found(1,2) = 1
      found(2,3) = -1

      call assertEqual(expected, found)
      call assertExceptionRaised('expected 1 but found: -1;  difference: |2|;  first difference at element [2, 3].')

   end subroutine testAssertEqual2D2D_unequal

   subroutine testAssertLessThan_falseA()
      call assertLessThan(1, 1)
      call assertExceptionRaised('expected 1 to be less than: 1.')
   end subroutine testAssertLessThan_falseA

   subroutine testAssertLessThan_falseB()
      call assertLessThan(2, 1)
      call assertExceptionRaised('expected 2 to be less than: 1.')
   end subroutine testAssertLessThan_falseB

   subroutine testAssertLessThan_true()
      call assertLessThan(1, 2)
   end subroutine testAssertLessThan_true
   
   subroutine testAssertLessThanOrEqual_false()
      call assertLessThanOrEqual(2, 1)
      call assertExceptionRaised('expected 2 to be less than or equal to: 1.')
   end subroutine testAssertLessThanOrEqual_false

   subroutine testAssertLessThanOrEqual_trueA()
      call assertLessThanOrEqual(2, 2)
   end subroutine testAssertLessThanOrEqual_trueA

   subroutine testAssertLessThanOrEqual_trueB()
     call assertLessThanOrEqual(1, 2)
   end subroutine testAssertLessThanOrEqual_trueB
   
   subroutine testAssertGreaterThan_falseA()
      call assertGreaterThan(1, 1)
      call assertExceptionRaised('expected 1 to be greater than: 1.')
   end subroutine testAssertGreaterThan_falseA

   subroutine testAssertGreaterThan_falseB()
      call assertGreaterThan(1, 2)
      call assertExceptionRaised('expected 1 to be greater than: 2.')
   end subroutine testAssertGreaterThan_falseB

   subroutine testAssertGreaterThan_true()
      call assertGreaterThan(2, 1)
   end subroutine testAssertGreaterThan_true
   
   subroutine testAssertGreaterThanOrEqual_false()
      call assertGreaterThanOrEqual(1, 2)
      call assertExceptionRaised('expected 1 to be greater than or equal to: 2.')
   end subroutine testAssertGreaterThanOrEqual_false

   subroutine testAssertGreaterThanOrEqual_trueA()
      call assertGreaterThanOrEqual(2, 2)
   end subroutine testAssertGreaterThanOrEqual_trueA

   subroutine testAssertGreaterThanOrEqual_trueB()
      call assertGreaterThanOrEqual(2, 1)
   end subroutine testAssertGreaterThanOrEqual_trueB
   
end module Test_AssertInteger_mod
