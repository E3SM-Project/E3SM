!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_FixtureTestCase_mod
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
module Test_FixtureTestCase_mod
   use TestSuite_mod
   use TestResult_mod, only: TestResult, newTestResult
   implicit none
   private

   public :: suite

contains

   function suite() result(aSuite)
      use TestSuite_mod, only: newTestSuite, TestSuite
      use TestMethod_mod, only: newTestMethod
      type (TestSuite) :: aSuite

      aSuite = newTestSuite('Test_TestCase')

!#define ADD(method) call aSuite%addTest(newTestMethod(REFLECT(method)))

      call aSuite%addTest( &
           &   newTestMethod('testRunWithFixture', &
           &                  testRunWithFixture))
      call aSuite%addTest( &
           &   newTestMethod('testBrokenTestCase', &
           &                  testBrokenTestCase))
      call aSuite%addTest( &
           &   newTestMethod('testBrokenSetUpCase', &
           &                  testBrokenSetUpCase))

   end function suite

   subroutine testRunWithFixture()
      use TestCase_mod
      use FixtureTestCase_mod, only: FixtureTestCase, newFixtureTestCase
      use FixtureTestCase_mod, only: delete
      use SerialContext_mod
      use Assert_mod, only: assertEqual
      type (FixtureTestCase) :: aTest
      type (TestResult) :: aTestResult

      aTestResult = newTestResult()
      aTest = newFixtureTestCase()
      call aTest%setSurrogate()
      call aTest%run(aTestResult, newSerialContext())
      call assertEqual('setUp run tearDown', aTest%runLog)
      call delete(aTest)

   end subroutine testRunWithFixture

   subroutine testBrokenTestCase()
      use TestCase_mod
      use BrokenTestCase_mod, only: BrokenTestCase
      use Assert_mod, only: assertEqual
      use SerialContext_mod
      type (BrokenTestCase) :: test
      type (TestResult) :: aTestResult

      call test%setSurrogate()
      call test%setName('foo')
      aTestResult = newTestResult()

      call test%run(aTestResult, newSerialContext())
      call assertEqual('setUp broken run tearDown', test%runLog)
      call assertEqual(1, aTestResult%failureCount())

   end subroutine testBrokenTestCase

   subroutine testBrokenSetUpCase()
      use TestCase_mod
      use BrokenSetUpCase_mod, only: BrokenSetUpCase
      use Assert_mod, only: assertEqual
      use SerialContext_mod
      type (BrokenSetUpCase) :: test
      type (TestResult) :: aTestResult

      call test%setSurrogate()
      call test%setName('foo')
      aTestResult = newTestResult()
      call test%run(aTestResult, newSerialContext())
      call assertEqual('broken setUp', test%runLog)
      call assertEqual(1, aTestResult%failureCount())

   end subroutine testBrokenSetUpCase

end module Test_FixtureTestCase_mod

