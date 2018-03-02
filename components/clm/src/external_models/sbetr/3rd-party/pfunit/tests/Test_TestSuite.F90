!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_TestSuite_mod
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
module Test_TestSuite_mod
   use TestSuite_mod, only: newTestSuite, TestSuite
   use TestResult_mod
   implicit none
   private

   public :: suite

   character(len=80) :: log

   ! Internal mock for TestResult
   type, extends(TestResult) :: Verbose
      character(len=80) :: log
   contains
      procedure :: run
   end type Verbose


contains

   function suite()
      use TestCase_mod, only: TestCase
      use TestMethod_mod, only: newTestMethod
      use TestSuite_mod, only: newTestSuite, TestSuite
      type (TestSuite) :: suite

      suite = newTestSuite('TestSuiteSuite')

!#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))

      call suite%addTest( &
           &   newTestMethod('testCountTestCases', &
           &                  testCountTestCases))
      call suite%addTest( &
           &   newTestMethod('testCountTestCasesNestedA', &
           &                  testCountTestCasesNestedA))
      call suite%addTest( &
           &   newTestMethod('testCountTestCasesNestedB', &
           &                  testCountTestCasesNestedB))
      call suite%addTest( &
           &   newTestMethod('testCountTestCasesNestedC', &
           &                  testCountTestCasesNestedC))
      call suite%addTest( &
           &   newTestMethod('testGetTestCases', &
           &                  testGetTestCases))

   end function suite

   subroutine testCountTestCases()
      use TestSuite_mod, only: newTestSuite, TestSuite
      use SimpleTestCase_mod, only: newSimpleTestCase
      use SimpleTestCase_mod, only: method1, method2
      use TestSuite_mod, only: newTestSuite, TestSuite
      use Assert_mod, only: assertEqual
      type (TestSuite) :: suite

      suite = newTestSuite('aSuite')
      call assertEqual(0, suite%countTestCases())
      call suite%addTest(newSimpleTestCase('method1', method1))
      call assertEqual(1, suite%countTestCases())
      call suite%addTest(newSimpleTestCase('method2', method2))
      call assertEqual(2, suite%countTestCases())

   end subroutine testCountTestCases

   subroutine testCountTestCasesNestedA()
      use TestSuite_mod, only: newTestSuite, TestSuite
      use Assert_mod, only: assertEqual

      type (TestSuite) :: innerSuite
      type (TestSuite) :: outerSuite

      innerSuite = newTestSuite('inner')
      outerSuite = newTestSuite('outer')
      call outerSuite%addTest(innerSuite)
      call assertEqual(0, outerSuite%countTestCases())

   end subroutine testCountTestCasesNestedA

   subroutine testCountTestCasesNestedB()
      use TestSuite_mod, only: newTestSuite, TestSuite
      use SimpleTestCase_mod, only: SimpleTestCase
      use Assert_mod, only: assertEqual
      type (TestSuite) :: innerSuite
      type (TestSuite) :: outerSuite

      type (SimpleTestCase) :: aTest

      call aTest%setName('aTest')
      innerSuite = newTestSuite('inner')
      outerSuite = newTestSuite('outer')
      call innerSuite%addTest(aTest)
      call outerSuite%addTest(innerSuite)
      call assertEqual(1, outerSuite%countTestCases())
      call assertEqual(1, innerSuite%countTestCases())

   end subroutine testCountTestCasesNestedB

   !
   ! Complex Suite nested structure:
   !   topSuite
   !      ->  suiteA
   !          -> Test1
   !          -> suiteC
   !             -> Test1
   !             -> Test2
   !      ->  suiteB
   !          -> Test1
   !          -> Test2
   !
   subroutine testCountTestCasesNestedC()
      use TestSuite_mod, only: newTestSuite, TestSuite
      use SimpleTestCase_mod, only: SimpleTestCase
      use Assert_mod, only: assertEqual
      type (TestSuite) :: suiteA, suiteB, suiteC, topSuite
      type (SimpleTestCase) :: aTest

      call aTest%setName('aTest')
      topSuite = newTestSuite('top')
      suiteA = newTestSuite('A')
      suiteB = newTestSuite('B')
      suiteC = newTestSuite('C')

      call suiteC%addTest(aTest)
      call suiteC%addTest(aTest)

      call suiteB%addTest(aTest)
      call suiteB%addTest(aTest)

      call suiteA%addTest(aTest)
      call suiteA%addTest(suiteC)

      call topSuite%addTest(suiteA)
      call topSuite%addTest(suiteB)

      call assertEqual(2+2+1, topSuite%countTestCases())

   end subroutine testCountTestCasesNestedC

   subroutine testGetTestCases()
      use Test_mod
      use TestCase_mod
      use TestMethod_mod
      use SerialContext_mod
      use Assert_mod

      type (TestSuite) :: top
      type (TestSuite) :: childA, childB
      type (Verbose) :: aResult
      type (TestCaseReference), allocatable :: testCases(:)
      integer :: i

      childA = newTestSuite('childA')
      call childA%addTest(newTestMethod('a1', myTestMethod))
      call childA%addTest(newTestMethod('a2', myTestMethod))
      call childA%addTest(newTestMethod('a3', myTestMethod))

      childB = newTestSuite('childB')
      call childB%addTest(newTestMethod('b1', myTestMethod))
      call childB%addTest(newTestMethod('b2', myTestMethod))

      top = newTestSuite('top')
      call top%addTest(childA)
      call top%addTest(childB)

      aResult%TestResult = newTestResult()
      aResult%log = ''

#if (defined(__INTEL_COMPILER) && (INTEL_13))
      testCases = top%getTestCases()
#else
      call top%getTestCases(testCases)
#endif
      do i = 1, size(testCases)
         call testCases(i)%test%run(aResult, newSerialContext())
      end do

      call assertEqual('::childA.a1::childA.a2::childA.a3::childB.b1::childB.b2', aResult%log)

   end subroutine testGetTestCases

   subroutine myTestMethod()
   end subroutine myTestMethod

   recursive subroutine run(this, test, context)
      use TestCase_mod
      use SurrogateTestCase_mod
      use ParallelContext_mod
      class (Verbose), intent(inout) :: this
      class (SurrogateTestCase) :: test
      class (ParallelContext), intent(in) :: context

      this%log = trim(this%log)//'::'//trim(test%getName())

   end subroutine run

end module Test_TestSuite_mod
