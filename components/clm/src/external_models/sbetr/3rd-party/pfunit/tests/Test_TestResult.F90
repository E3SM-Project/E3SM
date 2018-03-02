!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_TestResult_mod
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
module Test_TestResult_mod
   use TestSuite_mod, only: TestSuite, newTestSuite
   use TestResult_mod, only: TestResult, newTestResult
   use TestResult_mod, only: newTestResult, TestResult
   use TestCase_mod
   use SimpleTestCase_mod, only: newSimpleTestCase, SimpleTestCase
   implicit none
   private

   public :: suite

contains

   function suite()
      use TestSuite_mod, only: TestSuite, newTestSuite
      use TestResult_mod, only: TestResult, newTestResult
      use TestCase_mod
      use TestMethod_mod, only: newTestMethod
      type (TestSuite) :: suite

      suite = newTestSuite('TestResultSuite')

!#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))


      call suite%addTest( &
           &   newTestMethod('testGetNumRun', &
           &                  testGetNumRun))
      call suite%addTest( &
           &   newTestMethod('testGetNumFailed', &
           &                  testGetNumFailed))

      call suite%addTest( &
           &   newTestMethod('testAddListenerStart', &
           &                  testAddListenerStart))
      call suite%addTest( &
           &   newTestMethod('testAddListenerEnd', &
           &                  testAddListenerEnd))
      call suite%addTest( &
           &   newTestMethod('testAddListenerFailure', &
           &                  testAddListenerFailure))

   end function suite

   subroutine testGetNumRun()
      use Assert_mod, only: assertEqual
      use TestResult_mod, only: newTestResult, TestResult
!!$      use TestCase_mod
      type (TestResult) :: aResult
!!$      class(TestCase), pointer :: tstCase

      aResult = newTestResult()
!!$      call assertEqual(0, aResult%runCount())
!!$
!!$      tstCase => newSimpleTestCase(method1,'method1')
!!$
!!$      call aResult%startTest(tstCase%getSurrogate())
!!$      call aResult%endTest(tstCase%getSurrogate())
!!$      call assertEqual(1, aResult%runCount())
!!$
!!$      call aResult%startTest(tstCase%getSurrogate())
!!$      call aResult%endTest(tstCase%getSurrogate())
!!$      call assertEqual(2, aResult%runCount())

   end subroutine testGetNumRun

   subroutine testGetNumFailed()
      use Assert_mod, only: assertEqual
      use Exception_mod, only: newException
      use SimpleTestCase_mod, only: SimpleTestCase
      use SurrogateTestCase_mod
      use TestCase_mod

      type (TestResult) :: aResult
      
      type (SimpleTestCase) :: aTest
      call aTest%setSurrogate()
      aResult = newTestResult()
      call assertEqual(0, aResult%failureCount())

      call aResult%addFailure(aTest%getSurrogate(), [newException('fail')])
      call assertEqual(1, aResult%failureCount())

      call aResult%addFailure(aTest%getSurrogate(), [newException('fail again')])
      call assertEqual(2, aResult%failureCount())

   end subroutine testGetNumFailed

   subroutine testAddListenerEnd()
      use TestListener_mod
      use MockListener_mod
      use Assert_mod
      use SimpleTestCase_mod
      use SurrogateTestCase_mod
      use TestCase_mod

      type (TestResult) :: result
      type (MockListener) :: listener
      type (SimpleTestCase) :: tstCase
      
      result = newTestResult()
      call result%addListener(listener)
      tstCase = newSimpleTestCase('method1', method1)
      call result%endTest(tstCase%getSurrogate())
      call assertEqual('endTest() was called', listener%log)

   end subroutine testAddListenerEnd

   subroutine testAddListenerStart()
      use TestListener_mod
      use MockListener_mod
      use Assert_mod
      use SimpleTestCase_mod
      type (TestResult) :: result
      type (MockListener) :: listener
      
      type (SimpleTestCase) :: tstCase

      result = newTestResult()
      call result%addListener(listener)
      tstCase = newSimpleTestCase('method1', method1)
      call result%startTest(tstCase%getSurrogate())
      call assertEqual('startTest() was called', trim(listener%log))

   end subroutine testAddListenerStart

   subroutine testAddListenerFailure()
      use TestListener_mod
      use MockListener_mod
      use Assert_mod
      use Exception_mod
      use SimpleTestCase_mod
      use SurrogateTestCase_mod
      use TestCase_mod
      
      type (TestResult) :: result
      type (MockListener) :: listener
      type (Exception) :: anException
      
      class(TestCase), allocatable :: tstCase
      
      result = newTestResult()
      call result%addListener(listener)
      allocate(tstCase, source = newSimpleTestCase('method1', method1))
      call result%addFailure(tstCase%getSurrogate(), [anException])
      call assertEqual('addFailure() was called', listener%log)

   end subroutine testAddListenerFailure

end module Test_TestResult_mod
