!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_TestMethod_mod
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
module Test_TestMethod_mod
   use TestSuite_mod, only: TestSuite, newTestSuite
   implicit none
   private

   public :: suite

contains

   function suite()
      use TestSuite_mod, only: TestSuite, newTestSuite
      use TestMethod_mod, only: newTestMethod
      type (TestSuite) :: suite

      suite = newTestSuite('Test_TestMethod')
      call suite%addTest(newTestMethod('testMethodWasRun', testMethodWasRun))

   end function suite

   subroutine testMethodWasRun()
      use TestCase_mod
      use TestResult_mod, only: TestResult, newTestResult
      use TestMethod_mod, only: TestMethod, newTestMethod
      use Assert_mod, only: assertEqual
      use SerialContext_mod
      type (TestMethod) :: method
      type (TestResult) :: aResult

      method = newTestMethod(name = 'testWasRun', method = testWasRun)
      aResult = newTestResult()
      call method%run(aResult, newSerialContext())
      call assertEqual(1, aResult%runCount())
      call assertEqual(1, aResult%failureCount())

   end subroutine testMethodWasRun

   subroutine testWasRun()
      use Exception_mod, only: throw
      call throw('wasRun')
   end subroutine testWasRun

end module Test_TestMethod_mod

