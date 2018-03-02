program main
   use pFUnit_mod, only: initialize
   use pFUnit_mod, only: finalize
   use pFUnit_mod, only: TestResult
   use pFUnit_mod, only: ListenerPointer
   use pFUnit_mod, only: newResultPrinter
   use pFUnit_mod, only: ResultPrinter
   implicit none

   logical :: success

   call initialize()
   success = runTests()
   call finalize(success)

contains

   logical function runTests() result(success)
      use pFUnit_mod, only: newTestSuite
      use pFUnit_mod, only: TestSuite
      use pFUnit_mod, only: TestRunner, newTestRunner
#ifdef USE_MPI
      use MpiContext_mod
      use ParallelException_mod
#else
      use SerialContext_mod
#endif
      use ParallelContext_mod

      use Test_StringConversionUtilities_mod, only: StringConversionUtilitiesSuite => suite    ! (1)
#ifdef BUILD_ROBUST
      use Test_UnixProcess_mod, only: unixProcessSuite => suite                ! (1)
#endif
      use Test_Exception_mod, only: exceptionSuite => suite                ! (2)
      use Test_AssertBasic_mod, only: assertBasicSuite => suite            !
      use Test_Assert_mod, only: assertSuite => suite                      ! (3)
      use Test_AssertInteger_mod, only: assertIntegerSuite => suite        !

      use Test_AssertReal_mod, only: assertRealSuite => suite              ! (5)
      use Test_AssertComplex_mod, only: assertComplexSuite => suite              ! (5)

      use Test_TestResult_mod, only: testResultSuite => suite              ! (6)
      use Test_TestSuite_mod, only: testSuiteSuite => suite                ! (7)

      use Test_TestMethod_mod, only: testSimpleMethodSuite => suite  ! (8)
      use Test_SimpleTestCase_mod, only: testSimpleSuite => suite          ! (9)
      use Test_FixtureTestCase_mod, only: testFixtureSuite => suite        ! (10)

      use Test_BasicOpenMP_mod, only: testBasicOpenMpSuite => suite  ! (8)

      use Test_MockCall_mod, only: testMockCallSuite => suite      ! (11)
      use Test_MockRepository_mod, only: testMockRepositorySuite => suite      ! (11)
      use Test_XmlPrinter_mod, only: testXmlPrinterSuite => suite

#ifdef BUILD_ROBUST
      use Test_RobustRunner_mod, only: testRobustRunnerSuite => suite
#endif

#ifdef USE_MPI
      use Test_MpiContext_mod, only: MpiContextSuite => suite
      use Test_MpiException_mod, only: ParallelExceptionSuite => suite
      use Test_MpiTestCase_mod, only: MpiTestCaseSuite => suite
      use Test_MpiParameterizedTestCase_mod, only: MpiParameterizedTestCaseSuite => suite
#endif
      use iso_fortran_env, only: OUTPUT_UNIT

      type (TestSuite) :: allTests
      type (TestRunner) :: runner
      type (TestResult) :: tstResult

#ifdef INTEL_13
      type (ResultPrinter), target :: printer
#endif

!- MLR 2014-0224-1209 Why would intel 13 not like "listeners"?
      type (ListenerPointer), allocatable :: l1(:)
      allocate(l1(1))
#ifndef INTEL_13
      allocate(l1(1)%pListener, source=newResultPrinter(OUTPUT_UNIT))
#else
      printer = newResultPrinter(OUTPUT_UNIT)
      l1(1)%pListener => printer
#endif

      allTests = newTestSuite('allTests')
      runner = newTestRunner(l1)

#define ADD(suite) call allTests%addTest(suite())

      ADD(StringConversionUtilitiesSuite)
#ifdef BUILD_ROBUST
      ADD(UnixProcessSuite)
#endif
      ADD(exceptionSuite)

      ADD(assertBasicSuite)
      ADD(assertSuite)
      ADD(assertIntegerSuite)

      ADD(assertRealSuite)
      ADD(assertComplexSuite)

      ADD(testResultSuite)
      ADD(testSuiteSuite)

      ADD(testSimpleMethodSuite)
      ADD(testSimpleSuite)
      ADD(testFixtureSuite)

      ADD(testBasicOpenMpSuite)

      ADD(testMockCallSuite)
      ADD(testMockRepositorySuite)

      ADD(testXmlPrinterSuite)

#ifdef BUILD_ROBUST
      ADD(testRobustRunnerSuite)
#endif

#ifdef USE_MPI
      ADD(MpiContextSuite)
      ADD(ParallelExceptionSuite)
      ADD(MpiTestCaseSuite)
      ADD(MpiParameterizedTestCaseSuite)
#endif

#ifdef USE_MPI
      tstResult = runner%run(allTests, newMpiContext())
#else
      tstResult = runner%run(allTests, newSerialContext())
#endif
      success = tstResult%wasSuccessful()

  end function runTests

end program main

!if ! defined INTEL_13 
!...
!else
!      class (ListenerPointer), allocatable :: l1(:)
!      type (ResultPrinter), target :: printer
!      allocate(listeners1(1))
!      printer = newResultPrinter(OUTPUT_UNIT)
!      listeners1(1)%pListener=>printer
!endif


