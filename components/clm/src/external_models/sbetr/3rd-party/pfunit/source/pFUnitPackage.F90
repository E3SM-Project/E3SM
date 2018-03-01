!-------------------------------------------------------------------------------
! NASA/GSFC Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: pFUnit
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
!
! This module packages pFUnit entities while simultaneously inserting
! a prefix on all names.  Some developers may provide this explicit
! naming convention.  Others may choose to use the vanilla pFUnit module that
! has no such prefix.   
!
! The default prefix is "pf_", but just edit the #define below to suit
! your own preference.
!
!

#ifndef PFUNIT_PREFIX
#define PFUNIT_PREFIX pf_
#endif

#define TOKEN(a) a
#define RENAME(item) TOKEN(PFUNIT_PREFIX)TOKEN(item) => item

module pFUnit

   use pFUnit_mod, only: RENAME(initialize)
   use pFUnit_mod, only: RENAME(finalize)

   use SourceLocation_mod, only: RENAME(SourceLocation)
   use Exception_mod, only:  RENAME(throw), RENAME(catch), RENAME(catchNext)
   use Exception_mod, only:  RENAME(anyExceptions)
   use ParallelException_mod, only:  RENAME(anyExceptions)
   use Assert_mod, only: RENAME(assertFail)
   use Assert_mod, only: RENAME(assertTrue), RENAME(assertFalse)
   use Assert_mod, only: RENAME(assertSameShape)
   use Assert_mod, only: RENAME(assertEqual)
   use Assert_mod, only: RENAME(assertAny), RENAME(assertAll)
   use Assert_mod, only: RENAME(assertNone), RENAME(assertNotAll)
   use Assert_mod, only: RENAME(assertLessThan)
   use Assert_mod, only: RENAME(assertLessThanOrEqual)
   use Assert_mod, only: RENAME(assertGreaterThan)
   use Assert_mod, only: RENAME(assertGreaterThanOrEqual)
   use Assert_mod, only: RENAME(assertExceptionRaised)

   use Assert_mod, only: RENAME(assertIsNan)
   use Assert_mod, only: RENAME(assertIsFinite)

   ! workaround for ifort 13.0
   use Test_mod, only: Test

   use Test_mod, only: RENAME(Test)
   use TestCase_mod, only: RENAME(TestCase)
   use TestSuite_mod, only: RENAME(TestSuite)
   use TestSuite_mod, only: RENAME(newTestSuite)
   use TestMethod_mod, only: RENAME(TestMethod)
   use TestMethod_mod, only: RENAME(newTestMethod)
   use TestResult_mod, only: RENAME(TestResult)
   use BaseTestRunner_mod, only: RENAME(BaseTestRunner)
   use TestRunner_mod, only: RENAME(TestRunner)
   use TestRunner_mod, only: RENAME(newTestRunner)
#ifdef BUILD_ROBUST
   use RobustRunner_mod, only: RENAME(RobustRunner)
#endif

   use TestListener_mod, only: RENAME(ListenerPointer)
   use XmlPrinter_mod, only: RENAME(XmlPrinter)
   use DebugListener_mod, only: RENAME(DebugListener)

   use ParallelContext_mod, only: RENAME(ParallelContext)
   use SerialContext_mod, only: RENAME(SerialContext)
   use SerialContext_mod, only: RENAME(newSerialContext)
#ifdef USE_MPI
   use MpiContext_mod, only: RENAME(MpiContext)
   use MpiContext_mod, only: RENAME(newMpiContext)
   use MpiTestCase_mod, only: RENAME(MpiTestCase)
   use MpiTestParameter_mod, only: RENAME(MpiTestParameter)
   use MpiTestMethod_mod, only: RENAME(MpiTestMethod)
   use MpiTestMethod_mod, only: RENAME(newMpiTestMethod)
#endif

   use AbstractTestParameter_mod, only: RENAME(AbstractTestParameter)
   use ParameterizedTestCase_mod, only: RENAME(ParameterizedTestCase)

   implicit none
   public ! Nothing private in this module, just renaming exports.

   ! workaround for ifort 13.0
   private :: Test
contains

    function run() result(a)

      integer :: a
      a = 0

   end function run


end module pFUnit
