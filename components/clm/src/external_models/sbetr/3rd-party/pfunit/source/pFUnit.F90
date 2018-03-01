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
module pFUnit_mod
   use SourceLocation_mod
   use Exception_mod
   use ParallelException_mod
   use Expectation_mod
   use Test_mod
   use TestSuite_mod
   use TestCase_mod
   use TestMethod_mod
   use AbstractTestParameter_mod
   use ParameterizedTestCase_mod
   use TestResult_mod
   use TestRunner_mod
   use BaseTestRunner_mod
   use SubsetRunner_mod

   use TestListener_mod
   use XmlPrinter_mod
   use ResultPrinter_mod
   use DebugListener_mod

#ifdef BUILD_ROBUST
   use RobustRunner_mod
#endif
   use Assert_mod
!  AssertReal mod
   use ParallelContext_mod
   use SerialContext_mod
#ifdef USE_MPI
   use MpiContext_mod
   use MpiTestCase_mod
   use MpiTestParameter_mod
   use MpiTestMethod_mod
#endif
   implicit none
   private

   public :: initialize
   public :: finalize

   public :: SourceLocation
   public :: Test
   public :: TestSuite, newTestSuite
   public :: TestMethod, newTestMethod
   public :: TestResult
   public :: TestRunner, newTestRunner
   public :: BaseTestRunner
   public :: SubsetRunner

   public :: ListenerPointer
   public :: ResultPrinter
   public :: newResultPrinter
   public :: newXmlPrinter
   public :: DebugListener

#ifdef BUILD_ROBUST
   public :: RobustRunner
#endif
   public :: TestCase
   public :: AbstractTestParameter
   public :: ParameterizedTestCase
   public :: ParallelContext
   public :: SerialContext, newSerialContext
#ifdef USE_MPI
   public :: MpiContext, newMpiContext
   public :: MpiTestCase
   public :: MpiTestParameter
   public :: MpiTestMethod, newMpiTestMethod
#endif

   public :: assertFail
   public :: assertTrue, assertFalse
   public :: assertEqual
   public :: assertAny
   public :: assertAll
   public :: assertNone
   public :: assertNotAll
   public :: assertLessThan, assertLessThanOrEqual
   public :: assertGreaterThan, assertGreaterThanOrEqual
   public :: assertRelativelyEqual
   public :: assertExceptionRaised
   public :: assertSameShape
   public :: assertIsNan
   public :: assertIsFinite

   public :: throw, catchNext, catch, anyExceptions

   public :: Expectation, Subject, Predicate
   public :: wasCalled, wasNotCalled, wasCalledOnce

   ! Optional arguments for assertEqual
   public :: WhitespaceOptions
   public :: IGNORE_ALL, TRIM_ALL, KEEP_ALL, IGNORE_DIFFERENCES


#ifdef USE_MPI
   logical :: useMpi_
#endif

contains

   subroutine initialize(useMpi)
      logical, optional, intent(in) :: useMpi
#ifdef USE_MPI
      include 'mpif.h'
      integer :: error

      useMpi_ = .true.
      if (present(useMpi)) useMpi_ = useMpi

      if (useMpi_) then
         call mpi_init(error)
      end if
#endif
      call initializeGlobalExceptionList()

   end subroutine initialize

   subroutine finalize(successful)
#ifdef NAG
      use f90_unix_proc, only: exit
#endif
      logical, intent(in) :: successful

      logical :: allSuccessful
      logical :: amRoot

#ifdef USE_MPI
      integer :: error
      integer :: rank
      include 'mpif.h'

      allSuccessful = successful
      if (useMpi_) then
         call MPI_Comm_rank(MPI_COMM_WORLD, rank, error)
         amRoot = (rank == 0)
         call MPI_Bcast(allSuccessful, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, error)
         call mpi_finalize(error)
      else
         ! If using MPI-PFUNIT on serial code, ensure amRoot is set.
         amRoot = .true.
      end if
#else
      amRoot = .true.
      allSuccessful = successful
#endif

   if (.not. allSuccessful) then
#if defined(NAG) || defined(PGI)
      call exit(-1)
#else

      if (amRoot) then
         error stop '*** Encountered 1 or more failures/errors during testing. ***'
      else
         error stop
      end if
#endif
   end if

   end subroutine finalize

end module pFUnit_mod
