!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_MpiException_mod
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
module Test_MpiException_mod
   use Test_mod
   use TestCase_mod
   use Exception_mod
   use ParallelException_mod
   use MpiTestMethod_mod
   implicit none

   private
   
   public :: suite

contains

   function suite()
      use TestSuite_mod, only: TestSuite
      use TestSuite_mod, only: newTestSuite
      use TestMethod_mod, only: newTestMethod!, TestMethod

      type (TestSuite) :: suite

      suite = newTestSuite('ExceptionTests')

!#define ADD(method, npes) call suite%addTest(newMpiTestMethod(REFLECT(method), numProcesses=npes))
      call suite%addTest( &
           &   newMpiTestMethod('test_anyExceptions_none', &
           &                  test_anyExceptions_none,  &
           &                  numProcesses=3))
      call suite%addTest( &
           &   newMpiTestMethod('test_getNumExceptions', &
           &                  test_getNumExceptions,  &
           &                  numProcesses=4))
      call suite%addTest( &
           &   newMpiTestMethod('test_gather', &
           &                  test_gather,  &
           &                  numProcesses=3))

   end function suite

   subroutine test_anyExceptions_none(this)
      use Assert_mod
      use ParallelContext_mod
      class (MpiTestMethod), intent(inout) :: this
      class (ParallelContext), allocatable :: context

      allocate(context, source=this%getContext())
      call assertFalse(anyExceptions(context))
      
      if (this%getProcessRank() == 1) then
         call throw('some message')
      end if

      call assertTrue(anyExceptions(context))

      ! clear thrown exception
      if (this%getProcessRank() == 1) then
         call assertTrue(catch('some message'))
      end if


   end subroutine test_anyExceptions_none

   subroutine test_getNumExceptions(this)
      use Assert_mod
      use ParallelContext_mod
      class (MpiTestMethod), intent(inout) :: this

      call assertEqual(0, getNumExceptions(this%getContext()))

      select case (this%getProcessRank()) 
      case (0,2)
         call throw('some message')
      end select

      call assertEqual(2, getNumExceptions(this%getContext()))

      ! clear thrown exception
      select case (this%getProcessRank()) 
      case (0,2)
         call assertTrue(catch('some message'))
      end select

   end subroutine test_getNumExceptions

   subroutine test_gather(this)
      use Assert_mod
      use ParallelContext_mod
      class (MpiTestMethod), intent(inout) :: this

      select case (this%getProcessRank()) 
      case (0)
         call throw('exception 1')
      case (1)
         call throw('exception 2')
         call throw('exception 3')
      case (2)
         call throw('exception 4')
      end select

      call gather(this%getContext())

      select case (this%getProcessRank())
      case (0)
         ! remote exceptions now local with added suffix
         call assertTrue(catch('exception 1 (PE=0)'))
         call assertTrue(catch('exception 2 (PE=1)'))
         call assertTrue(catch('exception 3 (PE=1)'))
         call assertTrue(catch('exception 4 (PE=2)'))
      case (1:)
         ! local exceptions gone
         call assertEqual(0, getNumExceptions())
      end select

   end subroutine test_gather

end module Test_MpiException_mod
