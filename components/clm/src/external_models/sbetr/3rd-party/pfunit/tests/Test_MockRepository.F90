!#include "reflection.h"
#define HERE print*,__LINE__,__FILE__

!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: Test_MockRepository_mod
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


module SUT_mod
   implicit none
   private

   public :: SUT

   type SUT
      integer :: intComponent
   contains
      procedure :: method1
   end type SUT

contains

   subroutine method1(this)
      class (SUT), intent(in) :: this
   end subroutine method1

end module SUT_mod

module MockSUT_mod
   use MockRepository_mod
   use SUT_mod
   implicit none
   private

   public :: MockSUT
   public :: newMockSUT

   type, extends(SUT) :: MockSUT
      class (MockRepository), pointer :: mocker => null()
   contains
      procedure :: method1
      procedure :: verifyMocking
   end type MOCKSUT

contains

   function newMockSUT(repository) result(mock)
      type (MockSUT), allocatable :: mock
      class (MockRepository), target :: repository

      allocate(mock)
      mock%mocker => repository

   end function newMockSUT

!TODO - make FINAL routine once gfortran supports it
   subroutine verifyMocking(this)
      use Exception_mod
      class (MockSUT), intent(inout) :: this

      if (associated(this%mocker)) then
         call this%mocker%verifyMocking(this)
      end if

   end subroutine verifyMocking

   subroutine method1(this)
      class (MockSUT), intent(in) :: this
      call this%mocker%hasCalled(this, 'method1')
   end subroutine method1

end module MockSUT_mod

module Test_MockRepository_mod
   use TestSuite_mod
   use MockRepository_mod
   use Exception_mod
   use Assert_mod

   use SUT_mod
   use MockSUT_mod
   implicit none
   private

   public :: suite

   ! test that all registered objects are finalized (checked)

contains

!#define ADD(method) call suite%addTest(newTestMethod(REFLECT(method)))

   function suite()
      use TestSuite_mod, only: newTestSuite, TestSuite
      use TestMethod_mod, only: newTestMethod
      type (TestSuite) :: suite

      suite = newTestSuite('Test_MockRepository')

      call suite%addTest( &
           &   newTestMethod('testNoAction', &
           &                  testNoAction))
      call suite%addTest( &
           &   newTestMethod('testExpectMethod_NotCalled', &
           &                  testExpectMethod_NotCalled))
      call suite%addTest( &
           &   newTestMethod('testExpectMethod_IsCalled', &
           &                  testExpectMethod_IsCalled))
      call suite%addTest( &
           &   newTestMethod('testExpectMethod_CalledDifferentMethod', &
           &                  testExpectMethod_CalledDifferentMethod))

   end function suite

   subroutine testNoAction()
      class (MockRepository), pointer :: mocker
      type (SUT) :: object
      type (MockSUT) :: mockObject
      
      mocker => newMockRepository()
      mockObject = newMockSUT(mocker)
      call mocker%delete() ! move to a final

   end subroutine testNoAction

   subroutine testExpectMethod_NotCalled()

      call internalProcedure() ! verification is when object is final-ized
      call assertTrue(catch('Expected method not called: method1() on object of class MockSUT.'))

   contains

      subroutine internalProcedure()
         class (MockRepository), pointer :: mocker
         type (SUT) :: object
         type (MockSUT) :: mockObject

         mocker => newMockRepository()
         mockObject = newMockSUT(mocker)
         call mocker%expectCall(mockObject,'method1')
         call mockObject%verifyMocking()

      end subroutine internalProcedure

   end subroutine testExpectMethod_NotCalled

   subroutine testExpectMethod_IsCalled()

      call internalProcedure() ! verification is when object is final-ized

   contains

      subroutine internalProcedure()
         class (MockRepository), pointer :: mocker
         type (SUT) :: object
         type (MockSUT) :: mockObject

         mocker => newMockRepository()
         mockObject = newMockSUT(mocker)
         call mocker%expectCall(mockObject,'method1')
         call mockObject%method1()
         call mocker%delete() ! move to a final

      end subroutine internalProcedure

   end subroutine testExpectMethod_IsCalled

   subroutine testExpectMethod_CalledDifferentMethod()

      call internalProcedure() ! verification is when object is final-ized
      call assertTrue(catch('Expected method not called: method1() on object of class MockSUT.'), &
           & 'Failed to distinguish among method names.')

   contains

      subroutine internalProcedure()
         class (MockRepository), pointer :: mocker
         type (SUT) :: object
         type (MockSUT) :: mockObject

         mocker => newMockRepository()
         mockObject = newMockSUT(mocker)
         call mocker%expectCall(mockObject,'method2')
         call mockObject%method1()
         call mockObject%verifyMocking()

      end subroutine internalProcedure

   end subroutine testExpectMethod_CalledDifferentMethod

end module Test_MockRepository_mod
