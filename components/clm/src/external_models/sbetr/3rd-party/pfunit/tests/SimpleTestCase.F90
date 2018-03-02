!#include "reflection.h"
!-------------------------------------------------------------------------------
! NASA/GSFC, Advanced Software Technology Group
!-------------------------------------------------------------------------------
!  MODULE: SimpleTestCase_mod
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
module SimpleTestCase_mod
   use TestCase_mod, only: TestCase
   implicit none
   private

   public :: suite
   public :: newSimpleTestCase
   public :: SimpleTestCase
   public :: method1, method2
   public :: methodWith2Exceptions

   type, extends(TestCase) :: SimpleTestCase
      character(len=20), public :: runLog
      procedure(method), pointer :: testMethod => null()
   contains
      procedure :: runMethod
   end type SimpleTestCase

   abstract interface
      subroutine method(this)
        use Test_mod
        import SimpleTestCase
        class (SimpleTestCase), intent(inOut) :: this
      end subroutine method
   end interface

contains

   function suite()
     use TestSuite_mod, only: TestSuite, newTestSuite
      type (TestSuite) :: suite

      suite = newTestSuite('SimpleTestCase')

!#define ADD(method) call suite%addTest(newSimpleTestCase(REFLECT(method)))

      call suite%addTest( &
           &   newSimpleTestCase('method1', &
           &                      method1))
      call suite%addTest( &
           &   newSimpleTestCase('method2', &
           &                      method2))
      call suite%addTest( &
           &   newSimpleTestCase('methodWith2Exceptions', &
           &                      methodWith2Exceptions))
      
   end function suite

   function newSimpleTestCase(name, userMethod) result(this)
      type(SimpleTestCase) :: this
      character(len=*), intent(in) :: name
      procedure(method) :: userMethod

      this%testMethod => userMethod
      call this%setName(name)

    end function newSimpleTestCase

   recursive subroutine runMethod(this)
      class(SimpleTestCase), intent(inOut) :: this
      call this%testMethod()
   end subroutine runMethod

   subroutine method1(this)
      class (SimpleTestCase), intent(inOut) :: this
      this%runLog = 'run method1'
   end subroutine method1

   subroutine method2(this)
      class (SimpleTestCase), intent(inOut) :: this
      this%runLog = 'run method2'
   end subroutine method2

   subroutine methodWith2Exceptions(this)
      use Exception_mod, only: throw
      class (SimpleTestCase), intent(inOut) :: this

      call throw('failure A')
      call throw('failure B')
   end subroutine methodWith2Exceptions

   subroutine delete_(this)
      type (SimpleTestCase), intent(inOut) :: this
!!$      nullify(this%testMethod)
   end subroutine delete_

end module SimpleTestCase_mod

