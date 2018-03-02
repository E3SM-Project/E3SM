module TestCaseA_mod
   use pfunit_mod
   implicit none

   
!@testCase
   type, extends(TestCase) :: TestCaseA
      integer :: componentI
   contains
      procedure :: setUp
      procedure :: tearDown
   end type TestCaseA

contains

   subroutine setUp(this)
      class (TestCaseA), intent(inout) :: this
      this%componentI = 5
   end subroutine setUp

   subroutine tearDown(this)
      class (TestCaseA), intent(inout) :: this
      this%componentI = -1
   end subroutine tearDown

   ! First test
   !@test
   subroutine testA(this)
      class (TestCaseA), intent(inout) :: this
   end subroutine testA

   ! Second test
   !@test
   subroutine testB(this)
      class (TestCaseA), intent(inout) :: this
   end subroutine testB

end module TestCaseA_mod

module WrapTestCaseA_mod
   use pFUnit_mod
   use TestCaseA_mod
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(TestCaseA) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use TestCaseA_mod
        class (TestCaseA), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod) result(aTest)
      type (WrapUserTestCase) :: aTest
      character(len=*), intent(in) :: methodName
      procedure(userTestMethod) :: testMethod
      aTest%testMethodPtr => testMethod
      call aTest%setName(methodName)
   end function makeCustomTest

end module WrapTestCaseA_mod

function TestCaseA_mod_suite() result(suite)
   use pFUnit_mod
   use WrapTestCaseA_mod
   use TestCaseA_mod
   type (TestSuite) :: suite

   integer, allocatable :: npes(:)

   suite = newTestSuite('TestCaseA_mod_suite')

   call suite%addTest(makeCustomTest('testA', testA))

   call suite%addTest(makeCustomTest('testB', testB))


end function TestCaseA_mod_suite

