module TestCaseB_mod
   use pfunit_mod
   implicit none

   
!@testCase(constructor=newTestCaseB)
   type, extends(ParameterizedTestCase) :: TestCaseB
      integer, allocatable :: table(:)
      real :: phi
      real :: theta
   contains
      procedure :: setUp
      procedure :: tearDown
   end type TestCaseB

!@testParameter
   type, extends(AbstractTestParameter) :: B_Parameter
      real :: phi
      real :: theta
   contains
      procedure :: toString
   end type B_Parameter

!!$   interface TestCaseB
!!$      module procedure newTestCaseB
!!$   end interface TestCaseB

contains

   function newTestCaseB(testParameter) result(newTest)
      type (TestCaseB) :: newTest
      class (B_Parameter), intent(in) :: testParameter ! driver may subclass this
      newTest%phi = testParameter%phi
      newTest%theta = testParameter%theta
   end function newTestCaseB

   subroutine setUp(this)
      class (TestCaseB), intent(inout) :: this
      allocate(this%table(10))
   end subroutine setUp

   subroutine tearDown(this)
      class (TestCaseB), intent(inout) :: this
      deallocate(this%table)
   end subroutine tearDown

   ! First test
   !@test(testParameters={[B_Parameter(0.1,0.2),B_Parameter(0.3,0.1)]})
   subroutine testA(this)
      class (TestCaseB), intent(inout) :: this
   end subroutine testA

   ! Second test
   !@test(testParameters={[B_Parameter(0.1,0.2)]})
   subroutine testB(this)
      class (TestCaseB), intent(inout) :: this
   end subroutine testB

   function toString(this) result(str)
      class (B_Parameter), intent(in) :: this
      character(:), allocatable :: str

      str = 'no message'

   end function toString

end module TestCaseB_mod



module WrapTestCaseB_mod
   use pFUnit_mod
   use TestCaseB_mod
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(TestCaseB) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use TestCaseB_mod
        class (TestCaseB), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod, testParameter) result(aTest)
      type (WrapUserTestCase) :: aTest
      character(len=*), intent(in) :: methodName
      procedure(userTestMethod) :: testMethod
      type (B_Parameter), intent(in) :: testParameter
      aTest%TestCaseB = newTestCaseB(testParameter)

      aTest%testMethodPtr => testMethod
      call aTest%setName(methodName)
      call aTest%setTestParameter(testParameter)
   end function makeCustomTest

end module WrapTestCaseB_mod

function TestCaseB_mod_suite() result(suite)
   use pFUnit_mod
   use WrapTestCaseB_mod
   use TestCaseB_mod
   type (TestSuite) :: suite

   integer, allocatable :: npes(:)

   type (B_Parameter), allocatable :: testParameters(:)
   type (B_Parameter) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('TestCaseB_mod_suite')

   testParameters = [B_Parameter(0.1,0.2),B_Parameter(0.3,0.1)]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testA', testA, testParameter))
   end do

   testParameters = [B_Parameter(0.1,0.2)]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testB', testB, testParameter))
   end do


end function TestCaseB_mod_suite

