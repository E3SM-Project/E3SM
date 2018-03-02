module TestCaseC_mod
   use pfunit_mod
   implicit none

   
!@testCase(constructor=newTestCaseC)
   type, extends(MpiTestCase) :: TestCaseC
      integer, allocatable :: table(:)
      real :: phi
      real :: theta
   contains
      procedure :: setUp
      procedure :: tearDown
   end type TestCaseC

!@testParameter(constructor=newC_Parameter)
   type, extends(MpiTestParameter) :: C_Parameter
      real :: phi
      real :: theta
   contains
      procedure :: toString
   end type C_Parameter

   interface newC_Parameter
      module procedure newC_Parameter_phiTheta
      module procedure newC_Parameter_case
   end interface newC_Parameter

contains

   function newTestCaseC(testParameter) result(newTest)
      type (TestCaseC) :: newTest
      class (C_Parameter), intent(in) :: testParameter ! driver may subclass this
      newTest%phi = testParameter%phi
      newTest%theta = testParameter%theta
   end function newTestCaseC

   subroutine setUp(this)
      class (TestCaseC), intent(inout) :: this
      allocate(this%table(10))
   end subroutine setUp

   subroutine tearDown(this)
      class (TestCaseC), intent(inout) :: this
      deallocate(this%table)
   end subroutine tearDown

   ! First test
   !@test(npes=[1,3],testParameters={[newC_Parameter(1,0.1,0.1)]})
   subroutine testA(this)
      class (TestCaseC), intent(inout) :: this
   end subroutine testA

   ! Second test
   !@test(testParameters={paramGenerator()})
   subroutine testB(this)
      class (TestCaseC), intent(inout) :: this
   end subroutine testB

   ! Third test
   !@test(cases=[1,2])
   subroutine testC(this)
      class (TestCaseC), intent(inout) :: this
   end subroutine testC

   function newC_Parameter_phiTheta(npes, phi, theta) result(testParameter)
      type (C_Parameter) :: testParameter
      integer, intent(in) :: npes

      real, intent(in) :: phi
      real, intent(in) :: theta

      testParameter%phi = phi
      testParameter%theta = theta

      call testParameter%setNumProcessesRequested(npes)

   end function newC_Parameter_phiTheta

   elemental function newC_Parameter_case(i) result(testParameter)
      type (C_Parameter) :: testParameter
      integer, intent(in) :: i

      integer :: npes

      select case(i)
      case (1)
         testParameter%phi = 0.1
         testParameter%theta = 0.2
         npes = 3
      case (2)
         testParameter%phi = -0.1
         testParameter%theta = 0.2
         npes = 4
      end select

      call testParameter%setNumProcessesRequested(npes)

   end function newC_Parameter_case

   function paramGenerator() result(testParameters)
      type (C_Parameter), allocatable :: testParameters(:)

      testParameters = [ &
           & newC_Parameter(npes=1,phi=0.1,theta=0.2), &
           & newC_Parameter(npes=2,phi=0.1,theta=0.2), &
           & newC_Parameter(npes=3,phi=0.2,theta=0.2), &
           & newC_Parameter(npes=4,phi=0.1,theta=0.3) ]

   end function paramGenerator

   function toString(this) result(str)
      class (C_Parameter), intent(in) :: this
      character(:), allocatable :: str

      str = 'no message'

   end function toString

end module TestCaseC_mod



module WrapTestCaseC_mod
   use pFUnit_mod
   use TestCaseC_mod
   implicit none
   private

   public :: WrapUserTestCase
   public :: makeCustomTest
   type, extends(TestCaseC) :: WrapUserTestCase
      procedure(userTestMethod), nopass, pointer :: testMethodPtr
   contains
      procedure :: runMethod
   end type WrapUserTestCase

   abstract interface
     subroutine userTestMethod(this)
        use TestCaseC_mod
        class (TestCaseC), intent(inout) :: this
     end subroutine userTestMethod
   end interface

contains

   subroutine runMethod(this)
      class (WrapUserTestCase), intent(inout) :: this

      call this%testMethodPtr(this)
   end subroutine runMethod

   function makeCustomTest(methodName, testMethod, testParameter, npesRequested) result(aTest)
      type (WrapUserTestCase) :: aTest
      character(len=*), intent(in) :: methodName
      procedure(userTestMethod) :: testMethod
      type (C_Parameter), intent(in) :: testParameter
      integer, optional, intent(in) :: npesRequested

      aTest%TestCaseC = newTestCaseC(testParameter)

      aTest%testMethodPtr => testMethod
      call aTest%setName(methodName)
      call aTest%setTestParameter(testParameter)
     if (present(npesRequested)) then
         call aTest%setNumProcessesRequested(npesRequested) 
     end if

   end function makeCustomTest

end module WrapTestCaseC_mod

function TestCaseC_mod_suite() result(suite)
   use pFUnit_mod
   use WrapTestCaseC_mod
   use TestCaseC_mod
   type (TestSuite) :: suite

   integer, allocatable :: npes(:)

   type (C_Parameter), allocatable :: testParameters(:)
   type (C_Parameter) :: testParameter
   integer :: iParam 
   integer, allocatable :: cases(:) 
 
   suite = newTestSuite('TestCaseC_mod_suite')

   testParameters = [newC_Parameter(1,0.1,0.1)]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testA', testA, testParameter, npesRequested=1))
   end do
   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testA', testA, testParameter, npesRequested=3))
   end do

   testParameters = paramGenerator()

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testB', testB, testParameter))
   end do

   cases = [1,2]
   testParameters = [(newC_Parameter(cases(iCase)), iCase = 1, size(cases))]

   do iParam = 1, size(testParameters)
      testParameter = testParameters(iParam)
   call suite%addTest(makeCustomTest('testC', testC, testParameter))
   end do


end function TestCaseC_mod_suite

