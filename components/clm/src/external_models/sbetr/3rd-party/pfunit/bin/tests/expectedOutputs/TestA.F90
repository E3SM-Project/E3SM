module TestA_mod
   use pfunit_mod
   implicit none

contains

   ! First test
   !@test
   subroutine testMethodA()
   end subroutine testMethodA

   ! Second test
   !@test
   subroutine testMethodB()
   end subroutine testMethodB
   
   !@mpitest(npes=[1,3,5])
   subroutine testMethodC(this)
      class (MpiTestMethod), intent(inout) :: this
   end subroutine testMethodC

end module TestA_mod



module WrapTestA_mod
   use pFUnit_mod
   use TestA_mod
   implicit none
   private

contains


end module WrapTestA_mod

function TestA_mod_suite() result(suite)
   use pFUnit_mod
   use WrapTestA_mod
   use TestA_mod
   type (TestSuite) :: suite

   integer, allocatable :: npes(:)

   suite = newTestSuite('TestA_mod_suite')

   call suite%addTest(newTestMethod('testMethodA', testMethodA))

   call suite%addTest(newTestMethod('testMethodB', testMethodB))

   call suite%addTest(newMpiTestMethod('testMethodC', testMethodC, 1))
   call suite%addTest(newMpiTestMethod('testMethodC', testMethodC, 3))
   call suite%addTest(newMpiTestMethod('testMethodC', testMethodC, 5))


end function TestA_mod_suite

