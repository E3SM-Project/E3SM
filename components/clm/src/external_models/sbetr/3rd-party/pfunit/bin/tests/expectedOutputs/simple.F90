! First test
!@test
subroutine testMethodA()
end subroutine testMethodA

! Second test
!@test
subroutine testMethodB
end subroutine testMethodB

! An MPI test
!@mpitest(npes=[1,3,5])
subroutine testMethodC(this)
   use pfunit_mod
   class (MpiTestMethod), intent(inout) :: this
end subroutine testMethodC



module Wrapsimple
   use pFUnit_mod
   implicit none
   private

contains


end module Wrapsimple

function simple_suite() result(suite)
   use pFUnit_mod
   use Wrapsimple
   type (TestSuite) :: suite

   external testMethodA
   external testMethodB
   external testMethodC


   integer, allocatable :: npes(:)

   suite = newTestSuite('simple_suite')

   call suite%addTest(newTestMethod('testMethodA', testMethodA))

   call suite%addTest(newTestMethod('testMethodB', testMethodB))

   call suite%addTest(newMpiTestMethod('testMethodC', testMethodC, 1))
   call suite%addTest(newMpiTestMethod('testMethodC', testMethodC, 3))
   call suite%addTest(newMpiTestMethod('testMethodC', testMethodC, 5))


end function simple_suite

