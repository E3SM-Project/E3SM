!@before
subroutine initA()
end subroutine initA

!@after
subroutine finalA()
end subroutine finalA

! First test
!@test
subroutine testMethodA()
end subroutine testMethodA

! Second test
!@test
subroutine testMethodB
end subroutine testMethodB




module WrapbeforeAfter
   use pFUnit_mod
   implicit none
   private

contains


end module WrapbeforeAfter

function beforeAfter_suite() result(suite)
   use pFUnit_mod
   use WrapbeforeAfter
   type (TestSuite) :: suite

   external testMethodA
   external testMethodB

   external initA
   external finalA

   integer, allocatable :: npes(:)

   suite = newTestSuite('beforeAfter_suite')

   call suite%addTest(newTestMethod('testMethodA', testMethodA, initA, finalA))

   call suite%addTest(newTestMethod('testMethodB', testMethodB, initA, finalA))


end function beforeAfter_suite

