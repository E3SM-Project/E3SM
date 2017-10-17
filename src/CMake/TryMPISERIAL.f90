! Test for the mct mpiserial library
! in which mpi_cart is not defined (and likely never will be)
! Failure to compile means mpi-serial is being used.
!
program mpiserial_test
 implicit none
 include 'mpif.h'
 integer :: i
 select case(i)
 case(mpi_cart)
 end select
end program mpiserial_test
