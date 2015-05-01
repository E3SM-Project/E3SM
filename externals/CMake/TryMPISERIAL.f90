! Test for the mct mpiserial stubs library
! in which mpi_integer == mpi_real4 and the select will
! fail to compile. (PGI doesnt complain so add something 
! it will complain about
!
program mpiserial_test
 implicit none
 include 'mpif.h' 
 integer :: i
 select case(i)
 case(mpi_integer)
 case(mpi_real4)
 end select
end program mpiserial_test
