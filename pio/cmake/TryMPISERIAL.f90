! Test for the mct mpiserial stubs library
! in which mpi_integer == mpi_real4 and the select will
! fail to compile.
!
program mpiserial_test
 include 'mpif.h' 

 select case(i)
 case(mpi_integer)
 case(mpi_real4)
 end select
end program mpiserial_test
