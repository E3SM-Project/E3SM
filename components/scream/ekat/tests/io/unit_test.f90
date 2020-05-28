#include "ekat_config.f"

module unit_test_mod
  implicit none
contains

  function test_array_io() result(nerr) bind(c)
    use array_io_mod
    use iso_c_binding

#ifdef EKAT_DOUBLE_PRECISION
    integer, parameter :: c_real = c_double
#else
    integer, parameter :: c_real = c_float
#endif


    integer(kind=c_int) :: nerr
    integer ::  i, j
    real(kind=c_real), target :: a(10,3), b(10,3)
    logical :: ok

    character(kind=c_char, len=128), parameter :: &
         filename = c_char_"unit_test_f90_array_io.dat"//C_NULL_CHAR

    do j = 1,3
       do i = 1,10
          a(i,j) = 100*j + i
       end do
    end do
    nerr = 0
    ok = array_io_write(filename, c_loc(a), size(a))
    if (.not. ok) nerr = nerr + 1
    ok = array_io_file_exists(filename)
    if (.not. ok) nerr = nerr + 1
    ok = array_io_read(filename, c_loc(b), size(b))
    if (.not. ok) nerr = nerr + 1
    do j = 1,3
       do i = 1,10
          if (a(i,j) .ne. b(i,j)) nerr = nerr + 1
       end do
    end do
  end function test_array_io

end module unit_test_mod
