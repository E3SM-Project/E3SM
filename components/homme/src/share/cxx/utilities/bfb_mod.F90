#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

! These routines are meant only for BFB testing.

module bfb_mod
  implicit none
  
  interface
     subroutine czeroulpn(a_len, a, nbit) bind(c)
       use iso_c_binding, only: c_int, c_double
       integer(c_int), value, intent(in) :: a_len, nbit
       real(c_double), intent(inout) :: a(a_len)
     end subroutine czeroulpn

     subroutine cbfb_pow(a_len, a, e) bind(c)
       use iso_c_binding, only: c_int, c_double
       integer(c_int), value, intent(in) :: a_len
       real(c_double), value, intent(in) :: e
       real(c_double), intent(inout) :: a(a_len)
     end subroutine cbfb_pow

     ! These solve a diagonally dominant tridiagonal system A x = b, with
     ! (dl,d,du) the tridiags and x = b on input. See scream_tridag.hpp for
     ! performant solvers.
     ! The last argument is used to establish whether the inputs are float or double

     subroutine tridiag_diagdom_bfb_a1x1(n, dl, d, du, x, real_size) bind(c)
       use iso_c_binding, only: c_int, c_double
       integer(c_int), value, intent(in) :: n, real_size
       real(c_double), intent(inout) :: dl(n-1), d(n), du(n-2), x(n)
     end subroutine tridiag_diagdom_bfb_a1x1
  end interface
end module bfb_mod
