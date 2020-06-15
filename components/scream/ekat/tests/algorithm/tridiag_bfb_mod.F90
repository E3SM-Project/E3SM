#include "ekat_config.f"
#ifdef EKAT_DOUBLE_PRECISION
# define c_real c_double
#else
# define c_real c_float
#endif

! These routines are meant only for BFB testing. The solve a diagonally dominant
! tridiagonal system A x = b, with (dl,d,du) the tridiags and x = b on
! input. See scream_tridag.hpp for performant solvers.

module scream_tridag_bfb_mod
  implicit none
  
  interface
     subroutine tridiag_diagdom_bfb_a1x1(n, dl, d, du, x) bind(c)
       use iso_c_binding, only: c_int, c_real
       integer(c_int), value, intent(in) :: n
       real(c_real), intent(inout) :: dl(n), d(n), du(n), x(n)
     end subroutine tridiag_diagdom_bfb_a1x1

     subroutine tridiag_diagdom_bfb_a1xm(n, nrhs, dl, d, du, x) bind(c)
       use iso_c_binding, only: c_int, c_real
       integer(c_int), value, intent(in) :: n, nrhs
       real(c_real), intent(inout) :: dl(n), d(n), du(n), x(nrhs,n)
     end subroutine tridiag_diagdom_bfb_a1xm
  end interface
end module scream_tridag_bfb_mod
