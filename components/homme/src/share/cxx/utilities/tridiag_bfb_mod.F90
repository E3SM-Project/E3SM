#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

! These routines are meant only for BFB testing. The solve a diagonally dominant
! tridiagonal system A x = b, with (dl,d,du) the tridiags and x = b on
! input. See scream_tridag.hpp for performant solvers.

module tridag_bfb_mod
  implicit none
  
  interface
     subroutine tridiag_diagdom_bfb_a1x1(n, dl, d, du, x) bind(c)
       use iso_c_binding, only: c_int, c_double
       integer(c_int), value, intent(in) :: n
       real(c_double), intent(inout) :: dl(n), d(n), du(n), x(n)
     end subroutine tridiag_diagdom_bfb_a1x1

     subroutine tridiag_diagdom_bfb_a1xm(n, nrhs, dl, d, du, x) bind(c)
       use iso_c_binding, only: c_int, c_double
       integer(c_int), value, intent(in) :: n, nrhs
       real(c_double), intent(inout) :: dl(n), d(n), du(n), x(nrhs,n)
     end subroutine tridiag_diagdom_bfb_a1xm
  end interface
end module tridag_bfb_mod
