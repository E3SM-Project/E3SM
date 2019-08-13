module marsaglia
! public domain code
! made available from http://www.fortran.com/
! downloaded by pjr on 03/16/04
! converted to vector form, functions inlined by pjr,mvr on 05/10/2004

use shr_kind_mod,    only: r8 => shr_kind_r8

implicit none
save
private
public :: kissvec

contains
  subroutine kissvec(seed1,seed2,seed3,seed4,ran_arr)
! The  KISS (Keep It Simple Stupid) random number generator. Combines:
! (1) The congruential generator x(n)=69069*x(n-1)+1327217885, period 2^32.
! (2) A 3-shift shift-register generator, period 2^32-1,
! (3) Two 16-bit multiply-with-carry generators, period 597273182964842497>2^59
!  Overall period>2^123; 
!
    implicit none
    REAL(r8), DIMENSION (:), INTENT(INOUT)  :: ran_arr
    INTEGER, DIMENSION(:), INTENT(INOUT) :: seed1,seed2,seed3,seed4
    integer :: i,sz,kiss
    integer :: m, k, n

! inline function 
    m(k, n) = ieor (k, ishft (k, n) )

    sz = SIZE(ran_arr)
    DO i = 1, sz
       seed1(i) = 69069 * seed1(i) + 1327217885
       seed2(i) = m (m (m (seed2(i), 13), - 17), 5)
       seed3(i) = 18000 * iand (seed3(i), 65535) + ishft (seed3(i), - 16)
       seed4(i) = 30903 * iand (seed4(i), 65535) + ishft (seed4(i), - 16)
       kiss = seed1(i) + seed2(i) + ishft (seed3(i), 16) + seed4(i)
       ran_arr(i) = kiss*2.328306e-10_r8 + 0.5_r8
    end do
    
  end subroutine kissvec

end module marsaglia
