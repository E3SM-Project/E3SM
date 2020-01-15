#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

! These routines are meant only for BFB testing.

module bfb_mod
  use kinds, only: real_kind

  implicit none
  
  interface
     subroutine czeroulpn(a_len, a, nbit, replace) bind(c)
       use iso_c_binding, only: c_int, c_double
       integer(c_int), value, intent(in) :: a_len, nbit
       real(c_double), intent(inout) :: a(a_len)
       real(c_double), intent(inout) :: replace(a_len)
     end subroutine czeroulpn

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

  interface bfb_pow
    module procedure bfb_pow_0d
    module procedure bfb_pow_1d
    module procedure bfb_pow_2d
    module procedure bfb_pow_3d
  end interface bfb_pow

  public :: bfb_pow
  integer, public, parameter :: nmax = 10

contains

  recursive function int_pow (a,k) result(y)
    real (kind=real_kind), intent(in) :: a
    integer, intent(in) :: k

    real (kind=real_kind) :: y

    if (k<=0) then
      print *, "ERROR! Cannot handle non-positive exponents"
    endif

    if (k==1) then
      y = a
    else if (mod(k,2) .eq. 0) then
      y = int_pow(a*a,k/2)
    else
      y = a*int_pow(a,k-1)
    endif
  end function int_pow

  recursive function two_a (a) result(y)
    real (kind=real_kind), intent(in) :: a

    real (kind=real_kind) :: y

    y = 1 + a*(1 + (a-1)/2*(1 + (a-2)/3*(1 + (a-3)/4)));

  end function two_a

  recursive function bfb_pow_0d (val, e) result(res)

    real (kind=real_kind), intent(in)  :: val
    real (kind=real_kind), intent(in)  :: e

    real (kind=real_kind) :: res
#ifdef CUDA_BUILD
    real (kind=real_kind) :: x,x0,a,a0,tmp
    integer :: i,n,k

    if (val<0) then
      print *, "Negative base not allowed in bfb_pow!"
    endif

    ! Taylor formula [ (1+x)^a = sum (a choose k) x^k ] holds for -1<x<1, so 0<(1+x)<2
    ! For better convergence, we ask that -0.5<x<0.5, so 0.5<(1+x)<1.5, compute
    ! the factor (2^e)^k separately, and then combine the two
    ! Note: all of this is tailored (or, if you wish, taylored...eheh) for 0<e<2

    if (val .eq. 0) then
      if (e .eq. 0) then
        print *, "Cannot do 0^0, sorry."
      endif
      res = 0.0
    else if (val .ge. 2.0) then
      tmp = val
      do while (tmp .ge. 16.0)
        k = k+4
        tmp = tmp/16.0
      enddo
      do while (tmp .ge. 1.5)
        k = k+1
        tmp = tmp/2.0
      enddo

      res = int_pow(two_a(e),k)*bfb_pow_0d(tmp,e)
    else if (val .le. 0.5) then
      tmp = val
      do while(tmp .le. 0.0625)
        k = k+4
        tmp = tmp*16.0
      enddo
      do while (tmp .le. 0.5)
        k = k+1
        tmp = tmp*2.0
      enddo

      ! Taylor formula [ (1+x)^a = sum (a choose k) x^k ] holds for -1<x<1, so 0<(1+x)<2
      res = bfb_pow_0d(tmp,e)/int_pow(two_a(e),k)
    else

      x0 = val-1.0
      a0 = e

      ! For our uses of bfb_pow (0.2<e<1.5), 10 terms should yield a max error of about 1e-3
      res=0
      x = 1.0
      a = 1.0
      do n=0,nmax
        res = res + a*x
        x = x*x0
        a = a*((a0-n)/(n+1))
      enddo
    endif
#else
    res = val**e
#endif
  end function bfb_pow_0d

  function bfb_pow_1d (arr, e) result(res)
    real (kind=real_kind), intent(in) :: e
    real (kind=real_kind), intent(in), dimension(:) :: arr

    real (kind=real_kind), dimension(lbound(arr,dim=1):ubound(arr,dim=1)) :: res
    integer :: i,n

    n = SIZE(arr)

    do i=1,n
      res(i) = bfb_pow_0d(arr(i),e)
    enddo
  end function bfb_pow_1d

  function bfb_pow_2d (arr, e) result(res)
    real (kind=real_kind), intent(in) , dimension(:,:) :: arr
    real (kind=real_kind), intent(in) :: e

    real (kind=real_kind), dimension(lbound(arr,dim=1):ubound(arr,dim=1),&
                                     lbound(arr,dim=2):ubound(arr,dim=2)) :: res
    integer :: n

    n = PRODUCT(SHAPE(arr))
    res = RESHAPE(bfb_pow_1d(RESHAPE(arr,[n]),e),SHAPE(arr))
  end function bfb_pow_2d

  function bfb_pow_3d (arr, e) result(res)
    real (kind=real_kind), intent(in) , dimension(:,:,:) :: arr
    real (kind=real_kind), intent(in) :: e

    real (kind=real_kind), dimension(lbound(arr,dim=1):ubound(arr,dim=1),&
                                     lbound(arr,dim=2):ubound(arr,dim=2),&
                                     lbound(arr,dim=3):ubound(arr,dim=3)) :: res
    integer :: n

    n = PRODUCT(SHAPE(arr))
    res = RESHAPE(bfb_pow_1d(RESHAPE(arr,[n]),e),SHAPE(arr))
  end function bfb_pow_3d
end module bfb_mod
