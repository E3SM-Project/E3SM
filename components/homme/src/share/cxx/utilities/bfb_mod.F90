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

  interface
    function cxx_log(input) bind(C)
      use iso_c_binding, only: c_double

      !arguments:
      real(kind=c_double), value, intent(in) :: input

      ! return
      real(kind=c_double)            :: cxx_log
    end function cxx_log
  end interface

  interface bfb_pow
    module procedure bfb_pow_0d
    module procedure bfb_pow_1d
    module procedure bfb_pow_2d
    module procedure bfb_pow_3d
  end interface bfb_pow

  public :: bfb_pow, int_pow, power_of_two

contains

  recursive function int_pow (a,k) result(y)
    real (kind=real_kind), intent(in) :: a
    integer, intent(in) :: k

    real (kind=real_kind) :: y,x
    integer :: i,e,pow2
    integer, parameter :: imax = 30

    if (k<0) then
      print *, "ERROR! Cannot handle non-positive exponents"
    endif

    y = 1
    if (k /= 0) then
      e = k
      x = a
      pow2=1
      do i=0,imax
        if (mod(e,pow2*2) /= 0) then
          y = y * x
          e = e-pow2
        endif
        pow2 = pow2*2
        x = x*x
      enddo
    endif
  end function int_pow

  function power_of_two (a) result(y)
    real (kind=real_kind), intent(in) :: a

    ! 2^e approx with n-th order Taylor expansion of 2^x
    ! around x=1, where derivs are 2*[ln(2)]^k
    real (kind=real_kind) :: y,l2a
    real(kind=real_kind), parameter :: l2 = 0.693147180559945D0
    integer, parameter :: order = 10
    integer :: n

    l2a = l2*a
    ! Note: if order = 0, y=1, which is the initial value.
    y = 1.0
    do n=order,1,-1
      y = y * (l2a/n)
      y = y + 1.0D0
    enddo

  end function power_of_two

  recursive function bfb_pow_0d (val, a) result(y)

    real (kind=real_kind), intent(in)  :: val
    real (kind=real_kind), intent(in)  :: a

    real (kind=real_kind) :: y
#ifdef HOMMEXX_ENABLE_GPU_F90
    real (kind=real_kind) :: x,tmp,factor,e
    integer :: i,n,k
    integer, parameter :: order = 5

    if (val<0) then
      print *, "Negative base not allowed in bfb_pow!"
    endif

    ! Taylor formula [ (1+x)^a = sum (a choose k) x^k ] holds for -1<x<1, so 0<(1+x)<2
    ! For better convergence, we ask that -0.5<x<0.5, so 0.5<(1+x)<1.5, compute
    ! the factor (2^e)^k separately, and then combine the two
    ! Note: all of this is tailored (or, if you wish, taylored...eheh) for 0<e<2

    if (val .eq. 0) then
      if (a .eq. 0) then
        print *, "Cannot do 0^0, sorry."
      endif
      y = 0.0
    else
      factor = 1.0
      ! Taylor formula [ (1+x)^a = sum (a choose k) x^k ] holds for -1<x<1, so 0<(1+x)<2
      ! To converge faster, we enforce 0.5<1+x<1.5. We set x = 2^k * tmp, with tmp in
      ! the [0.5,1.5] bound. Then raise the two separately to the power e.
      tmp = val
      e = a

      if (e < 0) then
        e = -e
        tmp = 1.0D0 / tmp
      endif

      if (tmp .ge. 1.5) then
        k = 0
        do while (tmp .ge. 16.0)
          k = k+4
          tmp = tmp/16.0
        enddo
        do while (tmp .ge. 1.5)
          k = k+1
          tmp = tmp/2.0
        enddo

        factor = int_pow(power_of_two(e),k)
      else if (tmp .le. 0.5) then
        k = 0
        do while(tmp .le. 0.0625)
          k = k+4
          tmp = tmp*16.0
        enddo
        do while (tmp .le. 0.5)
          k = k+1
          tmp = tmp*2.0
        enddo

        factor = 1.0/int_pow(power_of_two(e),k)
      endif

      ! Note: if order = 0, y=1, which is the initial value.
      x = tmp - 1
      y = 1.0D0
      do n=order,1,-1
        y = y * (( (e-(n-1))/n ) * x)
        y = y + 1.0D0
      enddo

      y = factor*y
    endif
#else
    y = val**a
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
