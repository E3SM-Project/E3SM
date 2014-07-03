#include <shr_isnan.h>

module shr_infnan_mod

!!     Inf_NaN_Detection module 
!!     Copyright(c) 2003, Lahey Computer Systems, Inc.
!!     Copies of this source code, or standalone compiled files 
!!     derived from this source may not be sold without permission
!!     from Lahey Computers Systems. All or part of this module may be 
!!     freely incorporated into executable programs which are offered
!!     for sale. Otherwise, distribution of all or part of this file is
!!     permitted, provided this copyright notice and header are included.

!!     This module exposes four elemental functions:
!!
!!     isnan(x)    - test for a "not a number" value
!!
!!     isinf(x)    - test for either a positive or negative "infinite" value
!!
!!     isposinf(x) - test for a positive "infinite" value
!!
!!     isneginf(x) - test for a negative "infinite" value
!!
!!     Each function accepts a single or double precision real argument, and
!!     returns a true or false value to indicate the presence of the value 
!!     being tested for. If the argument is array valued, the function returns
!!     a conformable logical array, suitable for use with the ANY function, or
!!     as a logical mask.
!!
!!     Each function operates by transferring the bit pattern from a real 
!!     variable to an integer container. Unless testing for + or - infinity,
!!     the sign bit is cleared to zero. The value is exclusive ORed with
!!     the value being tested for. The integer result of the IEOR function is
!!     converted to a logical result by comparing it to zero.
!!
    use shr_kind_mod, only : r8  => SHR_KIND_R8, r4 => SHR_KIND_R4
    implicit none

    private

    public :: shr_infnan_isnan
    public :: shr_infnan_isinf
    public :: shr_infnan_isposinf
    public :: shr_infnan_isneginf

    ! Kind numbers for single and double precision integer containers
    integer, parameter :: Single = selected_int_kind(precision(1.0_r4))
    integer, parameter :: Double = selected_int_kind(precision(1.0_r8))

    ! Single precision IEEE values
    integer(Single), parameter :: sNaN    = Z"7FC00000"
    integer(Single), parameter :: sPosInf = Z"7F800000"
    integer(Single), parameter :: sNegInf = Z"FF800000"

    ! Double precision IEEE values
    integer(Double), parameter :: dNaN    = Z"7FF8000000000000"
    integer(Double), parameter :: dPosInf = Z"7FF0000000000000"
    integer(Double), parameter :: dNegInf = Z"FFF0000000000000"

    ! Locatation of single and double precision sign bit (Intel)
    ! Subtract one because bit numbering starts at zero
    integer, parameter :: SPSB = bit_size(sNaN) - 1
    integer, parameter :: DPSB = bit_size(dNaN) - 1

   interface shr_infnan_isnan
#ifdef NOFTN_INTRINSIC
      module procedure c_sisnan_scalar
      module procedure c_sisnan_1D
      module procedure c_sisnan_2D
      module procedure c_sisnan_3D
      module procedure c_sisnan_4D
      module procedure c_sisnan_5D
      module procedure c_sisnan_6D
      module procedure c_sisnan_7D
      module procedure c_disnan_scalar
      module procedure c_disnan_1D
      module procedure c_disnan_2D
      module procedure c_disnan_3D
      module procedure c_disnan_4D
      module procedure c_disnan_5D
      module procedure c_disnan_6D
      module procedure c_disnan_7D
#else
      module procedure sisnan
      module procedure disnan
#endif
   end interface   

   interface shr_infnan_isinf
      module procedure sisinf
      module procedure disinf
   end interface   
   
   interface shr_infnan_isposinf
      module procedure sisposinf
      module procedure disposinf
   end interface   
   
   interface shr_infnan_isneginf
      module procedure sisneginf
      module procedure disneginf
   end interface   


   integer :: shr_sisnan
   external :: shr_sisnan
   integer :: shr_disnan
   external :: shr_disnan

contains    

!
! If FORTRAN intrinsic's exist use them
!
#ifndef NOFTN_INTRINSIC

  ! Single precision test for NaN
  elemental function sisnan(x) result(res)
#ifdef SunOS
    use IEEE_ARITHMETIC, only : IEEE_IS_NAN
#endif
    implicit none
    real(r4), intent(in) :: x
    logical :: res
#ifdef AIX
    intrinsic :: IEEE_IS_NAN
#elif defined(ISNAN_INTRINSIC)
    intrinsic :: isnan
#endif

#if defined(AIX) || defined(SunOS)
    res = IEEE_IS_NAN(x)
#elif defined(ISNAN_INTRINSIC)
    res = isnan(x)
#endif

  end function  

  ! Double precision test for NaN
  elemental function disnan(d) result(res)
#ifdef SunOS
    use IEEE_ARITHMETIC, only : IEEE_IS_NAN
#endif
    implicit none
    real(r8), intent(in) :: d
    logical :: res
#ifdef AIX
    intrinsic :: IEEE_IS_NAN
#elif defined(ISNAN_INTRINSIC)
    intrinsic :: isnan
#endif

#if defined(AIX) || defined(SunOS)
    res = IEEE_IS_NAN(d)
#elif defined(ISNAN_INTRINSIC)
    res = isnan(d)
#endif

  end function  

!
! Otherwise link to a C function call that either uses the C90 isnan function or a x != x check
! To do this we can NOT use an elemental function, because the call to C is NOT pure.
! So we have to create functions for all the precision and array ranks
!
#else

  function c_sisnan_scalar(x) result(res)
    real(r4), intent(in) :: x
    logical :: res

    res = (shr_sisnan(x) /= 0)
  end function c_sisnan_scalar

  function c_sisnan_1D(x) result(res)
    real(r4), intent(in) :: x(:)
    logical :: res(size(x))

    integer :: i 

    do i = 1, size(x)
       res(i) = (shr_sisnan(x(i)) /= 0)
    end do
  end function c_sisnan_1D
  
  function c_sisnan_2D(x) result(res)
    real(r4), intent(in) :: x(:,:)
    logical :: res(size(x,1),size(x,2))

    integer :: i, j

    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j) = (shr_sisnan(x(i,j)) /= 0)
    end do
    end do
  end function c_sisnan_2D
  
  function c_sisnan_3D(x) result(res)
    real(r4), intent(in) :: x(:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3))

    integer :: i, j, k

    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k) = (shr_sisnan(x(i,j,k)) /= 0)
    end do
    end do
    end do
  end function c_sisnan_3D
  
  function c_sisnan_4D(x) result(res)
    real(r4), intent(in) :: x(:,:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3),size(x,4))

    integer :: i, j, k, m

    do m = 1, size(x,4)
    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k,m) = (shr_sisnan(x(i,j,k,m)) /= 0)
    end do
    end do
    end do
    end do
  end function c_sisnan_4D
  
  function c_sisnan_5D(x) result(res)
    real(r4), intent(in) :: x(:,:,:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5))

    integer :: i, j, k, m, n

    do n = 1, size(x,5)
    do m = 1, size(x,4)
    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k,m,n) = (shr_sisnan(x(i,j,k,m,n)) /= 0)
    end do
    end do
    end do
    end do
    end do
  end function c_sisnan_5D
  
  function c_sisnan_6D(x) result(res)
    real(r4), intent(in) :: x(:,:,:,:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6))

    integer :: i, j, k, m, n, o

    do o = 1, size(x,6)
    do n = 1, size(x,5)
    do m = 1, size(x,4)
    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k,m,n,o) = (shr_sisnan(x(i,j,k,m,n,o)) /= 0)
    end do
    end do
    end do
    end do
    end do
    end do
  end function c_sisnan_6D
  
  function c_sisnan_7D(x) result(res)
    real(r4), intent(in) :: x(:,:,:,:,:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6),size(x,7))

    integer :: i, j, k, m, n, o, p

    do p = 1, size(x,7)
    do o = 1, size(x,6)
    do n = 1, size(x,5)
    do m = 1, size(x,4)
    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k,m,n,o,p) = (shr_sisnan(x(i,j,k,m,n,o,p)) /= 0)
    end do
    end do
    end do
    end do
    end do
    end do
    end do
  end function c_sisnan_7D
  
  function c_disnan_scalar(x) result(res)
    real(r8), intent(in) :: x
    logical :: res

    res = (shr_disnan(x) /= 0)
  end function c_disnan_scalar

  function c_disnan_1D(x) result(res)
    real(r8), intent(in) :: x(:)
    logical :: res(size(x))

    integer :: i 

    do i = 1, size(x)
       res(i) = (shr_disnan(x(i)) /= 0)
    end do
  end function c_disnan_1D
  
  function c_disnan_2D(x) result(res)
    real(r8), intent(in) :: x(:,:)
    logical :: res(size(x,1),size(x,2))

    integer :: i, j

    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j) = (shr_disnan(x(i,j)) /= 0)
    end do
    end do
  end function c_disnan_2D
  
  function c_disnan_3D(x) result(res)
    real(r8), intent(in) :: x(:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3))

    integer :: i, j, k

    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k) = (shr_disnan(x(i,j,k)) /= 0)
    end do
    end do
    end do
  end function c_disnan_3D
  
  function c_disnan_4D(x) result(res)
    real(r8), intent(in) :: x(:,:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3),size(x,4))

    integer :: i, j, k, m

    do m = 1, size(x,4)
    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k,m) = (shr_disnan(x(i,j,k,m)) /= 0)
    end do
    end do
    end do
    end do
  end function c_disnan_4D
  
  function c_disnan_5D(x) result(res)
    real(r8), intent(in) :: x(:,:,:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5))

    integer :: i, j, k, m, n

    do n = 1, size(x,5)
    do m = 1, size(x,4)
    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k,m,n) = (shr_disnan(x(i,j,k,m,n)) /= 0)
    end do
    end do
    end do
    end do
    end do
  end function c_disnan_5D
  
  function c_disnan_6D(x) result(res)
    real(r8), intent(in) :: x(:,:,:,:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6))

    integer :: i, j, k, m, n, o

    do o = 1, size(x,6)
    do n = 1, size(x,5)
    do m = 1, size(x,4)
    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k,m,n,o) = (shr_disnan(x(i,j,k,m,n,o)) /= 0)
    end do
    end do
    end do
    end do
    end do
    end do
  end function c_disnan_6D
  
  function c_disnan_7D(x) result(res)
    real(r8), intent(in) :: x(:,:,:,:,:,:,:)
    logical :: res(size(x,1),size(x,2),size(x,3),size(x,4),size(x,5),size(x,6),size(x,7))

    integer :: i, j, k, m, n, o, p

    do p = 1, size(x,7)
    do o = 1, size(x,6)
    do n = 1, size(x,5)
    do m = 1, size(x,4)
    do k = 1, size(x,3)
    do j = 1, size(x,2)
    do i = 1, size(x,1)
       res(i,j,k,m,n,o,p) = (shr_disnan(x(i,j,k,m,n,o,p)) /= 0)
    end do
    end do
    end do
    end do
    end do
    end do
    end do
  end function c_disnan_7D

#endif
  
  ! Single precision test for Inf
  elemental function sisinf(x) result(res)
    real(r4), intent(in) :: x
    logical :: res
    res = ieor(ibclr(transfer(x,sPosInf),SPSB), sPosInf) == 0
  end function  

  ! Double precision test for Inf
  elemental function disinf(d) result(res)
    real(r8), intent(in) :: d
    logical :: res
    res = ieor(ibclr(transfer(d,dPosInf),DPSB), dPosInf) == 0
  end function  
  
  ! Single precision test for +Inf
  elemental function sisposinf(x) result(res)
    real(r4), intent(in) :: x
    logical :: res
    res = ieor(transfer(x,sPosInf), sPosInf) == 0
  end function  

  ! Double precision test for +Inf
  elemental function disposinf(d) result(res)
    real(r8), intent(in) :: d
    logical :: res
    res = ieor(transfer(d,dPosInf), dPosInf) == 0
  end function  
  
  ! Single precision test for -Inf
  elemental function sisneginf(x) result(res)
    real(r4), intent(in) :: x
    logical :: res
    res = ieor(transfer(x,sNegInf), sNegInf) == 0
  end function  

  ! Double precision test for -Inf
  elemental function disneginf(d) result(res)
    real(r8), intent(in) :: d
    logical :: res
    res = ieor(transfer(d,dNegInf), dNegInf) == 0
  end function  

end module shr_infnan_mod


