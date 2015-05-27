module MathfuncMod
#include "shr_assert.h"

!Module contains mathematical functions for some elementary manipulations
!Created by Jinyun Tang
   use shr_kind_mod, only: r8 => shr_kind_r8
   use clm_varctl      , only : iulog
   use abortutils,   only: endrun   
   use shr_log_mod ,   only: errMsg => shr_log_errMsg
   !contains subroutines for array manipulation
   !and some useful functions
implicit none
   save
   private
   public :: cumsum
   public :: swap
   public :: minmax
   public :: cumdif
   public :: diff
   public :: safe_div
   public :: dot_sum
   public :: addone
   public :: asc_sort_vec
   public :: is_bounded
   interface cumsum
     module procedure cumsum_v, cumsum_m
   end interface cumsum
   interface swap
     module procedure swap_i, swap_r, swap_rv
   end interface swap
contains
!-------------------------------------------------------------------------------   
   function heviside(x)result(ans)
   implicit none
   real(r8), intent(in) :: x
   real(r8) :: ans
   
   if(x>0._r8)then
      ans = 1._r8
   else
      ans = 0._r8
   endif
   end function heviside


!-------------------------------------------------------------------------------
   subroutine swap_i(a,b)
   implicit none
   integer, intent(inout) :: a, b
   integer :: c
   
   c = a
   a = b
   b = c
   
   end subroutine swap_i
!-------------------------------------------------------------------------------
   subroutine swap_r(a,b)
   implicit none
   real(r8), intent(inout) :: a, b
   real(r8) :: c
   
   c = a
   a = b
   b = c
   
   end subroutine swap_r
!-------------------------------------------------------------------------------
   subroutine swap_rv(a,b)
   implicit none
   real(r8), dimension(:), intent(inout) :: a, b
   real(r8), dimension(size(a)) :: c
   
   integer :: n
   
   if(size(a)/=size(b))then
      write(iulog,*)'the input vectors are not of same size in swap_rv'
      write(iulog,*)'clm model is stopping'      
      call endrun()
   endif
   
   c = a
   a = b
   b = c
   
   end subroutine swap_rv
!-------------------------------------------------------------------------------
   function minmax(x)result(ans)
   !returnd the minimum and maximum of the input vector
   implicit none
   real(r8), dimension(:), intent(in) :: x
   
   integer :: n, j
   real(r8) :: ans(2)
   n = size(x)
   ans(1) = x(1)
   ans(2) = x(1)
   
   do j = 2, n
      if(ans(1)>x(j))then
         ans(1) = x(j)
      endif
      
      if(ans(2)<x(j))then
         ans(2)=x(j)
      endif   
   enddo
   return
   
   end function minmax
!-------------------------------------------------------------------------------   
   subroutine cumsum_v(x, y)
   implicit none
   real(r8), dimension(:), intent(in)  :: x
   real(r8), dimension(:), intent(out) :: y
   
   integer :: n
   integer :: j
   SHR_ASSERT_ALL((size(x)   == size(y)),        errMsg(__FILE__,__LINE__))
   
   n = size(x)
   
   y(1)=x(1)
   do j = 2, n
     y(j) = y(j-1)+x(j)
   enddo
   
   end subroutine cumsum_v
!-------------------------------------------------------------------------------   
   subroutine cumsum_m(x, y, idim)
   implicit none
   real(r8), dimension(:,:), intent(in)  :: x
   real(r8), dimension(:,:), intent(out) :: y
   integer , optional,     intent(in)  :: idim
   
   integer :: n
   integer :: j
   integer :: idim_loc
   
   if(present(idim))idim_loc=idim
   
   SHR_ASSERT_ALL((size(x,1)   == size(y,1)),        errMsg(__FILE__,__LINE__))
   SHR_ASSERT_ALL((size(x,2)   == size(y,2)),        errMsg(__FILE__,__LINE__))
   
   
   if(idim_loc == 1)then
     !summation along dimension 1
     n = size(x,2)
     do j = 1, n
       call cumsum_v(x(:,j),y(:,j))
     enddo
   else
     !summation along dimension 2
     n = size(x,1)
     do j = 1, n
       call cumsum_v(x(j,:),y(j,:))
     enddo
     
   endif
   
   end subroutine cumsum_m
   
!-------------------------------------------------------------------------------   
   subroutine cumdif(x, y)
   implicit none
   real(r8), dimension(:), intent(in)  :: x
   real(r8), dimension(:), intent(out) :: y
   
   integer :: n
   integer :: j
   
   SHR_ASSERT_ALL((size(x)   == size(y)),        errMsg(__FILE__,__LINE__))
   n = size(x)
   call diff(x,y(2:n))
   y(1)=x(1)
   
   end subroutine cumdif
!-------------------------------------------------------------------------------   
   
   subroutine diff(x,y)
   implicit none
   real(r8), dimension(:), intent(in)  :: x
   real(r8), dimension(:), intent(out) :: y
   
   integer :: n
   integer :: j
   SHR_ASSERT_ALL((size(x)   == size(y)+1),        errMsg(__FILE__,__LINE__))
   
   n = size(x)
   do j = 2, n
     y(j-1) = x(j)-x(j-1)
   enddo
   end subroutine diff   
  
!------------------------------------------------------------------------------- 
   function safe_div(a,b)result(ans)
   !
   !DESCRIPTION
   !avoid division by zero when calculate a/b
   implicit none
   real(r8), intent(in) :: a
   real(r8), intent(in) :: b
   
   real(r8) :: ans
   if(abs(b)<1.e-40_r8)then 
     ans = a * b / (b**2._r8+1.e-40_r8)
   else
     ans = a/b
   endif
   return
   end function safe_div

!--------------------------------------------------------------------------------
  function dot_sum(x,y)result(ans)
  !
  !DESCRIPTIONS
  ! calculate the dot product 
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  real(r8), dimension(:), intent(in) :: x
  real(r8), dimension(:), intent(in) :: y

  integer  :: n, j
  real(r8) :: ans
  SHR_ASSERT_ALL((size(x)           == size(y)), errMsg(__FILE__,__LINE__)) 

  n = size(x)
  ! use subroutine from blas
  !DOUBLE PRECISION FUNCTION ddot(N,DX,INCX,DY,INCY)
  !
  ans=dot_product(x,y)
  !ans = 0._r8
  !do j = 1, n
  !  ans = ans + x(j)*y(j)
  !enddo 
  end function dot_sum
!--------------------------------------------------------------------------------  
  function addone(a)result(ans)
  !
  !return a variable with a + 1
  use shr_kind_mod, only: r8 => shr_kind_r8  
  implicit none
  integer, intent(inout) :: a
  
  integer :: ans
  a = a + 1
  ans = a
  end function
  
!--------------------------------------------------------------------------------
  subroutine asc_sort_vec(zvec)
  !
  ! sort an array into ascending order
  implicit none
  real(r8), dimension(:), intent(inout) :: zvec
  
  
  integer :: n, j, k
  logical :: lswap
  
  n = size(zvec)
  
  do j = 1, n
    lswap=.false.
    do k = 2, n-j+1
      if(zvec(k)<zvec(k-1))then
        lswap=.true.
        call swap_r(zvec(k),zvec(k-1))
      endif
    enddo
    if(.not. lswap)exit
  enddo
  
  end subroutine asc_sort_vec
  
!--------------------------------------------------------------------------------  
  function is_bounded(x, xl, xr)result(ans)
  !
  ! test if x is bounded within xl and xr
  implicit none
  real(r8), intent(in) :: x, xl, xr 
  logical :: ans
  if(x>=xl .and. x<=xr)then
    ans = .true.
  else
    ans = .false.
  endif
  end function is_bounded  
end module MathfuncMod
