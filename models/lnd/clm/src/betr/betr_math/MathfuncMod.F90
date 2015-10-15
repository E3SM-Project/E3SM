module MathfuncMod
#include "shr_assert.h"
  ! !DESCRIPTION:
  ! mathematical functions for some elementary manipulations
  ! History: Created by Jinyun Tang
  !
  ! !USES:
  use shr_kind_mod    , only : r8 => shr_kind_r8
  use clm_varctl      , only : iulog
  use abortutils      , only : endrun
  use shr_log_mod     , only : errMsg => shr_log_errMsg
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
  public :: minp
  public :: pd_decomp
  public :: num2str
  interface cumsum
     module procedure cumsum_v, cumsum_m
  end interface cumsum
  interface swap
     module procedure swap_i, swap_r, swap_rv
  end interface swap
contains
  !-------------------------------------------------------------------------------
  function heviside(x)result(ans)
    !
    ! !DESCRIPTION:
    !  heviside function
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: x
    ! !LOCAL VARIABLES:
    real(r8) :: ans

    if(x>0._r8)then
       ans = 1._r8
    else
       ans = 0._r8
    endif
  end function heviside


  !-------------------------------------------------------------------------------
  subroutine swap_i(a,b)
    !
    ! !DESCRIPTION:
    ! swap two integers
    implicit none
    ! !ARGUMENTS:
    integer, intent(inout) :: a, b

    ! !LOCAL VARIABLES:
    integer :: c

    c = a
    a = b
    b = c

  end subroutine swap_i
  !-------------------------------------------------------------------------------
  subroutine swap_r(a,b)
    !
    ! !DESCRIPTION:
    ! swap two real numbers
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(inout) :: a, b

    ! !LOCAL VARIABLES:
    real(r8) :: c

    c = a
    a = b
    b = c

  end subroutine swap_r
  !-------------------------------------------------------------------------------
  subroutine swap_rv(a,b)
    !
    ! !DESCRIPTION:
    ! swap two vectors
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(inout) :: a, b
    ! !LOCAL VARIABLES:
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
    !
    ! !DESCRIPTION:
    !returnd the minimum and maximum of the input vector
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: x

    ! !LOCAL VARIABLES:
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
    !
    ! !DESCRIPTION:
    ! cumulative sum of a vector x
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in)  :: x  !input vector
    real(r8), dimension(:), intent(out) :: y  !sum

    ! !LOCAL VARIABLES:
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
    !
    ! !DESCRIPTION:
    ! do cumulative summation for maxtrix x along dimnension idim
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:,:), intent(in)  :: x    !input array
    real(r8), dimension(:,:), intent(out) :: y    !output cum sum
    integer , optional,     intent(in)    :: idim !dimension to be summed

    ! !LOCAL VARIABLES:
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
    !
    ! !DESCRIPTION:
    ! do nearest neighbor finite difference
    !
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in)  :: x   !input array
    real(r8), dimension(:), intent(out) :: y   !output dif

    ! !LOCAL VARIABLES:
    integer :: n
    integer :: j

    SHR_ASSERT_ALL((size(x)   == size(y)),        errMsg(__FILE__,__LINE__))
    n = size(x)
    call diff(x,y(2:n))
    y(1)=x(1)

  end subroutine cumdif
  !-------------------------------------------------------------------------------

  subroutine diff(x,y)
    !
    ! !DESCRIPTION:
    ! do nearest neighbor forward difference
    !
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in)  :: x  !input array
    real(r8), dimension(:), intent(out) :: y  !output array

    integer :: n
    integer :: j
    SHR_ASSERT_ALL((size(x)   == size(y)+1),        errMsg(__FILE__,__LINE__))

    n = size(x)
    do j = 2, n
       y(j-1) = x(j)-x(j-1)
    enddo
  end subroutine diff

  !-------------------------------------------------------------------------------
  function safe_div(a,b,eps)result(ans)
    !
    ! !DESCRIPTION:
    ! avoid division by zero when calculate a/b
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in)           :: a   !numerator
    real(r8), intent(in)           :: b   !denominator
    real(r8), optional, intent(in) :: eps !screening threshold
    !
    ! !LOCAL VARIABLES:
    real(r8) :: ans
    real(r8) :: loc_eps
    if(present(eps))then
       loc_eps=eps
    else
       loc_eps=1.e-40_r8
    endif
    if(abs(b)<loc_eps)then
       ans = a * b / (b**2._r8+loc_eps)
    else
       ans = a/b
    endif
    return
  end function safe_div

  !--------------------------------------------------------------------------------
  function dot_sum(x,y)result(ans)
    !
    ! !DESCRIPTION:
    ! calculate the dot product
    !
    ! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: x
    real(r8), dimension(:), intent(in) :: y
    ! !LOCAL VARIABLES:
    integer  :: n, j
    real(r8) :: ans
    SHR_ASSERT_ALL((size(x)           == size(y)), errMsg(__FILE__,__LINE__))

    n = size(x)
    ! use subroutine from blas
    !DOUBLE PRECISION FUNCTION ddot(N,DX,INCX,DY,INCY)
    !
    ans=dot_product(x,y)

  end function dot_sum
  !--------------------------------------------------------------------------------
  function addone(a)result(ans)
    ! !DESCRIPTION:
    ! return a variable with a + 1
    !
    ! !USES:
    use shr_kind_mod, only: r8 => shr_kind_r8
    implicit none
    ! !ARGUMENTS:
    integer, intent(inout) :: a
    ! !LOCAL VARIABLES:
    integer :: ans

    a = a + 1
    ans = a
  end function addone

  !--------------------------------------------------------------------------------
  subroutine asc_sort_vec(zvec)
    !
    ! !DESCRIPTION:
    ! sort an array into ascending order
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(inout) :: zvec
    ! !LOCAL VARIABLES:
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
    ! !DESCRIPTION:
    ! test if x is bounded within xl and xr
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: x, xl, xr

    ! !LOCAL VARIABLES:
    logical :: ans
    if(x>=xl .and. x<=xr)then
       ans = .true.
    else
       ans = .false.
    endif
  end function is_bounded

  !--------------------------------------------------------------------------------
  function minp(p,v)result(ans)
    !
    ! !DESCRIPTION:
    !find the minimum of the nonzero p entries, with the entry determined by
    !nonzero values of v

    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: p
    real(r8), dimension(:), intent(in) :: v
    ! !LOCAL VARIABLES:
    integer  :: j, sz
    real(r8) :: ans      !(<=1._r8)

    SHR_ASSERT_ALL((size(p)           == size(v)), errMsg(__FILE__,__LINE__))

    sz = size(p)
    ans = 1._r8
    do j = 1, sz
       if(v(j)/=0._r8)then
          ans = min(ans, p(j))
       endif
    enddo
  end function minp

  !--------------------------------------------------------------------------------
  subroutine pd_decomp(m, n, A, AP, AD)
    !
    ! !DESCRIPTION:
    !separate a input matrix A into AP and AD with positive
    !and negative entries respectively.

    implicit none
    ! !ARGUMENTS:
    integer  , intent(in) :: n, m
    real(r8) , intent(in) :: A(1: ,  1: )
    real(r8) , intent(out):: AP(1: , 1: )
    real(r8) , intent(out):: AD(1: , 1: )

    ! !LOCAL VARIABLES:
    integer :: i, j


    SHR_ASSERT_ALL((ubound(A)           == (/m,n/)), errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(AP)          == (/m,n/)), errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(AD)          == (/m,n/)), errMsg(__FILE__,__LINE__))

    AP(:,:) = 0._r8
    AD(:,:) = 0._r8

    where(A>0._r8)
       AP=A
    elsewhere
       AD=A
    endwhere
  end subroutine pd_decomp
  !--------------------------------------------------------------------------------
  
  function num2str(a,fmt)result(ans)
    !
    ! !DESCRIPTION:
    !turn a number into a string using the specified format
    implicit none
    ! !ARGUMENTS:
    integer, intent(in) :: a
    character(len=*), intent(in) :: fmt

    ! !LOCAL VARIABLES:
    character(len=32) :: ans
    character(len=32) :: str

    write(str,fmt)a
    ans =  trim(adjustl(str))
  end function num2str
end module MathfuncMod
