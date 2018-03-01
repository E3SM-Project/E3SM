module MathfuncMod
#include "bshr_assert.h"
  ! !DESCRIPTION:
  ! mathematical functions for some elementary manipulations
  ! History: Created by Jinyun Tang
  !
  ! !USES:
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use betr_ctrl     , only : iulog  => biulog
  use bshr_log_mod  , only : errMsg => shr_log_errMsg

  implicit none

  private

  character(len=*), parameter :: mod_filename = &
       __FILE__

  public :: heviside
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
  public :: fpmax
  public :: bisnan
  public :: apvb
  public :: countelm
  interface apvb
    module procedure apvb_v, apvb_s
  end interface apvb
  interface cumsum
     module procedure cumsum_v, cumsum_m
  end interface cumsum
  interface swap
     module procedure swap_i, swap_r, swap_rv
  end interface swap

  !law of minimum flux based back tracing tools
  !for ODE integration, Tang and Riley, BG, 2015.
  type, public :: lom_type
  contains
    procedure, public :: calc_state_pscal
    procedure, public :: calc_reaction_rscal
    procedure, public :: apply_reaction_rscal
  end type lom_type
contains
  !-------------------------------------------------------------------------------
  function countelm(ibeg, iend, igap)result(ans)
  implicit none
  integer, intent(in) :: ibeg   !begin of the counter
  integer, intent(in) :: iend   !end of the counter
  integer, optional, intent(in) :: igap   !gap between two consecutive numbers

  integer :: ans

  if(present(igap))then
    ans = (iend-ibeg)/igap + 1
  else
    ans = (iend-ibeg) + 1
  endif
  end function countelm
  !-------------------------------------------------------------------------------
  subroutine apvb_s(a, brr, sgn)
  !
  ! DESCRIPTION
  ! compute a = a + sum(brr)
  implicit none
  real(r8), intent(inout) :: a
  real(r8), intent(in) :: brr
  real(r8), optional, intent(in) :: sgn
  real(r8) :: sgn_loc

  sgn_loc = 1._r8
  if(present(sgn))sgn_loc = sgn

  a = a + sgn_loc * brr

  end subroutine apvb_s

  !-------------------------------------------------------------------------------
  subroutine apvb_v(a, brr, sgn)
  !
  ! DESCRIPTION
  ! compute a = a + sum(brr)
  implicit none
  real(r8), intent(inout) :: a
  real(r8), dimension(:), intent(in) :: brr
  real(r8), optional, intent(in) :: sgn
  real(r8) :: sgn_loc

  sgn_loc = 1._r8
  if(present(sgn))sgn_loc = sgn

  a = a + sgn_loc * sum(brr)

  end subroutine apvb_v
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
  subroutine swap_rv(a,b, betr_status)
    !
    ! !DESCRIPTION:
    ! swap two vectors
    use BetrStatusType         , only : betr_status_type
    use betr_constants         , only : betr_errmsg_len
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(inout) :: a, b
    type(betr_status_type), intent(out)   :: betr_status
    ! !LOCAL VARIABLES:
    real(r8), dimension(size(a)) :: c
    character(len=betr_errmsg_len) :: msg
    integer :: n

    if(size(a)/=size(b))then
       write(msg,*)'the input vectors are not of same size in swap_rv'
       msg = trim(msg)//new_line('A')//'stop in '//errmsg(mod_filename, __LINE__)
       call betr_status%set_msg(msg=msg,err=-1)
       return
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
  subroutine cumsum_v(bstatus, x, y)
    !
    ! !DESCRIPTION:
    ! cumulative sum of a vector x
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in)  :: x  !input vector
    real(r8), dimension(:), intent(out) :: y  !sum
    type(betr_status_type), intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    integer :: n
    integer :: j

    call bstatus%reset()
    SHR_ASSERT_ALL((size(x)   == size(y)), errMsg(mod_filename,__LINE__), bstatus)
    if(bstatus%check_status())return
    n = size(x)

    y(1)=x(1)
    do j = 2, n
       y(j) = y(j-1)+x(j)
    enddo

  end subroutine cumsum_v
  !-------------------------------------------------------------------------------
  subroutine cumsum_m(bstatus, x, y, idim)
    !
    ! !DESCRIPTION:
    ! do cumulative summation for maxtrix x along dimnension idim
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:,:), intent(in)  :: x    !input array
    real(r8), dimension(:,:), intent(out) :: y    !output cum sum
    integer , optional,     intent(in)    :: idim !dimension to be summed
    type(betr_status_type), intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    integer :: n
    integer :: j
    integer :: idim_loc

    call bstatus%reset()
    if(present(idim))idim_loc=idim

    SHR_ASSERT_ALL((size(x,1)   == size(y,1)),        errMsg(mod_filename,__LINE__), bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((size(x,2)   == size(y,2)),        errMsg(mod_filename,__LINE__), bstatus)
    if(bstatus%check_status())return
    if(idim_loc == 1)then
       !summation along dimension 1
       n = size(x,2)
       do j = 1, n
         call cumsum_v(bstatus, x(:,j),y(:,j))
         if(bstatus%check_status())return
       enddo
    else
       !summation along dimension 2
       n = size(x,1)
       do j = 1, n
          call cumsum_v(bstatus, x(j,:),y(j,:))
          if(bstatus%check_status())return
       enddo
    endif

  end subroutine cumsum_m

  !-------------------------------------------------------------------------------
  subroutine cumdif(x, y, bstatus)
    !
    ! !DESCRIPTION:
    ! do nearest neighbor finite difference
    !
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in)  :: x   !input array
    real(r8), dimension(:), intent(out) :: y   !output dif
    type(betr_status_type), intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    integer :: n
    integer :: j
    call bstatus%reset()
    SHR_ASSERT_ALL((size(x)   == size(y)),  errMsg(mod_filename,__LINE__), bstatus)
    if(bstatus%check_status())return

    n = size(x)
    call diff(x,y(2:n), bstatus)
    if(bstatus%check_status())return
    y(1)=x(1)

  end subroutine cumdif
  !-------------------------------------------------------------------------------

  subroutine diff(x,y, bstatus)
    !
    ! !DESCRIPTION:
    ! do nearest neighbor forward difference
    !
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in)  :: x  !input array
    real(r8), dimension(:), intent(out) :: y  !output array
    type(betr_status_type), intent(out) :: bstatus

    integer :: n
    integer :: j

    call bstatus%reset()
    SHR_ASSERT_ALL((size(x)   == size(y)+1),        errMsg(mod_filename,__LINE__), bstatus)
    if(bstatus%check_status())return

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
    real(r8),           intent(in) :: a   !numerator
    real(r8),           intent(in) :: b   !denominator
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
  function dot_sum(x,y, bstatus)result(ans)
    !
    ! !DESCRIPTION:
    ! calculate the dot product
    !
    ! !USES:
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: x
    real(r8), dimension(:), intent(in) :: y
    type(betr_status_type), optional, intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    integer  :: n, j
    real(r8) :: ans
    type(betr_status_type) :: bstatus1

    call bstatus1%reset()
    SHR_ASSERT_ALL((size(x)  == size(y)), errMsg(mod_filename,__LINE__), bstatus1)
    if(present(bstatus))then
      call bstatus%reset()
      call bstatus%set_msg(bstatus1%print_msg(),bstatus1%print_err())
      if(bstatus%check_status())return
    endif
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
    ! don't use it directly to refer an array element in the following
    ! arr(addone(a))
    ! rather you should use
    ! id=addone(a); arr(id)
    ! !USES:

    implicit none
    ! !ARGUMENTS:
    integer, intent(inout) :: a
    ! !LOCAL VARIABLES:
    integer :: ans

    ans = a + 1
    a = ans
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
  function minp(p,v, bstatus)result(ans)
    !
    ! !DESCRIPTION:
    !find the minimum of the nonzero p entries, with the entry determined by
    !nonzero values of v
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: p
    real(r8), dimension(:), intent(in) :: v
    ! !LOCAL VARIABLES:
    integer  :: j, sz
    real(r8) :: ans      !(<=1._r8)
    type(betr_status_type), intent(out) :: bstatus

    call bstatus%reset()
    SHR_ASSERT_ALL((size(p)   == size(v)), errMsg(mod_filename,__LINE__), bstatus)
    if(bstatus%check_status())return
    sz = size(p)
    ans = 1._r8
    do j = 1, sz
       if(v(j)/=0._r8)then
          ans = min(ans, p(j))
       endif
    enddo
  end function minp

  !--------------------------------------------------------------------------------
  subroutine pd_decomp(m, n, A, AP, AD, bstatus)
    !
    ! !DESCRIPTION:
    !separate a input matrix A into AP and AD with positive
    !and negative entries respectively.
    use BetrStatusType, only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    integer  , intent(in) :: n, m
    real(r8) , intent(in) :: A(1: ,  1: )
    real(r8) , intent(out):: AP(1: , 1: )
    real(r8) , intent(out):: AD(1: , 1: )
    type(betr_status_type), intent(out) :: bstatus

    ! !LOCAL VARIABLES:
    integer :: i, j

    call bstatus%reset()
    SHR_ASSERT_ALL((ubound(A)           == (/m,n/)), errMsg(mod_filename,__LINE__), bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(AP)          == (/m,n/)), errMsg(mod_filename,__LINE__), bstatus)
    if(bstatus%check_status())return
    SHR_ASSERT_ALL((ubound(AD)          == (/m,n/)), errMsg(mod_filename,__LINE__), bstatus)
    if(bstatus%check_status())return

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
    integer,          intent(in) :: a
    character(len=*), intent(in) :: fmt

    ! !LOCAL VARIABLES:
    character(len=32) :: ans
    character(len=32) :: str

    write(str,fmt)a
    ans =  trim(adjustl(str))
  end function num2str

  !-------------------------------------------------------------------------------
  subroutine calc_state_pscal(this, nprimvars, dtime, ystate, p_dt,  d_dt, pscal, lneg, bstatus)
    !
    ! !DESCRIPTION:
    ! calcualte limiting factor from each primary state variable
    !
    use BetrstatusType     , only : betr_status_type
    use betr_constants     , only : betr_errmsg_len
    implicit none
    ! !ARGUMENTS:
    class(lom_type), intent(in) :: this
    integer,  intent(in)  :: nprimvars
    real(r8), intent(in)  :: dtime
    real(r8), intent(in)  :: ystate(1:nprimvars)
    real(r8), intent(in)  :: p_dt(1:nprimvars)
    real(r8), intent(in)  :: d_dt(1:nprimvars)
    real(r8), intent(out) :: pscal(1:nprimvars)
    logical,  intent(out) :: lneg
    type(betr_status_type), intent(out) :: bstatus
    character(len=betr_errmsg_len) :: msg

    ! !LOCAL VARIABLES:
    real(r8) :: yt
    integer  :: j
    real(r8),parameter :: p_par=0.999_r8

    call bstatus%reset()
    lneg =.false.

    do j = 1, nprimvars
       yt = ystate(j) + (p_dt(j)+d_dt(j))*dtime
       if(yt<0._r8)then
          pscal(j) = -(p_dt(j)*dtime+ystate(j))/(dtime*d_dt(j))*p_par
          lneg=.true.
          if(pscal(j)<0._r8)then
             msg = 'ngeative p in calc_state_pscal'//errmsg(mod_filename, __LINE__)
             call bstatus%set_msg(msg,err=-1)
             if(bstatus%check_status())return
          endif
       else
          pscal(j) = 1._r8
       endif
    enddo
  end subroutine calc_state_pscal

  !-------------------------------------------------------------------------------
  subroutine calc_reaction_rscal(this, nprimvars, nr, pscal, cascade_matrixd, rscal, bstatus)
    !
    ! !DESCRIPTION:
    ! calcualte limiting factor for each reaction
    ! !USES:
    use BetrstatusType     , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    class(lom_type), intent(in) :: this
    integer , intent(in) :: nprimvars
    integer , intent(in) :: nr
    real(r8), intent(in) :: pscal(1:nprimvars)
    real(r8), intent(in) :: cascade_matrixd(1:nprimvars, 1:nr)
    real(r8), intent(out):: rscal(1:nr)
    type(betr_status_type), intent(out) :: bstatus
    ! !LOCAL VARIABLES:
    integer :: j

    call bstatus%reset()
    do j = 1, nr
      rscal(j) = minp(pscal,cascade_matrixd(1:nprimvars, j), bstatus)
      if(bstatus%check_status())return
    enddo

  end subroutine calc_reaction_rscal

  !-------------------------------------------------------------------------------
  subroutine apply_reaction_rscal(this, nr, rscal, reaction_rates)
    !
    ! !DESCRIPTION:
    ! reduce reaction rates using input scalar
    !
    implicit none
    ! !ARGUMENTS:
    class(lom_type), intent(in) :: this
    integer , intent(in)    :: nr
    real(r8), intent(in)    :: rscal(1:nr)
    real(r8), intent(inout) :: reaction_rates(1:nr)
    ! !LOCAL VARIABLES:
    integer :: j

    do j = 1, nr
       reaction_rates(j) = reaction_rates(j)*rscal(j)
    enddo
  end subroutine  apply_reaction_rscal

  !-------------------------------------------------------------------------------

  function fpmax(inval)result(ans)
  !
  ! DESCRIPTION
  ! return positive values
  implicit none
  real(r8), intent(in) :: inval

  real(r8) :: ans

  ans = max(inval, 0._r8)
  return
  end function fpmax

  !-------------------------------------------------------------------------------
  function bisnan(inval)result(ans)

  !DESCRIPTION
  !determine if the variable is nan
  implicit none
  real(r8), intent(in) :: inval

  logical :: ans

  ans = (inval/=inval)
  end function bisnan
end module MathfuncMod
