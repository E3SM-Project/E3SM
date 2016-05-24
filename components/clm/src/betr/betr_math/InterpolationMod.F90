module InterpolationMod
#include "shr_assert.h"
  !
  ! !DESCRIPTION:
  ! subroutines to do polynomial interpolation
  ! author: Jinyun Tang, Sep, 2014

  ! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
  use abortutils            , only : endrun
  use shr_log_mod ,   only: errMsg => shr_log_errMsg
  implicit none
  private
  save

  public :: Lagrange_interp
  public :: pchip_polycc
  public :: pchip_interp
contains

  !-------------------------------------------------------------------------------
  subroutine Lagrange_interp(pn, x, y, xi, yi)
    !
    ! !DESCRIPTION:
    ! do order pn lagrangian interpolation
    implicit none
    ! !ARGUMENTS:
    integer,  intent(in)                :: pn    !order of interpolation
    real(r8), dimension(:), intent(in)  :: x   !location of data
    real(r8), dimension(:), intent(in)  :: y   !value of data
    real(r8), dimension(:), intent(in)  :: xi  !target points to be evaluated
    real(r8), dimension(:), intent(out) :: yi !target values

    ! !LOCAL VARIABLES:
    integer :: k, ni, nx
    integer :: pos, disp, disp1

    SHR_ASSERT_ALL((ubound(x) == ubound(y)), errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(xi) == ubound(yi)), errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((ubound(x) >= pn+1), errMsg(__FILE__,__LINE__))

    ni  = size(xi)
    nx = size(x)
    disp=int((pn+1)*0.5_r8+1.e-8_r8)
    !get the half size of the local window
    if(mod(pn,2)==0)then
       disp1=disp
    else
       disp1=disp-1
    endif

    do k = 1, ni
       ! find the position of z in array x
       pos = find_idx(x, xi(k))
       if(pos == -100) then
          !left boundary
          yi(k) = y(1)
       elseif(pos == -200) then
          !right boundary
          yi(k) = y(nx)
       else
          ! call function Lagrange
          if (pos <= disp1) then
             yi(k) = Lagrange_poly(pn, x(1:pn+1), y(1:pn+1), xi(k))
          else if (pos >= nx-disp) then
             yi(k) = Lagrange_poly(pn, x(nx-pn:nx), y(nx-pn:nx), xi(k))
          else
             yi(k) = Lagrange_poly(pn, x(pos-disp1:pos+disp), y(pos-disp1:pos+disp), xi(k))
          end if
       endif
    enddo

  end subroutine Lagrange_interp

  !-------------------------------------------------------------------------------
  function Lagrange_poly(pn, xvect, yvect, z)result(Pz)
    !
    ! !DESCRIPTION:
    ! do lagrangian interpolation at order pn
    !
    implicit none
    ! !ARGUMENTS:
    integer, intent(in)                :: pn   ! Order of Interpolation Polynomial
    real(r8), dimension(:), intent(in) :: xvect, yvect  ! vectors of known data: x,y-values
    real(r8), intent(in)               :: z   ! the target point "z"

    ! !LOCAL VARIABLES:
    integer   :: i, j, n
    real(r8)  :: L(pn+1)    ! Lagrange cardinal function
    real(r8)  :: Pz    ! target value

    SHR_ASSERT_ALL((size(xvect) == size(yvect)), errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((size(xvect) == pn+1), errMsg(__FILE__,__LINE__))

    ! n = number of data points:length of each data vector
    n = size(xvect)
    ! Initializations of Pz and L
    Pz = 0._r8 ! initializing the polynomia value at z
    L(:)  = 1._r8  ! initalizing the vector of cardinal functions to 1
    ! Performing the interpolation
    do i = 1, n
       do j = 1, n
          if (i /= j) then
             ! part of L(i)
             L(i) = ( (z - xvect(j)) / (xvect(i) - xvect(j)) )* L(i)
          end if
       end do
       Pz = Pz + L(i)*yvect(i) ! update Pz ~ f(z)
    end do
  end function Lagrange_poly
  !------------------------------------------------------------
  function find_idx(xvect, x)result(k)
    !
    ! !DESCRIPTION:
    ! locate the position of x in xvect
    !
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: xvect       ! vector of x-values
    real(r8), intent(in)               :: x

    integer :: i, k, n

    ! array dimension
    n = size(xvect)



    if(x<xvect(1))then
       k = -100  !beyond left boundary
    elseif(x>xvect(n))then
       k = -200  !beyond right boundary
    elseif(x==xvect(1))then
       k=1
    elseif(x==xvect(n))then
       k=n-1
    else
       ! find index k so that x[k] < x < x[k+1]
       do i = 1, n-1
          if ((xvect(i) <= x) .and. (x < xvect(i+1))) then
             k = i
             exit
          end if
       end do
    endif


  end function find_idx
  !------------------------------------------------------------
  subroutine pchip_polycc(x, fx, di, region)
    !
    ! DESCRIPTION
    ! Given the data, generate the coefficients of the monotonic cubic
    ! polynomials
    ! Ref, Fritsch and Carlson, 1980
    !
    ! !USES:
    use MathfuncMod, only : diff
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: x
    real(r8), dimension(:), intent(in) :: fx
    real(r8), dimension(:), intent(out):: di
    integer,  optional,     intent(in) :: region

    ! !LOCAL VARIABLES:
    real(r8), allocatable :: h(:)
    real(r8), allocatable :: df(:)
    real(r8), allocatable :: slp(:)
    real(r8) :: alpha, beta, tao, rr
    integer :: region_loc
    integer :: n, j

    SHR_ASSERT_ALL((size(x) == size(fx)), errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((size(x) == size(di)), errMsg(__FILE__,__LINE__))
    region_loc=2
    if(present(region))region_loc=region

    n = size(x)
    allocate(h(n-1))
    allocate(df(n-1))
    allocate(slp(n-1))
    !get interval length
    call diff(x, h)
    !get function step
    call diff(fx, df)

    !get slope
    do j = 1, n-1
       slp(j)=df(j)/h(j)
    enddo

    !get di
    di(:) = 0._r8

    j = 1
    di(j)=(fx(j+1)+fx(j+2)-2*fx(1))/(2*h(j)+h(j+1))
    do j = 2, n-1
       di(j)=(fx(j+1)-fx(j-1))/(h(j)+h(j-1))
    enddo
    j = n
    di(j)=(2._r8*fx(j)-(fx(j-1)+fx(j-2)))/(2._r8*h(j-1)+h(j-2))

    !enforce the sign condition
    if(slp(1)*di(1)<=0._r8)then
       di(1)=0._r8
    endif

    if(slp(n-1)*di(n)<=0)then
       di(n)=0._r8
    endif

    !enforce the range 2 constraint

    do j = 1, n-1
       if(abs(slp(j))<=1.e-16_r8)then
          di(j)=0._r8
          di(j+1)=0._r8
       else
          alpha=di(j)/slp(j)
          beta =di(j+1)/slp(j)
          select case (region_loc)
          case (1)
             rr=beta/alpha
             if(rr>1._r8)then
                if(beta>3._r8)then
                   beta=3._r8
                   alpha=beta/rr
                   di(j)=slp(j)*alpha
                   di(j+1)=slp(j)*beta
                endif
             else
                if(alpha>3._r8)then
                   alpha=3._r8
                   beta=rr*alpha
                   di(j)=slp(j)*alpha
                   di(j+1)=slp(j)*beta
                endif
             endif
          case (2)
             tao=3._r8/sqrt(alpha*alpha+beta*beta)
             if(tao<1._r8)then
                di(j)=tao*di(j)
                di(j+1)=tao*di(j+1)
             endif
          case (3)
             if(alpha+beta>3._r8)then
                if(alpha>0._r8)then
                   rr=beta/alpha
                   alpha=3._r8/(1._r8+rr)
                   beta=alpha*rr
                   di(j)=slp(j)*alpha
                   di(j+1)=slp(j)*beta
                else
                   beta=3._r8
                   di(j+1)=slp(j)*beta
                endif
             endif
          case (4)
             if(alpha>0._r8)then
                rr=beta/alpha
                if(rr>=1._r8)then
                   if(2._r8*alpha+beta>3._r8)then
                      alpha=3._r8/(2._r8+rr)
                      beta=alpha*rr
                      di(j)=slp(j)*alpha
                      di(j+1)=slp(j)*beta;
                   endif
                else
                   if(alpha+2._r8*beta>3._r8)then
                      alpha=3._r8/(1._r8+2._r8*rr)
                      beta=alpha*rr
                      di(j)=slp(j)*alpha
                      di(j+1)=slp(j)*beta
                   endif
                endif
             else
                if(beta>3._r8)then
                   beta=3._r8
                   di(j+1)=slp(j)*beta
                endif
             endif
          case default
             call endrun(msg='an constraint region must be specified for pchip_polycc '//errMsg(__FILE__, __LINE__))
          end select

       endif
    enddo

    deallocate(h)
    deallocate(df)
    deallocate(slp)
  end subroutine pchip_polycc
  !------------------------------------------------------------

  subroutine pchip_interp(x, fx, di, xi, yi)

    ! !DESCRIPTION:
    ! do monotonic cubic spline interpolation
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: x
    real(r8), dimension(:), intent(in) :: fx
    real(r8), dimension(:), intent(in) :: di
    real(r8), dimension(:), intent(in) :: xi
    real(r8), dimension(:), intent(out) :: yi

    ! !LOCAL VARIABLES:
    real(r8) :: h, t1, t2
    real(r8) :: h1x, h2x, h3x, h4x
    integer  :: n, j
    integer  :: id

    SHR_ASSERT_ALL((size(x) == size(fx)),  errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((size(x) == size(di)),  errMsg(__FILE__,__LINE__))
    SHR_ASSERT_ALL((size(xi) == size(yi)), errMsg(__FILE__,__LINE__))

    n=size(xi)  !total number of points to be interpolated

    yi(:)=0._r8
    do j = 1, n

       id=find_idx(x,xi(j))
       h=x(id+1)-x(id)
       t1=(x(id+1)-xi(j))/h
       t2=(xi(j)-x(id))/h

       h1x=phi(t1)
       h2x=phi(t2)
       h3x=-h*psi(t1)
       h4x=h*psi(t2)
       yi(j)=fx(id)*h1x+fx(id+1)*h2x+di(id)*h3x+di(id+1)*h4x
    enddo

  contains

    function phi(t) result(fval)
      implicit none
      real(r8), intent(in) :: t

      real(r8) :: fval

      fval=(3._r8-2._r8*t)*t*t

    end function phi


    function psi(t) result(fval)
      implicit none
      real(r8), intent(in) :: t

      real(r8) :: fval
      fval=t*t*(t-1._r8)
    end function psi
  end subroutine pchip_interp

end module InterpolationMod
