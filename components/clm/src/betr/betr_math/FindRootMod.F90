module FindRootMod
  !
  ! !DESCRIPTION:
  !  Functions to solve simple equations
  !  History: created by Jinyun Tang, 2013

  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use abortutils          , only : endrun
  use clm_varctl          , only : iulog
  use MathfuncMod         , only : is_bounded
  implicit none
  interface hybrid_findroot
     module procedure hybrid_findroot_np, hybrid_findroot_p
  end interface hybrid_findroot

  interface brent
     module procedure brent_np, brent_p
  end interface brent

contains
  !-------------------------------------------------------------------------------
  function quadrootbnd(a,b,c, xl, xr)result(x)
    !
    ! !DESCRIPTION:
    ! return a root of the qudratic equation
    ! within bound xl and xr

    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: a, b, c
    real(r8), intent(in) :: xl, xr

    ! !LOCAL VARIABLES:
    real(r8) :: x
    real(r8) :: delta
    character(len=32) :: subname ='quadrootbnd'

    delta = b * b -4._r8 * a * c
    if(delta>=0._r8)then
       x = (-b + sqrt(delta))/2._r8
       if(is_bounded(x,xl,xr))then
          return
       else
          x = (-b - sqrt(delta))/2._r8
          if(is_bounded(x,xl,xr))then
             return
          else
             write(iulog,*)'no bounded solution for the given quadratic equation'
             call endrun(msg=errmsg(__FILE__, __LINE__))
          endif
       endif
    else
       write(iulog,*)'no real solution for the given quadratic equation'
       call endrun(msg=errmsg(__FILE__, __LINE__))
    endif
    return
  end function quadrootbnd

  !-------------------------------------------------------------------------------
  
  function quadproot(a,b,c)result(x)
    !
    ! !DESCRIPTION:
    ! return positive root of the qudratic equation

    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: a, b, c

    ! !LOCAL VARIABLES:
    real(r8) :: x
    real(r8) :: delta
    character(len=32) :: subname ='quadproot'

    delta = b * b -4._r8 * a * c
    if(delta>=0._r8)then
       x = (-b + sqrt(delta))/2._r8
    else
       write(iulog,*)'no positive solution for the given quadratic equation'
       call endrun(msg=errmsg(__FILE__, __LINE__))
    endif
    return
  end function quadproot
  !===============================================================================

  function cubicrootbnd(a,b,c,d, xl, xr)result(x)
    !
    ! !DESCRIPTION:
    ! return positive root of the cubic equation
    !
    ! !USES:
    use clm_varcon , only : rpi
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: a, b, c, d
    real(r8), intent(in) :: xl, xr

    ! !LOCAL VARIABLES:
    real(r8) :: x
    real(r8) :: p, q
    real(r8) :: b1, c1, d1
    real(r8) :: n, u, f, y
    real(r8) :: delta
    character(len=32) :: subname ='cubicrootbnd'

    b1 = b/a
    c1 = c/a
    d1 = d/a

    p = c1 - b1 * b1 /3._r8
    q = d1 - b1 / 3._r8 * (c1 - 2._r8 * b1**2._r8/9._r8)

    delta =-4._r8 * p**3._r8 - 27._r8 * q ** 2._r8
    if(delta<0._r8)then
       write(iulog,*)'no real solution for the given cubic equation'
       call endrun(msg=errmsg(__FILE__, __LINE__))
    else
       n = sqrt(-4._r8*p/3._r8)
       f = -q/2._r8 * (-p/3._r8)**(-1.5_r8)
       u = acos(f)/3._r8

       y = n * cos(u)
       x = y - b1 / 3._r8
       if(is_bounded(x,xl,xr))then
          return
       else
          y = n * max(cos(u), cos(u-rpi*2._r8/3._r8))
          x = y - b1 /3._r8
          if(is_bounded(x,xl,xr))then
             return
          else
             y = n * cos(u-rpi*2._r8/3._r8)
             x = y - b1 / 3._r8
             if(is_bounded(x,xl,xr))then
                return
             else
                write(iulog,*)'no bounded solution for the given cubic equation'
                call endrun(msg=errmsg(__FILE__, __LINE__))
             endif
          endif
       endif

    endif
  end function cubicrootbnd

  !===============================================================================
  function cubicproot(a,b,c,d)result(x)
    !
    ! !DESCRIPTION:
    ! return positive root of the cubic equation
    !
    ! !USES:
    use clm_varcon , only : rpi
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: a, b, c, d
    ! !LOCAL VARIABLES:
    real(r8) :: x
    real(r8) :: p, q
    real(r8) :: b1, c1, d1
    real(r8) :: n, u, f, y
    real(r8) :: delta
    b1 = b/a
    c1 = c/a
    d1 = d/a

    p = c1 - b1 * b1 /3._r8
    q = d1 - b1 / 3._r8 * (c1 - 2._r8 * b1**2._r8/9._r8)

    delta =-4._r8 * p**3._r8 - 27._r8 * q ** 2._r8
    if(delta<0._r8)then
       write(iulog,*)'no real solution for the given cubic equation'
       call endrun(msg=errmsg(__FILE__, __LINE__))
    else
       n = sqrt(-4._r8*p/3._r8)
       f = -q/2._r8 * (-p/3._r8)**(-1.5_r8)
       u = acos(f)/3._r8

       if(u<=rpi/3._r8)then
          y = n * cos(u)
       elseif(u>rpi/3._r8 .and. u < rpi/2._r8)then
          !return the maximum of the two non-negative solutions
          y = n * max(cos(u), cos(u-rpi*2._r8/3._r8))
       else
          y = n * cos(u-rpi*2._r8/3._r8)
       endif
       x = y - b1 / 3._r8
    endif
  end function cubicproot

  !===============================================================================

  subroutine LUsolvAxr(a,r, n)
    ! !DESCRIPTION:
    !solve linear equation Ax=r, using the LU decomposition
    implicit none
    ! !ARGUMENTS:
    real(r8) , intent(inout) :: a(n,n)
    real(r8) , intent(inout) :: r(n)
    integer  , intent(in)    :: n

    ! !LOCAL VARIABLES:
    real(r8) :: d(n)
    integer :: indx(n)

    !do lu decomposition
    call ludcmp(a,indx,d,n)

    !solve for the equation

    call lubksb(a,indx,r,n)
  end subroutine LUsolvAxr
  !===============================================================================
  
  subroutine lubksb(a,indx,b,n)
    !
    ! !DESCRIPTION:
    ! Solves the set of N linear equations A X = B. Here the N x N matrix a is input, not
    ! as the original matrix A, but rather as its LU decomposition, determined by the routine
    ! ludcmp. indx is input as the permutation vector of length N returned by ludcmp. b is
    ! input as the right-hand-side vector B, also of length N, and returns with the solution vector
    ! X. a and indx are not modified by this routine and can be left in place for successive calls
    ! with different right-hand sides b. This routine takes into account the possibility that b will
    ! begin with many zero elements, so it is efficient for use in matrix inversion.

    implicit none
    ! !ARGUMENTS:
    real(r8) ,    intent(in) :: a(n,n)
    integer  ,    intent(in) :: indx(n)
    real(r8) , intent(inout) :: b(n)
    integer  ,    intent(in) :: n

    ! !LOCAL VARIABLES:
    integer :: i,ii,ll
    real(r8) :: summ

    ii=0 !When ii is set to a positive value, it will become the index
    !          of the first nonvanishing element of b. We now do
    !          the forward substitution, equation (2.3.6). The only new
    !          wrinkle is to unscramble the permutation as we go.
    do i=1,n
       ll=indx(i)
       summ=b(ll)
       b(ll)=b(i)
       if (ii /= 0) then
          summ=summ-dot_product(a(i,ii:i-1),b(ii:i-1))
       else if (summ /= 0.0) then
          ii=i !A nonzero element was encountered, so from now on we will
       end if !have to do the dot product above.
       b(i)=summ
    end do
    do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
  end subroutine lubksb
  !===============================================================================


  subroutine ludcmp(a,indx,d,n)
    ! !DESCRIPTION:
    !
    !   LU docomposition
    ! adapted from Numerical recipe, chptB2
    ! Given an N by N input matrix a, this routine replaces it by the LU decomposition of a
    ! rowwise permutation of itself. On output, a is arranged as in equation (2.3.14); indx is an
    ! output vector of length N that records the row permutation effected by the partial pivoting;
    ! d is output as �1 depending on whether the number of row interchanges was even or odd,
    ! respectively. This routine is used in combination with lubksb to solve linear equations or
    ! invert a matrix.
    !
    ! !USES:
    use MathfuncMod, only : swap
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(inout) :: a(n,n)
    integer   , intent(out) :: indx(n)
    real(r8)  , intent(out) :: d(n)
    integer  ,  intent(in)  :: n

    ! !LOCAL VARIABLES:
    real(r8), dimension(size(a,1)) :: vv !vv stores the implicit scaling of each row.
    real(r8), parameter :: TINY=1.0e-20  !A small number.
    integer :: j,imax


    d=1.0 !No row interchanges yet.
    vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling
    if (any(vv == 0.0)) then
       write(6,*)'singular matrix in ludcmp' !information.
       stop
    endif
    !There is a row of zeros.
    vv=1.0 / vv !Save the scaling.
    do j=1,n
       imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
       if (j /= imax) then !Do we need to interchange rows?
          call swap(a(imax,:),a(j,:)) !Yes, do so...
          d=-d !...and change the parity of d.
          vv(imax)=vv(j) !Also interchange the scale factor.
       end if
       indx(j)=imax
       if (a(j,j) == 0.0) a(j,j)=TINY
       !       If the pivot element is zero the matrix is singular (at least to the precision of the algorithm).
       !       For some applications on singular matrices, it is desirable to substitute TINY
       !       for zero.
       a(j+1:n,j)=a(j+1:n,j)/a(j,j) !Divide by the pivot element.
       a(j+1:n,j+1:n)=a(j+1:n,j+1:n)-outerprod(a(j+1:n,j),a(j,j+1:n))
       !Reduce remaining submatrix.
    end do
  end subroutine ludcmp
  !===============================================================================
  function imaxloc(arr)
    !
    ! !DESCRIPTION:
    ! locate the maximum in a vector

    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: arr

    ! !LOCAL VARIABLES:
    integer  :: imaxloc
    integer, dimension(1) :: imax

    imax=maxloc(arr(:))
    imaxloc=imax(1)

  end function imaxloc

  !===============================================================================

  function outerprod(a,b)
    ! !DESCRIPTION:
    ! do out product of two vectors

    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:), intent(in) :: a,b

    ! !LOCAL VARIABLES:
    real(r8), dimension(size(a),size(b)) :: outerprod

    outerprod = spread(a,dim=2,ncopies=size(b)) * &
         spread(b,dim=1,ncopies=size(a))

  end function outerprod
  !-------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: brent
  !
  ! !INTERFACE:

  subroutine brent_p(x, x1,x2,f1, f2, macheps, tol, pp, func)

    !
    !!DESCRIPTION:
    !Use Brent's method to find the root of a single variable function func, which is known to exist between x1 and x2.
    !The found root will be updated until its accuracy is tol.

    !!REVISION HISTORY:
    !Dec 14/2012: Jinyun Tang, modified from numerical recipes in F90 by press et al. 1188-1189
    !
    !!USES:


    !
    !!ARGUMENTS:
    implicit none
    real(r8), intent(in) :: x1, x2, f1, f2    !minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in) :: macheps           !machine precision
    integer,  intent(in) :: pp                !index argument used by subroutine func
    real(r8), intent(in) :: tol               !the error tolerance
    real(r8), intent(out):: x                !indepedent variable of the single value function func(x)

    interface
       subroutine func(x,f, pp)
         use shr_kind_mod       , only : r8 => shr_kind_r8
         implicit none
         real(r8), intent(in)  :: x
         real(r8), intent(out) :: f
         integer,  intent(in) :: pp
       end subroutine func
    end interface

    ! !LOCAL VARIABLES:
    integer, parameter :: ITMAX = 40            !maximum number of iterations
    integer :: iter
    real(r8)  :: a,b,c,d,e,fa,fb,fc,p,q,r,s,xm,tol1


    a=x1
    b=x2
    fa=f1
    fb=f2
    if((fa > 0._r8 .and. fb > 0._r8).or.(fa < 0._r8 .and. fb < 0._r8))then
       write(iulog,*) 'root must be bracketed for brent'
       write(iulog,*) 'a=',a,' b=',b,' fa=',fa,' fb=',fb
       call endrun(msg=errmsg(__FILE__, __LINE__))
    endif
    c=b
    fc=fb
    iter = 0
    do
       if(iter==ITMAX)exit
       iter=iter+1
       if((fb > 0._r8 .and. fc > 0._r8) .or. (fb < 0._r8 .and. fc < 0._r8))then
          c=a   !Rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
       endif
       if( abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2._r8*macheps*abs(b)+0.5_r8*tol  !Convergence check.
       xm=0.5_r8*(c-b)
       if(abs(xm) <= tol1 .or. fb == 0.)then
          x=b
          return
       endif
       if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa !Attempt inverse quadratic interpolation.
          if(a == c) then
             p=2._r8*xm*s
             q=1._r8-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2._r8*xm*q*(q-r)-(b-a)*(r-1._r8))
             q=(q-1._r8)*(r-1._r8)*(s-1._r8)
          endif
          if(p > 0._r8) q=-q !Check whether in bounds.
          p=abs(p)
          if(2._r8*p < min(3._r8*xm*q-abs(tol1*q),abs(e*q))) then
             e=d !Accept interpolation.
             d=p/q
          else
             d=xm  !Interpolation failed, use bisection.
             e=d
          endif
       else !Bounds decreasing too slowly, use bisection.
          d=xm
          e=d
       endif
       a=b !Move last best guess to a.
       fa=fb
       if(abs(d) > tol1) then !Evaluate new trial root.
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif
       call func(b,fb, pp)
       if(fb==0._r8)exit
    enddo
    if(iter==ITMAX)write(iulog,*) 'brent exceeding maximum iterations', b, fb
    x=b

  end subroutine brent_p
  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: brent
  !
  ! !INTERFACE:
  
  subroutine brent_np(x, x1,x2,f1, f2, macheps, tol,func)

    !
    !!DESCRIPTION:
    !Use Brent's method to find the root to a single variable function func, which is known to exist between x1 and x2.
    !The found root will be updated until its accuracy is tol.

    !!REVISION HISTORY:
    !Dec 14/2012: Jinyun Tang, modified from numerical recipes in F90 by press et al. 1188-1189
    !
    !!USES:

    !
    !!ARGUMENTS:
    implicit none
    real(r8), intent(in) :: x1, x2, f1, f2    !minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in) :: macheps           !machine precision
    real(r8), intent(in) :: tol               !the error tolerance
    real(r8), intent(out):: x
    interface
       subroutine func(x,f)
         use shr_kind_mod       , only : r8 => shr_kind_r8
         implicit none
         real(r8), intent(in)  :: x
         real(r8), intent(out) :: f
       end subroutine func
    end interface

    ! !CALLED FROM:
    ! whenever it is needed

    integer, parameter :: ITMAX = 40            !maximum number of iterations
    integer, parameter :: iulog = 6
    integer :: iter
    real(r8)  :: a,b,c,d,e,fa,fb,fc,p,q,r,s,xm,tol1


    a=x1
    b=x2
    fa=f1
    fb=f2
    if((fa > 0._r8 .and. fb > 0._r8).or.(fa < 0._r8 .and. fb < 0._r8))then
       write(iulog,*) 'root must be bracketed for brent'
       write(iulog,*) 'a=',a,' b=',b,' fa=',fa,' fb=',fb
       call endrun(msg=errmsg(__FILE__, __LINE__))
    endif
    c=b
    fc=fb
    iter = 0
    do
       if(iter==ITMAX)exit
       iter=iter+1
       if((fb > 0._r8 .and. fc > 0._r8) .or. (fb < 0._r8 .and. fc < 0._r8))then
          c=a   !Rename a, b, c and adjust bounding interval d.
          fc=fa
          d=b-a
          e=d
       endif
       if( abs(fc) < abs(fb)) then
          a=b
          b=c
          c=a
          fa=fb
          fb=fc
          fc=fa
       endif
       tol1=2._r8*macheps*abs(b)+0.5_r8*tol  !Convergence check.
       xm=0.5_r8*(c-b)
       if(abs(xm) <= tol1 .or. fb == 0.)then
          x=b
          return
       endif
       if(abs(e) >= tol1 .and. abs(fa) > abs(fb)) then
          s=fb/fa !Attempt inverse quadratic interpolation.
          if(a == c) then
             p=2._r8*xm*s
             q=1._r8-s
          else
             q=fa/fc
             r=fb/fc
             p=s*(2._r8*xm*q*(q-r)-(b-a)*(r-1._r8))
             q=(q-1._r8)*(r-1._r8)*(s-1._r8)
          endif
          if(p > 0._r8) q=-q !Check whether in bounds.
          p=abs(p)
          if(2._r8*p < min(3._r8*xm*q-abs(tol1*q),abs(e*q))) then
             e=d !Accept interpolation.
             d=p/q
          else
             d=xm  !Interpolation failed, use bisection.
             e=d
          endif
       else !Bounds decreasing too slowly, use bisection.
          d=xm
          e=d
       endif
       a=b !Move last best guess to a.
       fa=fb
       if(abs(d) > tol1) then !Evaluate new trial root.
          b=b+d
       else
          b=b+sign(tol1,xm)
       endif
       call func(b, fb)
       if(fb==0._r8)exit
    enddo
    if(iter==ITMAX)write(iulog,*) 'brent exceeding maximum iterations', b, fb
    x=b

  end subroutine brent_np

  !------------------------------------------------------------------------------------
  subroutine hybrid_findroot_p(x0, p, iter, func)
    !
    !! DESCRIPTION:
    ! use a hybrid solver to find the root of equation
    ! f(x) = x- h(x),
    !  s.t. f(x) = 0.
    !the hybrid approach combines the strength of the newton secant approach (find the solution domain)
    !and the bisection approach implemented with the Brent's method to guarrantee convergence.
    !
    !! REVISION HISTORY:
    !Apr 14/2013: created by Jinyun Tang

    implicit none
    ! !ARGUMENTS:
    real(r8) , intent(inout) :: x0           !solution's initial guess
    integer    ,  intent(in) :: p            !index used in the function
    integer    , intent(out) :: iter         !number of used iterations
    interface
       subroutine func(x,f,p)
         use shr_kind_mod       , only : r8 => shr_kind_r8
         implicit none
         real(r8), intent(in)  :: x
         real(r8), intent(out) :: f
         integer,  intent(in)  :: p
       end subroutine func
    end interface

    ! !LOCAL VARIABLES:
    real(r8) :: a, b
    real(r8) :: fa, fb
    real(r8) :: x1, f0, f1
    real(r8) :: x, dx
    real(r8), parameter :: eps = 1.e-2_r8      !relative accuracy
    real(r8), parameter :: eps1= 1.e-4_r8
    integer,  parameter :: itmax = 40          !maximum number of iterations
    real(r8) :: tol,minx,minf


    call func(x0, f0, p)
    if(f0 == 0._r8)return

    minx=x0
    minf=f0
    x1 = x0 * 0.99_r8
    call func(x1,f1, p)

    if(f1==0._r8)then
       x0 = x1
       return
    endif
    if(f1<minf)then
       minx=x1
       minf=f1
    endif

    !first use the secant approach, then use the brent approach as a backup
    iter = 0
    do
       iter = iter + 1
       dx = - f1 * (x1-x0)/(f1-f0)
       x = x1 + dx
       tol = abs(x) * eps
       if(abs(dx)<tol)then
          x0 = x
          exit
       endif
       x0 = x1
       f0 = f1
       x1 = x
       call func(x1,f1, p)
       if(f1<minf)then
          minx=x1
          minf=f1
       endif
       if(abs(f1)<=eps1)then
          x0 = x1
          exit
       endif

       !if a root zone is found, use the brent method for a robust backup strategy
       if(f1 * f0 < 0._r8)then
          call brent_p(x, x0,x1,f0,f1, eps, tol, p, func)
          x0=x
          exit
       endif
       if(iter>itmax)then
          !in case of failing to converge within itmax iterations
          !stop at the minimum function
          !this happens because of some other issues besides the stomatal conductance calculation
          !and it happens usually in very dry places and more likely with c4 plants.
          call func(minx,f1, p)
          exit
       endif
    enddo
  end subroutine hybrid_findroot_p


  !------------------------------------------------------------------------------------
  subroutine hybrid_findroot_np(x0, iter, func)
    !
    !! DESCRIPTION:
    ! use a hybrid solver to find the root of equation
    ! f(x) = x- h(x),
    ! s.t. f(x) = 0.
    !the hybrid approach combines the strength of the newton secant approach (find the solution domain)
    !and the bisection approach implemented with the Brent's method to guarrantee convergence.

    !
    !! REVISION HISTORY:
    !Apr 14/2013: created by Jinyun Tang

    implicit none
    ! !ARGUMENTS:
    real(r8), intent(inout) :: x0           !solution's initial guess
    integer,  intent(out) :: iter
    interface
       subroutine func(x,f)
         use shr_kind_mod       , only : r8 => shr_kind_r8
         implicit none
         real(r8), intent(in)  :: x
         real(r8), intent(out) :: f
       end subroutine func
    end interface

    ! !LOCAL VARIABLES:
    real(r8) :: a, b
    real(r8) :: fa, fb
    real(r8) :: x1, f0, f1
    real(r8) :: x, dx
    real(r8), parameter :: eps = 1.e-2_r8      !relative accuracy
    real(r8), parameter :: eps1= 1.e-4_r8
    integer,  parameter :: itmax = 40          !maximum number of iterations
    real(r8) :: tol,minx,minf

    call func(x0, f0)
    if(f0 == 0._r8)return

    minx=x0
    minf=f0
    x1 = x0 * 0.99_r8
    call func(x1,f1)

    if(f1==0._r8)then
       x0 = x1
       return
    endif
    if(f1<minf)then
       minx=x1
       minf=f1
    endif

    !first use the secant approach, then use the brent approach as a backup
    iter = 0
    do
       iter = iter + 1
       dx = - f1 * (x1-x0)/(f1-f0)
       x = x1 + dx
       tol = abs(x) * eps
       if(abs(dx)<tol)then
          x0 = x
          exit
       endif
       x0 = x1
       f0 = f1
       x1 = x
       call func(x1,f1)
       if(f1<minf)then
          minx=x1
          minf=f1
       endif
       if(abs(f1)<=eps1)then
          x0 = x1
          exit
       endif

       !if a root zone is found, use the brent method for a robust backup strategy
       if(f1 * f0 < 0._r8)then
          call brent_np(x, x0,x1,f0,f1, eps, tol, func)
          x0=x
          exit
       endif
       if(iter>itmax)then
          !in case of failing to converge within itmax iterations
          !stop at the minimum function
          !this happens because of some other issues besides the stomatal conductance calculation
          !and it happens usually in very dry places and more likely with c4 plants.
          call func(minx,f1)
          exit
       endif
    enddo
  end subroutine hybrid_findroot_np

  !--------------------------------------------------------------------------
  SUBROUTINE gaussian_solve(a,b,error)
    ! !DESCRIPTION:
    ! This subroutine solves the linear system Ax = b
    !  Copyright 1994, Miles Ellis, Ivor Philips and Tom Lahey
    !  Copyright 1994, Addison-Wesley Publishers Ltd.
    !  Copyright 1994, Addison-Wesley Publishing Company Inc.
    !  Permission is granted for the use of this code for the purpose of teaching
    !  and/or learning the Fortran 90 language provided that the above copyright
    !  notices are included in any copies made.
    !  Neither the authors nor the publishers accept any responsibility for
    !  any results obtained by use of this code.
    !  modified by Jinyun Tang

    ! !ARGUMENTS:
    real(r8), dimension(:,:), intent(inout) :: a     !coefficients of A
    real(r8), dimension(:)  , intent(inout) :: b     !right handside and solution on returning.
    integer, intent(out)                    :: error ! indicates if errors are found

    ! Reduce the equations by Gaussian elimination
    call gaussian_elimination(a,b,error)

    ! If reduction was successful, calculate solution by
    ! back substitution
    if (error == 0) call back_substitution(a,b,error)

  end subroutine gaussian_solve

  !---------------------------------------------------------------------------

  subroutine gaussian_elimination(a,b,error)

    ! !DESCRIPTION:
    ! This subroutine performs Gaussian elimination on a
    ! system of linear equations
    !  Copyright 1994, Miles Ellis, Ivor Philips and Tom Lahey
    !  Copyright 1994, Addison-Wesley Publishers Ltd.
    !  Copyright 1994, Addison-Wesley Publishing Company Inc.
    !  Permission is granted for the use of this code for the purpose of teaching
    !  and/or learning the Fortran 90 language provided that the above copyright
    !  notices are included in any copies made.
    !  Neither the authors nor the publishers accept any responsibility for
    !  any results obtained by use of this code.

    ! !USES:
    use MathfuncMod, only : swap
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:,:), intent(inout) :: a !contains the coefficients
    real(r8), dimension(:)  , intent(inout) :: b !contains the right-hand side
    integer,  intent(out)                   :: error

    ! !LOCAL VARIABLES:
    real(r8), dimension(size(a,1)) :: temp_array ! Automatic array
    integer, dimension(1) :: ksave
    integer :: i, j, k, n
    real(r8) :: temp, m

    ! Validity checks
    n = size(a,1)

    if (n == 0) then
       error = -1              ! There is no problem to solve
       return
    endif

    if (n /= size(a,2))then
       error = -2              ! a is not square
       return
    endif

    if (n/=size(b))then
       error = -3              ! Size of b does not match a
       return
    endif

    ! Dimensions of arrays are OK, so go ahead with Gaussian
    ! elimination
    error = 0

    do i = 1, n-1
       ! Find row with largest value of |a(j,i)|, j=i, ..., n
       ksave = maxloc(abs(a(i:n, i)))

       ! Check whether largest |a(j,i)| is near zero
       k = ksave(1) + i - 1

       if ( abs(a(k,i)) <= 1.e-5_r8 ) then
          error = -4            ! No solution possible
          return
       endif

       !Interchange row i and row k, if necessary
       if(k /= i) then
          call swap(a(i,:),a(k,:))
          ! Interchange corresponding elements of b
          call swap(b(i), b(k))
       endif

       ! Subtract multiples of row i from subsequent rows to
       ! zero all subsequent coefficients of x sub i
       do j = i + 1, n
          m = a(j,i)/a(i,i)
          a(j,:) = a(j,:) - m*a(i,:)
          b(j) = b(j) - m*b(i)
       enddo
    enddo

  end subroutine gaussian_elimination


  !-----------------------------------------------------------------
  SUBROUTINE back_substitution(a,b,error)

    ! !DESCRIPTION:
    !
    ! This subroutine performs back substition once a system
    ! of equations has been reduced by Gaussian elimination
    !  Copyright 1994, Miles Ellis, Ivor Philips and Tom Lahey
    !  Copyright 1994, Addison-Wesley Publishers Ltd.
    !  Copyright 1994, Addison-Wesley Publishing Company Inc.
    !  Permission is granted for the use of this code for the purpose of teaching
    !  and/or learning the Fortran 90 language provided that the above copyright
    !  notices are included in any copies made.
    !  Neither the authors nor the publishers accept any responsibility for
    !  any results obtained by use of this code.
    !  Modified by Jinyun Tang, Apr, 2013

    implicit none

    ! !ARGUMENTS:
    real(r8), dimension(:,:), intent(in)  :: a     !contains the coefficients
    real(r8), dimension(:), intent(inout) :: b     !contains the right-hand side coefficients, will contain the solution on exit
    integer , intent(out)                 :: error ! will be set non-zero if an error is found

    ! !LOCAL VARIABLES:
    real(r8) :: sum
    integer :: i,j,n


    error = 0
    n = size(b)
    ! Solve for each variable in turn
    do i = n,1,-1
       ! Check for zero coefficient
       if ( abs(a(i,i)) <= 1.e-5_r8 ) then
          error = -4
          return
       endif
       sum = b(i)
       do j = i+1,n
          sum = sum - a(i,j)*b(j)
       enddo
       b(i) = sum/a(i,i)
    enddo

  end subroutine back_substitution
end module FindRootMod
