module FindRootMod
  !
  ! !DESCRIPTION:
  !  Functions to solve simple equations
  !  History: created by Jinyun Tang, 2013

  ! !USES:
  use bshr_kind_mod , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
  use betr_ctrl     , only : iulog  => biulog, do_betr_output
  use MathfuncMod   , only : is_bounded

  implicit none

  character(len=*), parameter :: mod_filename = &
       __FILE__

contains
  !-------------------------------------------------------------------------------
  subroutine quadrootbnd(a,b,c, xl, xr, x, bstatus)
    !
    ! !DESCRIPTION:
    ! return a root of the qudratic equation
    ! within bound xl and xr
    use betr_constants  , only : betr_errmsg_len
    use BetrStatusType  , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: a, b, c
    real(r8), intent(in) :: xl, xr
    real(r8), intent(out):: x
    type(betr_status_type), intent(out)   :: bstatus

    ! !LOCAL VARIABLES:
    real(r8) :: delta
    character(len=betr_errmsg_len) :: msg
    character(len=32) :: subname ='quadrootbnd'

    call bstatus%reset()
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
             msg = 'no bounded solution for the given quadratic equation '//errmsg(mod_filename, __LINE__)
             call bstatus%set_msg(msg=msg, err=-1)
          endif
       endif
    else
       msg='no real solution for the given quadratic equation '//errmsg(mod_filename, __LINE__)
       call bstatus%set_msg(msg=msg, err=-1)
    endif
    return
  end subroutine quadrootbnd

  !-------------------------------------------------------------------------------

  subroutine quadproot(a,b,c,x, bstatus)
    !
    ! !DESCRIPTION:
    ! return positive root of the qudratic equation
    use BetrStatusType  , only : betr_status_type
    use betr_constants  , only : betr_errmsg_len
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: a, b, c
    real(r8), intent(out):: x
    type(betr_status_type), intent(out) :: bstatus

    ! !LOCAL VARIABLES:
    real(r8) :: delta
    character(len=32) :: subname ='quadproot'
    character(len=betr_errmsg_len) :: msg
    call bstatus%reset()
    delta = b * b -4._r8 * a * c
    if(delta>=0._r8)then
       x = (-b + sqrt(delta))/2._r8
    else
       msg='no positive solution for the given quadratic equation '//errmsg(mod_filename, __LINE__)
       call bstatus%set_msg(msg=msg, err=-1)
    endif
    return
  end subroutine quadproot
  !===============================================================================

  subroutine cubicrootbnd(a,b,c,d, xl, xr, x, bstatus)
    !
    ! !DESCRIPTION:
    ! return positive root of the cubic equation
    !
    ! !USES:
    use betr_varcon     , only : rpi  => brpi
    use BetrStatusType  , only : betr_status_type
    use betr_constants  , only : betr_errmsg_len
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: a, b, c, d
    real(r8), intent(in) :: xl, xr
    real(r8), intent(out):: x
    type(betr_status_type), intent(out) :: bstatus

    ! !LOCAL VARIABLES:
    real(r8) :: p, q
    real(r8) :: b1, c1, d1
    real(r8) :: n, u, f, y
    real(r8) :: delta
    character(len=betr_errmsg_len) :: msg
    character(len=32) :: subname ='cubicrootbnd'

    call bstatus%reset()
    b1 = b/a
    c1 = c/a
    d1 = d/a

    p = c1 - b1 * b1 /3._r8
    q = d1 - b1 / 3._r8 * (c1 - 2._r8 * b1**2._r8/9._r8)

    delta =-4._r8 * p**3._r8 - 27._r8 * q ** 2._r8
    if(delta<0._r8)then
       msg='no real solution for the given cubic equation '//errmsg(mod_filename, __LINE__)
       call bstatus%set_msg(msg=msg,err=-1)
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
                msg='no bounded solution for the given cubic equation '//errmsg(mod_filename, __LINE__)
                call bstatus%set_msg(msg=msg,err=-1)
             endif
          endif
       endif
    endif
  end subroutine cubicrootbnd

  !===============================================================================
  subroutine cubicproot(a,b,c,d, x, bstatus)
    !
    ! !DESCRIPTION:
    ! return positive root of the cubic equation
    !
    ! !USES:
    use betr_varcon    , only : rpi  => brpi
    use BetrStatusType , only : betr_status_type
    use betr_constants , only : betr_errmsg_len
    implicit none
    ! !ARGUMENTS:
    real(r8), intent(in) :: a, b, c, d
    real(r8), intent(out) :: x
    type(betr_status_type), intent(out) :: bstatus

    ! !LOCAL VARIABLES:
    real(r8) :: p, q
    real(r8) :: b1, c1, d1
    real(r8) :: n, u, f, y
    real(r8) :: delta
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    b1 = b/a
    c1 = c/a
    d1 = d/a

    p = c1 - b1 * b1 /3._r8
    q = d1 - b1 / 3._r8 * (c1 - 2._r8 * b1**2._r8/9._r8)
    delta =-4._r8 * p**3._r8 - 27._r8 * q ** 2._r8
    if(delta<0._r8)then
       msg='no real solution for the given cubic equation '//errmsg(mod_filename, __LINE__)
       call bstatus%set_msg(msg=msg,err=-1)
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
  end subroutine cubicproot

  !===============================================================================

  subroutine LUsolvAxr(a,r, n, bstatus)
    ! !DESCRIPTION:
    !solve linear equation Ax=r, using the LU decomposition
    use BetrStatusType    , only : betr_status_type
    implicit none
    ! !ARGUMENTS:
    integer  , intent(in)    :: n
    real(r8) , intent(inout) :: a(n,n)
    real(r8) , intent(inout) :: r(n)
    type(betr_status_type), intent(out)   :: bstatus

    ! !LOCAL VARIABLES:
    real(r8) :: d(n)
    integer :: indx(n)

    call bstatus%reset()
    !do lu decomposition
    call ludcmp(a,indx,d,n, bstatus)
    if(bstatus%check_status())return
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
    integer  ,    intent(in) :: n
    real(r8) ,    intent(in) :: a(n,n)
    integer  ,    intent(in) :: indx(n)
    real(r8) , intent(inout) :: b(n)

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
       else if (summ /= 0.0_r8) then
          ii=i !A nonzero element was encountered, so from now on we will
       end if !have to do the dot product above.
       b(i)=summ
    end do
    do i=n,1,-1 !Now we do the backsubstitution, equation (2.3.7).
       b(i) = (b(i)-dot_product(a(i,i+1:n),b(i+1:n)))/a(i,i)
    end do
  end subroutine lubksb
  !===============================================================================


  subroutine ludcmp(a,indx,d,n, bstatus)
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
    use MathfuncMod       , only : swap
    use BetrStatusType    , only : betr_status_type
    use betr_constants    , only : betr_errmsg_len
    implicit none
    ! !ARGUMENTS:
    integer  ,  intent(in)  :: n
    real(r8), intent(inout) :: a(n,n)
    integer   , intent(out) :: indx(n)
    real(r8)  , intent(out) :: d(n)
    type(betr_status_type), intent(out)   :: bstatus

    ! !LOCAL VARIABLES:
    real(r8), dimension(size(a,1)) :: vv !vv stores the implicit scaling of each row.
    real(r8), parameter :: TINY=1.0e-20  !A small number.
    integer :: j,imax
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    d=1.0_r8 !No row interchanges yet.
    vv=maxval(abs(a),dim=2) !Loop over rows to get the implicit scaling
    if (any(vv == 0.0_r8)) then
       write(msg,*)'singular matrix in ludcmp'//errMsg(mod_filename, __LINE__) !information.
       call bstatus%set_msg(msg=msg,err=-1)
       return
    endif
    !There is a row of zeros.
    vv=1.0_r8 / vv !Save the scaling.
    do j=1,n
       imax=(j-1)+imaxloc(vv(j:n)*abs(a(j:n,j))) !Find the pivot row.
       if (j /= imax) then !Do we need to interchange rows?
          call swap(a(imax,:),a(j,:), bstatus) !Yes, do so...
          d=-d !...and change the parity of d.
          vv(imax)=vv(j) !Also interchange the scale factor.
       end if
       indx(j)=imax
       if (a(j,j) == 0.0_r8) a(j,j)=TINY
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

  !------------------------------------------------------------------------------
  !BOP
  !
  ! !IROUTINE: brent
  !
  ! !INTERFACE:

  subroutine brent(x, x1,x2,f1, f2, macheps, tol, func_data, func, bstatus)

    !
    !!DESCRIPTION:
    !Use Brent's method to find the root to a single variable function func, which is known to exist between x1 and x2.
    !The found root will be updated until its accuracy is tol.

    !!REVISION HISTORY:
    !Dec 14/2012: Jinyun Tang, modified from numerical recipes in F90 by press et al. 1188-1189
    !
    !!USES:
    use func_data_type_mod, only : func_data_type
    use betr_constants    , only : betr_errmsg_len
    use BetrstatusType    , only : betr_status_type
    !
    !!ARGUMENTS:
    implicit none
    real(r8), intent(in) :: x1, x2, f1, f2    !minimum and maximum of the variable domain to search for the solution ci_func(x1) = f1, ci_func(x2)=f2
    real(r8), intent(in) :: macheps           !machine precision
    real(r8), intent(in) :: tol               !the error tolerance
    type(func_data_type), intent(in) :: func_data ! data passed to subroutine func
    real(r8), intent(out):: x
    type(betr_status_type), intent(out) :: bstatus

    interface
       subroutine func(x, func_data, f)
         use bshr_kind_mod        , only : r8 => shr_kind_r8
         use func_data_type_mod, only : func_data_type
         implicit none
         real(r8), intent(in)  :: x
         type(func_data_type), intent(in) :: func_data ! data passed to subroutine func
         real(r8), intent(out) :: f
       end subroutine func
    end interface

    ! !CALLED FROM:
    ! whenever it is needed
    integer, parameter :: ITMAX = 40            !maximum number of iterations
    integer, parameter :: iulog = 6
    integer :: iter
    real(r8)  :: a,b,c,d,e,fa,fb,fc,p,q,r,s,xm,tol1
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    a=x1
    b=x2
    fa=f1
    fb=f2
    if((fa > 0._r8 .and. fb > 0._r8).or.(fa < 0._r8 .and. fb < 0._r8) .and. do_betr_output)then
       write(msg,*) 'root must be bracketed for brent', new_line('A')//'a=',a,' b=',b,' fa=',fa,' fb=',fb
       call bstatus%set_msg(msg=msg, err=-1)
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
       if(abs(xm) <= tol1 .or. fb == 0._r8)then
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
       call func(b, func_data, fb)
       if(fb==0._r8)exit
    enddo
    if(iter==ITMAX)then
      write(msg,*) 'brent exceeding maximum iterations', b, fb
      call bstatus%set_msg(msg=msg, err=1)
    endif
    x=b

  end subroutine brent

  !------------------------------------------------------------------------------------
  subroutine hybrid_findroot(x0, iter, func, bstatus)
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
    use func_data_type_mod, only : func_data_type
    use BetrStatusType    , only : betr_status_type
    implicit none

    ! !ARGUMENTS:
    real(r8), intent(inout) :: x0           !solution's initial guess
    integer,  intent(out) :: iter
    type(betr_status_type), intent(out) :: bstatus
    interface
       subroutine func(x, func_data, f)
         use bshr_kind_mod        , only : r8 => shr_kind_r8
         use func_data_type_mod, only : func_data_type
         implicit none
         real(r8), intent(in)  :: x
         type(func_data_type), intent(in) :: func_data ! data passed to subroutine func
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

    type(func_data_type) :: func_data ! dummy data passed to subroutine func

    call bstatus%reset()
    call func(x0, func_data, f0)
    if(f0 == 0._r8)return

    minx=x0
    minf=f0
    x1 = x0 * 0.99_r8
    call func(x1, func_data, f1)

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
       call func(x1, func_data, f1)
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
          call brent(x, x0,x1,f0,f1, eps, tol, func_data, func, bstatus)
          if(bstatus%check_status())return
          x0=x
          exit
       endif
       if(iter>itmax)then
          !in case of failing to converge within itmax iterations
          !stop at the minimum function
          !this happens because of some other issues besides the stomatal conductance calculation
          !and it happens usually in very dry places and more likely with c4 plants.
          call func(minx, func_data, f1)
          exit
       endif
    enddo
  end subroutine hybrid_findroot

  !--------------------------------------------------------------------------
  SUBROUTINE gaussian_solve(a,b,bstatus)
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
    use BetrStatusType   , only : betr_status_type
    ! !ARGUMENTS:
    real(r8), dimension(:,:), intent(inout) :: a     !coefficients of A
    real(r8), dimension(:)  , intent(inout) :: b     !right handside and solution on returning.
    type(betr_status_type)  , intent(out)   :: bstatus

    call bstatus%reset()
    ! Reduce the equations by Gaussian elimination
    call gaussian_elimination(a,b,bstatus)

    ! If reduction was successful, calculate solution by
    ! back substitution
    if (.not. bstatus%check_status()) call back_substitution(a,b,bstatus)

  end subroutine gaussian_solve

  !---------------------------------------------------------------------------

  subroutine gaussian_elimination(a,b,bstatus)

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
    ! re-written by Jinyun Tang, Apr, 2013.
    ! !USES:
    use MathfuncMod      , only : swap
    use BetrStatusType   , only : betr_status_type
    use betr_constants   , only : betr_errmsg_len
    implicit none
    ! !ARGUMENTS:
    real(r8), dimension(:,:), intent(inout) :: a !contains the coefficients
    real(r8), dimension(:)  , intent(inout) :: b !contains the right-hand side
    type(betr_status_type)  , intent(out)   :: bstatus

    ! !LOCAL VARIABLES:
    real(r8), dimension(size(a,1)) :: temp_array ! Automatic array
    integer, dimension(1) :: ksave
    integer :: i, j, k, n
    real(r8) :: temp, m
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    ! Validity checks
    n = size(a,1)

    if (n == 0) then
       msg = 'There is no problem to solve'//errMsg(mod_filename,__LINE__)
       call bstatus%set_msg(msg=msg, err= -1)
       return
    endif

    if (n /= size(a,2))then
       msg = 'a is not square'//errMsg(mod_filename,__LINE__)
       call bstatus%set_msg(msg=msg, err= -2)
       return
    endif

    if (n/=size(b))then
       msg = ' Size of b does not match a'//errMsg(mod_filename,__LINE__)
       call bstatus%set_msg(msg=msg, err= -3)
       return
    endif

    ! Dimensions of arrays are OK, so go ahead with Gaussian
    ! elimination

    do i = 1, n-1
       ! Find row with largest value of |a(j,i)|, j=i, ..., n
       ksave = maxloc(abs(a(i:n, i)))

       ! Check whether largest |a(j,i)| is near zero
       k = ksave(1) + i - 1

       if ( abs(a(k,i)) <= 1.e-5_r8 ) then
          msg = 'No solution possible'//errMsg(mod_filename,__LINE__)
          call bstatus%set_msg(msg=msg, err= -4)
          return
       endif

       !Interchange row i and row k, if necessary
       if(k /= i) then
          call swap(a(i,:),a(k,:), bstatus)
          if(bstatus%check_status())return
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
  SUBROUTINE back_substitution(a,b,bstatus)

    ! !DESCRIPTION:
    !
    ! This subroutine performs back substitution once a system
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
    use BetrStatusType   , only : betr_status_type
    use betr_constants   , only : betr_errmsg_len
    implicit none

    ! !ARGUMENTS:
    real(r8), dimension(:,:), intent(in)  :: a     !contains the coefficients
    real(r8), dimension(:), intent(inout) :: b     !contains the right-hand side coefficients, will contain the solution on exit
    type(betr_status_type)  , intent(out) :: bstatus

    ! !LOCAL VARIABLES:
    real(r8) :: sum
    integer :: i,j,n
    character(len=betr_errmsg_len) :: msg

    call bstatus%reset()
    n = size(b)
    ! Solve for each variable in turn
    do i = n,1,-1
       ! Check for zero coefficient
       if ( abs(a(i,i)) <= 1.e-5_r8 ) then
          call bstatus%set_msg(msg='zero coefficient'//errMsg(mod_filename,__LINE__), err=-4)
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
