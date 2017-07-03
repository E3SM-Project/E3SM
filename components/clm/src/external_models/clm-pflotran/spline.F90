  module spline_module

#include "petsc/finclude/petscsys.h"
  
  public
  
  contains

! ************************************************************************** !

subroutine spline(x,y,n,y2)


!----------------------description-------------------------------------!
!
!     cubic spline second derivative.
!
!     press, w.h., b.p. flannery, s.a. teukolsky, and w.t. vetterling.
!     1986.  numerical recipes, the art of scientific computing,
!     cambridge university press, cambridge.  pp. 86-89.

  use PFLOTRAN_Constants_module

      implicit none
      
      PetscInt :: i,n,k
      PetscReal :: x(n),y(n),y2(n),u(n)
      PetscReal :: sig,p,qn,un

      y2(1) = 0.d0
      u(1) = 0.d0
      do i = 2,n-1
        sig = (x(i)-x(i-1))/(x(i+1)-x(i-1))
        p = sig*y2(i-1)+2.d0
        y2(i) = (sig-1.d0)/p
        u(i) = (6.d0*((y(i+1)-y(i))/(x(i+1)-x(i)) - &
        (y(i)-y(i-1))/(x(i)-x(i-1)))/ &
        (x(i+1)-x(i-1)) - sig*u(i-1))/p
      enddo
      qn = 0.d0
      un = 0.d0
      y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.d0)
      do k = n-1,1,-1
        y2(k) = y2(k)*y2(k+1)+u(k)
      enddo

      return
      end subroutine spline

! ************************************************************************** !

subroutine splint(xa,ya,y2a,n,x,y)

!     cubic spline interpolation.

!     press, w.h., b.p. flannery, s.a. teukolsky, and w.t. vetterling.
!     1986.  numerical recipes, the art of scientific computing,
!     cambridge university press, cambridge.  pp. 86-89.


      implicit none
      PetscInt :: n,k,klo,khi
      PetscReal :: xa(n),ya(n),y2a(n)
      PetscReal :: h,a,b,x,y

      klo = 1
      khi = n
   10 continue
      if (khi-klo.gt.1) then
        k = (khi+klo)/2
        if (xa(k).gt.x) then
          khi = k
        else
          klo = k
        endif
        goto 10
      endif
      h = xa(khi)-xa(klo)
      a = (xa(khi)-x)/h
      b = (x-xa(klo))/h
      y = a*ya(klo)+b*ya(khi)+ &
          ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.d0
!     y1a = (ya(khi)-ya(klo))/h - &
!           ((3.d+0*(a**2)-1.d+0)/6.d+0)*h*y2a(klo) + &
!           ((3.d+0*(b**2)-1.d+0)/6.d+0)*h*y2a(khi)

      return
      end subroutine splint

! ************************************************************************** !

subroutine locate(xx,n,x,j)

!     given an array xx of length n, and given a value x, returns a
!     value j such that x is between xx(j) and xx(j+1).  xx must be
!     monotonic, either increasing or decreasing.  j=0 or j=n is 
!     returned to indicate that x is out of range.
!
!     press, w.h., b.p. flannery, s.a. teukolsky, and w.t. vetterling.
!     1986.  numerical recipes, the art of scientific computing.
!     cambridge university press, cambridge.

      implicit none
      PetscInt :: jl,ju,jm,j,n
      PetscReal :: xx(n)
      PetscReal :: x
      
      jl = 0
      ju = n+1
   10 continue
      if (ju-jl.gt.1) then
        jm = (ju+jl)/2
        if ((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm))) then
          jl = jm
        else
          ju = jm
        endif
        goto 10
      endif
      j = jl

      return
      end subroutine locate

  end module spline_module
