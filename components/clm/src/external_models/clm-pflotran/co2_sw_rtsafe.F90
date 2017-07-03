  module co2_sw_rtsafe_module

  use PFLOTRAN_Constants_module

  implicit none

#include "petsc/finclude/petscsys.h"

  contains

! ************************************************************************** !

  FUNCTION rtsafe(funcd,x1,x2,xacc)

  IMPLICIT NONE
  PetscReal, INTENT(IN) :: x1,x2,xacc
  PetscReal :: rtsafe

! INTERFACE
!   SUBROUTINE funcd(x,fval,fderiv)
!   IMPLICIT NONE
!   PetscReal, INTENT(IN) :: x
!   PetscReal, INTENT(OUT) :: fval,fderiv
!   END SUBROUTINE funcd
! END INTERFACE

  external funcd
  INTEGER, PARAMETER :: MAXIT=100
  INTEGER :: j
  PetscReal :: df,dx,dxold,f,fh,fl,temp,xh,xl

  call funcd(x1,fl,df)
  call funcd(x2,fh,df)
  if ((fl > 0.0 .and. fh > 0.0) .or. &
    (fl < 0.0 .and. fh < 0.0)) &
    print *, 'root must be bracketed in rtsafe'
  if (fl == 0.0) then
    rtsafe=x1
    RETURN
  else if (fh == 0.0) then
    rtsafe=x2
    RETURN
  else if (fl < 0.0) then
    xl=x1
    xh=x2
  else
    xh=x1
    xl=x2
  end if
  rtsafe=0.5*(x1+x2)
  dxold=dabs(x2-x1)
  dx=dxold
  call funcd(rtsafe,f,df)
  do j=1,MAXIT
    if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) >= 0.0 .or. &
      dabs(2.0*f) > dabs(dxold*df) ) then
      dxold=dx
      dx=0.5*(xh-xl)
      rtsafe=xl+dx
      if (xl == rtsafe) RETURN
    else
      dxold=dx
      dx=f/df
      temp=rtsafe
      rtsafe=rtsafe-dx
      if (temp == rtsafe) RETURN
    end if
    if (dabs(dx) < xacc) RETURN
    call funcd(rtsafe,f,df)
    if (f < 0.0) then
      xl=rtsafe
    else
      xh=rtsafe
    end if
  end do
  print *,'rtsafe: exceeded maximum iterations'
  END FUNCTION rtsafe

! ************************************************************************** !

subroutine bracket(func,x1,x2)
  
  implicit none
  
  PetscInt :: i,ifind
  PetscReal :: fac,f1,f2,x1,x2,df
  
  external func
  
  fac = 1.2d0
  call func(x1,f1,df)
  call func(x2,f2,df)
  ifind = 1
  do i = 1, 200
    if (f1*f2 < 0.d0) return
    if (dabs(f1) < dabs(f2)) then
      x1 = x1+fac*(x1-x2)
      call func(x1,f1,df)
     else
       x2 = x2+fac*(x2-x1)
       call func(x2,f2,df)
     endif
  enddo
  ifind = 0
  print *,'root bracket failed',x1,x2,f1,f2
  return
end subroutine bracket
end module co2_sw_rtsafe_module
