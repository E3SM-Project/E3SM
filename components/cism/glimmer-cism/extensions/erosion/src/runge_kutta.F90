!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   runge_kutta.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2013
!   Glimmer-CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of Glimmer-CISM.
!
!   Glimmer-CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   Glimmer-CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with Glimmer-CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

! module doing RK4 integration
! adopted from num rec (page 558)

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module rk4module
  use glimmer_global, only : dp

  integer :: rkunit = 6           ! messages are written to this stream

  ! interface for function to be integrated
  interface 
     subroutine derivs(t,x,dxdt)
       use glimmer_global, only : dp
       implicit none
       real(kind=dp), intent(in) ::t
       real(kind=dp), intent(in), dimension(:)  :: x
       real(kind=dp), intent(out), dimension(:) :: dxdt
     end subroutine derivs
  end interface

  
contains
  subroutine odeint(x, t1, t2, eps, h1, hmin, nok, nbad, derivs)
    implicit none
    real(kind=dp), intent(inout), dimension(:) :: x
    real(kind=dp), intent(in) :: t1, t2
    real(kind=dp), intent(in) :: eps
    real(kind=dp), intent(in) :: h1, hmin
    integer, intent(out) :: nok,nbad
    external derivs

    ! local variables
    integer :: i
    real(kind=dp), dimension(size(x)) :: xscal, dxdt
    real(kind=dp) :: h, hdid, hnext,t
    
    ! parameters
    integer, parameter :: maxsteps=10000
    real(kind=dp), parameter :: tiny=1.e-30

    h = sign(h1,t2-t1)
    t = t1
    nok = 0
    nbad = 0
    do i=1,maxsteps ! take at most maxsteps
       call derivs(t, x, dxdt)
       xscal = abs(x)+abs(h*dxdt)+tiny
       
       ! cut down stepsize if step can overshoot end
       if ((t+h-t2)*(t+h-t1).gt.0.) h=t2-t
       
       call rkqc(x,dxdt,t,h,eps,xscal,hdid,hnext,derivs)
       if (hdid.eq.h) then
          nok = nok+1
       else
          nbad = nbad+1
       end if

       ! are we done?
       if ((t-t2)*(t2-t1).ge.0) then
          return
       end if
       if (abs(hnext).lt.hmin) then
          write(rkunit,*) 'Error (',__FILE__,__LINE__,'): stepsize smaller than minimum ', hnext
          exit
       end if
       h = hnext
    end do
    write(rkunit,*) 'Error (',__FILE__,__LINE__,'): Too many steps'
  end subroutine odeint

  subroutine rkqc(x, dxdt, t, htry, eps, xscal, hdid, hnext, derivs)
    implicit none
    real(kind=dp), intent(inout), dimension(:) :: x, dxdt    ! variable vector x
    real(kind=dp), intent(inout) :: t                        ! independant variable t
    real(kind=dp), intent(in)    :: htry                     ! step size to be tried
    real(kind=dp), intent(in)    :: eps                      ! required accuracy
    real(kind=dp), intent(in), dimension(:)    :: xscal      ! vector against which the error is scaled
    real(kind=dp), intent(out) :: hdid, hnext                ! the actual step size achieved and the estimated next step size
    external derivs

    ! local variables
    integer :: n
    real(kind=dp), dimension(size(x)) ::  xtemp, xsav, dxsav
    real(kind=dp) :: tsav
    real(kind=dp) :: h, hh
    real(kind=dp) :: errmax
    
    ! some parameters
    real(kind=dp), parameter :: pgrow = -0.2
    real(kind=dp), parameter :: pshrink = -0.25
    real(kind=dp), parameter :: fcor = 1./15.
    real(kind=dp), parameter :: safety = 0.9
    real(kind=dp), parameter :: errcon = 6.e-4

    n = size(x)
    ! saving inital values
    tsav = t
    xsav = x
    dxsav = dxdt

    ! set step size to initial trial size
    h = htry
    do
       ! take two half steps
       hh=0.5*h
       call rk4(xsav, dxsav, tsav, hh, xtemp, derivs)
       t = tsav + hh
       call derivs(t, xtemp, dxdt)
       call rk4(xtemp,dxdt,t,hh,x,derivs)
       t = tsav + h
       if (t .eq. tsav) then
          write(rkunit,*) 'Error (',__FILE__,__LINE__,'): Stepsize is not segnificant'
          stop
       end if
       ! take the large step
       call rk4(xsav, dxsav, tsav, h, xtemp, derivs)
       errmax = 0.
       ! calculate error estimates
       xtemp = x-xtemp
       errmax = maxval(abs(xtemp/xscal))
       ! and scale relative to required tolerance
       errmax = errmax/eps
       if (errmax.gt.1) then !truncation error too large, try smaller step
          h = safety * h * (errmax**pshrink)
       else ! step succeeded, compute size of next step
          hdid = h
          if (errmax.gt.errcon) then
             hnext = safety*h*(errmax**pgrow)
          else
             hnext = 4.*h
          end if
          exit
       end if
    end do
    
    ! mop up fith-order truncation error
    x = x+xtemp*fcor
  end subroutine rkqc

  subroutine rk4(x, dxdt, t, h, xout, derivs)
    implicit none
    real(kind=dp), intent(in), dimension(:) :: x, dxdt    ! variable vector x and derivatives
    real(kind=dp), intent(in) :: t                        ! independant variable t
    real(kind=dp), intent(in) :: h                        ! step size
    real(kind=dp), intent(out), dimension(:) :: xout      ! results
    external derivs
    
    !local variables
    real(kind=dp), dimension(size(x)) :: xtemp, dxtemp, dxm
    real(kind=dp) hh, h6, th

    hh=0.5*h
    h6=h/6.
    th=t+h

    ! first step
    xtemp = x+hh*dxdt

    ! second step
    call derivs(th,xtemp,dxtemp)
    xtemp = x+hh*dxtemp

    ! third step
    call derivs(th,xtemp,dxm)
    xtemp = x+h*dxm
    dxm = dxtemp + dxm

    ! fourth step
    call derivs(t+h,xtemp, dxtemp)
    
    ! accumulate increments with propper weights
    xout = x+h6*(dxdt+dxtemp+2.*dxm)
  end subroutine rk4

end module rk4module
