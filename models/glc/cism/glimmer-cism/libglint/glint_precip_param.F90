!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_precip_param.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glint_precip_param

  !*FD The Roe/Lindzen precip downscaling parameterization

  ! Author: Ian Rutt
  ! Date:   19/11/03

  use glimmer_global

  implicit none

  private satvap,calc_w0,error_func

contains

  subroutine glint_precip(precip,xwind,ywind,temp,topo,dx,dy,fixed_a)

    !*FD Calculates the precipitation field over 
    !*FD the ice sheet using the parameterization given in 
    !*FD Roe (2002)\footnote{{\em J. Glaciol.} {\bf 48,} no.160 pp.70--80}
    !*FD Note that:
    !*FD \begin{itemize}
    !*FD \item All arrays are on the ice model grid, and must have the same shape.
    !*FD \item There is some confusion in the lit between Roe (2002) and Roe and Lindzen (2001)
    !*FD       over the dimensions of some quantities. I have used the combination
    !*FD       that is most consistent.
    !*FD \item For the value of a, which R\&L take to be $2.5\times 10^{-11}
    !*FD        \mathrm{m}^2\,\mathrm{s\,kg}^{-1}$, I use
    !*FD       the precip rate divided by the saturated vapour pressure, unless 
    !*FD       \texttt{fixed\_a} is set.
    !*FD \item The equation used is:
    !*FD \[P=be_{\mathrm{sat}}(T_\mathrm{s})\left[\frac{|x_0|}{2}+\frac{|x_0|}{2}
    !*FD \mathrm{erf}(|x_0|/\alpha)+\frac{\alpha}{2\sqrt{\pi}}\exp(-(1/\alpha)^2x_0^2)\right]\]
    !*FD with $x_0=\frac{a}{b}-w_0$,
    !*FD $w_0$ is the mean vertical velocity
    !*FD $b=5.9\times 10^{-9}\,\mathrm{kg\,m\,s^2}$,
    !*FD $\alpha=0.0115\,\mathrm{ms^{-1}}$,
    !*FD and either $a=2.5\times 10^{-11}\,\mathrm{m^2\,s\,kg^{-1}}$, or
    !*FD $a=P/e_{\mathrm{sat}}(T)$, depending on the value of \texttt{fixed\_a}.
    !*FD \end{itemize}

    implicit none

    ! Arguments

    real(sp),dimension(:,:),intent(inout) :: precip  !*FD Precipitation field (mm/a) 
                                                     !*FD used for input and output. on
                                                     !*FD input it contains the large-scale 
                                                     !*FD field calculated by interpolation. 
                                                     !*FD On output, it contains the field calculated by 
                                                     !*FD this subroutine and used for the mass-balance.
    real(rk), dimension(:,:),intent(in)    :: xwind   !*FD Annual mean wind field: $x$-component (m/s)
    real(rk), dimension(:,:),intent(in)    :: ywind   !*FD Annual mean wind field: $y$-component (m/s)
    real(sp), dimension(:,:),intent(in)    :: topo    !*FD Surface topography (m)
    real(sp),dimension(:,:), intent(in)    :: temp    !*FD Mean annual surface temperature field,
                                                      !*FD corrected for height ($^{\circ}$C)
    real(rk),                intent(in)    :: dx      !*FD $x$ grid spacing (m)
    real(rk),                intent(in)    :: dy      !*FD $y$ grid spacing (m)
    logical,optional,        intent(in)    :: fixed_a !*FD Set to fix $\mathtt{a}=2.5\times 10^{-11} 
                                                      !*FD \mathrm{m}^2\,\mathrm{s\,kg}^{-1}$ over the 
                                                      !*FD whole domain, else scale \texttt{a} as described below. 
                                                      !*FD If not present, assumed \texttt{.false.}.

    ! Internal variables

    integer :: nx,ny,i,j
    real(rk),dimension(size(precip,1),size(precip,2)) :: w0
    real(rk),parameter :: pi=3.141592654
    real(rk),parameter :: b=5.9e-9       ! m s^2 kg (these are Roe and Lindzens dims) 
    real(rk),parameter :: alpha=0.0115   ! ms^-1
    real(rk),parameter :: pc=3.17098e-11 ! precipitation conversion factor: mm/a -> m/s
    real(rk)           :: x0,a
    logical :: fa                        ! ever-present proxy for fixed_a

    ! Beginning of code

    nx=size(precip,1) ; ny=size(precip,2)

    if (present(fixed_a)) then
       fa=fixed_a
    else
       fa=.false.
    endif

    w0=calc_w0(xwind,ywind,topo,dx,dy)   ! calculate the mean vertical velocity

    do i=1,nx
       do j=1,ny
          if (fa) then
             a=2.5e-11                             ! fixed a
          else
             a=precip(i,j)*pc/satvap(real(temp(i,j),rk))    ! calculate a from the precip and satvap
          endif
          x0=a/b+w0(i,j)                        ! for convenience, calc x0
          precip(i,j)=b*satvap(real(temp(i,j),rk))* &    ! the function from Roe and Lindzen (corrected)
               (0.5*abs(x0)+(x0**2/2*abs(x0))*error_func((1.0/alpha)*abs(x0)) &
               +(alpha/(2*sqrt(pi)))*exp(-(1.0/alpha)**2*x0**2))
          a=a
       enddo
    enddo

    precip=precip/pc   ! convert back to mm/a 

  end subroutine glint_precip

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Private subroutines and functions follow
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function satvap(temp,kelvin)

    !*FD Calculates the saturated vapour pressure using the 
    !*FD Clausius-Clapyron equation (from Roe 2002)
    !*FD Note that by default the units of temperature are degC,
    !*FD unless the kelvin flag is present and set, in which
    !*FD case they are Kelvin.
    !*RV Saturated vapour pressure in Pascals.

    implicit none

    ! Arguments

    real(rk),        intent(in) :: temp    !*FD Temperature ($^{\circ}$C or K)
    logical,optional,intent(in) :: kelvin  !*FD Set if temperature is in Kelvin

    ! Internal variables

    real(rk) :: ts
    real(rk),parameter :: e0=6.112          ! This is in millibars, but multiplied by 100 below.
    real(rk),parameter :: c1=17.67,c2=243.5 ! These names from Roe(2002)

    ! Beginning of code

    ! First check to see if kelvin is present and adjust accordingly

    if (present(kelvin)) then
       if (kelvin) then 
          ts=temp-273.15
       else
          ts=temp
       endif
    else
       ts=temp
    endif

    ! finally, the function itself

    satvap=100*e0*exp(c1*ts/(c2+ts))

  end function satvap

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function calc_w0(u,v,h0,dx,dy)

    !*FD Calculates the mean vertical velocity field, 
    !*FD based on the horizontal flow, and the topography.
    !*FD 
    !*FD Note that:
    !*FD \begin{itemize}
    !*FD \item All input arrays must be of the same rank and size.
    !*FD \item The vertical velocity is calculated from the divergence
    !*FD       of the horizontal wind over the topography.
    !*FD       \[w_0=u\frac{\partial h_0}{\partial x}+
    !*FD       v\frac{\partial h_0}{\partial y} \]
    !*FD \item Differentiation is done with conventional centred-differences,
    !*FD       except at the corners, where uncentred differences are employed.
    !*FD \end{itemize}

    !*RV An array of the same size as \texttt{u} is 
    !*RV returned. The units are m/s.

    implicit none

    ! Arguments

    real(rk),dimension(:,:),intent(in) :: u   !*FD The $x$ component of the mean wind field (m/s)
    real(rk),dimension(:,:),intent(in) :: v   !*FD The $y$ component of the mean wind field (m/s)
    real(sp),dimension(:,:),intent(in) :: h0  !*FD The topography (m)
    real(rk),               intent(in) :: dx  !*FD The $x$ grid-spacing (m)
    real(rk),               intent(in) :: dy  !*FD The $y$ grid-spacing (m)

    ! Returned array

    real(rk),dimension(size(u,1),size(u,2)) :: calc_w0

    ! Internal variables

    integer :: i,j,nx,ny

    ! Beginning of code

    nx=size(u,1) ; ny=size(u,2)

    ! main block loop

    do i=2,nx-1
       do j=2,ny-1
          calc_w0(i,j)=u(i,j)*(h0(i+1,j)-h0(i-1,j))/(2*dx) +v(i,j)*(h0(i,j+1)-h0(i,j-1))/(2*dy)
       enddo
    enddo

    ! top and bottom rows

    do i=2,nx-1
       calc_w0(i,1) =u(i,1) *(h0(i+1,1) -h0(i-1,1)) /(2*dx) +v(i,1) *(h0(i,2) -h0(i,1))   /dy
       calc_w0(i,ny)=u(i,ny)*(h0(i+1,ny)-h0(i-1,ny))/(2*dx) +v(i,ny)*(h0(i,ny)-h0(i,ny-1))/dy
    enddo

    ! left and right columns

    do j=2,ny-1
       calc_w0(1,j) =u(1,j) *(h0(2,j) -h0(1,j))   /dx +v(1,j) *(h0(1,j+1) -h0(1,j-1)) /(2*dy)
       calc_w0(nx,j)=u(nx,j)*(h0(nx,j)-h0(nx-1,j))/dx +v(nx,j)*(h0(nx,j+1)-h0(nx,j-1))/(2*dy)
    enddo

    ! corners

    calc_w0(1,1)  =u(1,1)  *(h0(2,1)  -h0(1,1))    /dx +v(1,1)  *(h0(1,2)  -h0(1,1))    /dy
    calc_w0(1,ny) =u(1,ny) *(h0(2,ny) -h0(1,ny))   /dx +v(1,ny) *(h0(1,ny) -h0(1,ny-1)) /dy
    calc_w0(nx,1) =u(nx,1) *(h0(nx,1) -h0(nx-1,1)) /dx +v(nx,1) *(h0(nx,2) -h0(nx,1))   /dy
    calc_w0(nx,ny)=u(nx,ny)*(h0(nx,ny)-h0(nx-1,ny))/dx +v(nx,ny)*(h0(nx,ny)-h0(nx,ny-1))/dy

  end function calc_w0

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(rk) function error_func(y)

    !*FD The error function
    !*FD
    !*FD The error function may be approximated by:
    !*FD \[ \mathrm{erf}(y)=1-(\gamma_1 t+\gamma_2 t^2+\gamma_3 t^3)\exp(-y^2)\]
    !*FD with
    !*FD \[t=\frac{1}{1+\gamma_0 y}\]  
    !*FD and
    !*FD \[\gamma_0=0.47047, \]
    !*FD \[\gamma_1=0.3480242, \]
    !*FD \[\gamma_2=-0.0958798, \]
    !*FD \[\gamma_3=0.7478556 \]
    !*FD (from Abramowitz and Stegun 1965)
    !*FD However, this doesn't seem to be right for $\mathtt{y}<0$, but ok for $\mathtt{y}\geq 0$.
    !*FD Since the input is always $>0$, this isn't a problem here.
    !*RV The value of the error function at \texttt{y}.

    implicit none

    real(rk),intent(in) :: y !*FD The independent variable

    real(rk) :: t
    real(rk),parameter :: g0=0.47047
    real(rk),parameter :: g1=0.3480242
    real(rk),parameter :: g2=-0.0958798
    real(rk),parameter :: g3=0.7478556

    t=1/(1+g0*y)

    error_func=1-(g1*t+g2*t**2+g3*t**3)*exp(-y**2)

  end function error_func

end module glint_precip_param
