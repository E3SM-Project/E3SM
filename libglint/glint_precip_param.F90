!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glint_precip_param.F90 - part of the Community Ice Sheet Model (CISM)  
!                                                              
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!
!   Copyright (C) 2005-2014
!   CISM contributors - see AUTHORS file for list of contributors
!
!   This file is part of CISM.
!
!   CISM is free software: you can redistribute it and/or modify it
!   under the terms of the Lesser GNU General Public License as published
!   by the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   CISM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   Lesser GNU General Public License for more details.
!
!   You should have received a copy of the Lesser GNU General Public License
!   along with CISM. If not, see <http://www.gnu.org/licenses/>.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

module glint_precip_param

  ! The Roe/Lindzen precip downscaling parameterization

  ! Author: Ian Rutt
  ! Date:   19/11/03

  use glimmer_global, only: dp

  implicit none

  private satvap,calc_w0,error_func

contains

  subroutine glint_precip(precip,xwind,ywind,temp,topo,dx,dy,fixed_a)

    ! Calculates the precipitation field over 
    ! the ice sheet using the parameterization given in 
    ! Roe (2002)\footnote{{\em J. Glaciol.} {\bf 48,} no.160 pp.70--80}
    ! Note that:
    ! \begin{itemize}
    ! \item All arrays are on the ice model grid, and must have the same shape.
    ! \item There is some confusion in the lit between Roe (2002) and Roe and Lindzen (2001)
    !       over the dimensions of some quantities. I have used the combination
    !       that is most consistent.
    ! \item For the value of a, which R\&L take to be $2.5\times 10^{-11}
    !        \mathrm{m}^2\,\mathrm{s\,kg}^{-1}$, I use
    !       the precip rate divided by the saturated vapour pressure, unless 
    !       \texttt{fixed\_a} is set.
    ! \item The equation used is:
    ! \[P=be_{\mathrm{sat}}(T_\mathrm{s})\left[\frac{|x_0|}{2}+\frac{|x_0|}{2}
    ! \mathrm{erf}(|x_0|/\alpha)+\frac{\alpha}{2\sqrt{\pi}}\exp(-(1/\alpha)^2x_0^2)\right]\]
    ! with $x_0=\frac{a}{b}-w_0$,
    ! $w_0$ is the mean vertical velocity
    ! $b=5.9\times 10^{-9}\,\mathrm{kg\,m\,s^2}$,
    ! $\alpha=0.0115\,\mathrm{ms^{-1}}$,
    ! and either $a=2.5\times 10^{-11}\,\mathrm{m^2\,s\,kg^{-1}}$, or
    ! $a=P/e_{\mathrm{sat}}(T)$, depending on the value of \texttt{fixed\_a}.
    ! \end{itemize}

    use glimmer_physcon, only: pi
    implicit none

    ! Arguments

    real(dp),dimension(:,:),intent(inout) :: precip  ! Precipitation field (mm/a) 
                                                     ! used for input and output. on
                                                     ! input it contains the large-scale 
                                                     ! field calculated by interpolation. 
                                                     ! On output, it contains the field calculated by 
                                                     ! this subroutine and used for the mass-balance.
    real(dp), dimension(:,:),intent(in)    :: xwind   ! Annual mean wind field: $x$-component (m/s)
    real(dp), dimension(:,:),intent(in)    :: ywind   ! Annual mean wind field: $y$-component (m/s)
    real(dp), dimension(:,:),intent(in)    :: topo    ! Surface topography (m)
    real(dp),dimension(:,:), intent(in)    :: temp    ! Mean annual surface temperature field,
                                                      ! corrected for height ($^{\circ}$C)
    real(dp),                intent(in)    :: dx      ! $x$ grid spacing (m)
    real(dp),                intent(in)    :: dy      ! $y$ grid spacing (m)
    logical,optional,        intent(in)    :: fixed_a ! Set to fix $\mathtt{a}=2.5\times 10^{-11} 
                                                      ! \mathrm{m}^2\,\mathrm{s\,kg}^{-1}$ over the 
                                                      ! whole domain, else scale \texttt{a} as described below. 
                                                      ! If not present, assumed \texttt{.false.}.

    ! Internal variables

    integer :: nx,ny,i,j
    real(dp),dimension(size(precip,1),size(precip,2)) :: w0
!!    real(dp),parameter :: pi=3.141592654   ! use value from glimmer_physcon
    real(dp),parameter :: b = 5.9d-9         ! m s^2 kg (these are Roe and Lindzens dims) 
    real(dp),parameter :: alpha = 0.0115d0   ! ms^-1
    real(dp),parameter :: pc = 3.17098d-11   ! precipitation conversion factor: mm/a -> m/s
    real(dp)           :: x0,a
    logical :: fa                        ! ever-present proxy for fixed_a

    ! Beginning of code

    nx=size(precip,1) ; ny=size(precip,2)

    if (present(fixed_a)) then
       fa = fixed_a
    else
       fa = .false.
    endif

    w0 = calc_w0(xwind,ywind,topo,dx,dy)   ! calculate the mean vertical velocity

    do i=1,nx
       do j=1,ny
          if (fa) then
             a = 2.5d-11                             ! fixed a
          else
             a = precip(i,j)*pc/satvap(real(temp(i,j),dp))    ! calculate a from the precip and satvap
          endif
          x0 = a/b + w0(i,j)                               ! for convenience, calc x0
          precip(i,j) = b*satvap(real(temp(i,j),dp))* &    ! the function from Roe and Lindzen (corrected)
               (0.5d0*abs(x0)+(x0**2/2.d0*abs(x0))*error_func((1.d0/alpha)*abs(x0)) &
               +(alpha/(2.d0*sqrt(pi)))*exp(-(1.d0/alpha)**2*x0**2))
          a = a
       enddo
    enddo

    precip = precip/pc   ! convert back to mm/a 

  end subroutine glint_precip

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Private subroutines and functions follow
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function satvap(temp,kelvin)

    ! Calculates the saturated vapour pressure using the 
    ! Clausius-Clapyron equation (from Roe 2002)
    ! Note that by default the units of temperature are degC,
    ! unless the kelvin flag is present and set, in which
    ! case they are Kelvin.
    !*RV Saturated vapour pressure in Pascals.

    use glimmer_physcon, only: trpt
    implicit none

    ! Arguments

    real(dp),        intent(in) :: temp    ! Temperature ($^{\circ}$C or K)
    logical,optional,intent(in) :: kelvin  ! Set if temperature is in Kelvin

    ! Internal variables

    real(dp) :: ts
    real(dp),parameter :: e0 = 6.112d0             ! This is in millibars, but multiplied by 100 below.
    real(dp),parameter :: c1 = 17.67, c2=243.5d0   ! These names from Roe(2002)

    ! Beginning of code

    ! First check to see if kelvin is present and adjust accordingly

    if (present(kelvin)) then
       if (kelvin) then 
          ts = temp - trpt   !trpt = 273.15
       else 
          ts = temp
       endif
    else
       ts = temp
    endif

    ! finally, the function itself

    satvap = 100.d0 * e0 * exp(c1*ts/(c2+ts))

  end function satvap

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  function calc_w0(u,v,h0,dx,dy)

    ! Calculates the mean vertical velocity field, 
    ! based on the horizontal flow, and the topography.
    ! 
    ! Note that:
    ! \begin{itemize}
    ! \item All input arrays must be of the same rank and size.
    ! \item The vertical velocity is calculated from the divergence
    !       of the horizontal wind over the topography.
    !       \[w_0=u\frac{\partial h_0}{\partial x}+
    !       v\frac{\partial h_0}{\partial y} \]
    ! \item Differentiation is done with conventional centred-differences,
    !       except at the corners, where uncentred differences are employed.
    ! \end{itemize}

    !*RV An array of the same size as \texttt{u} is 
    !*RV returned. The units are m/s.

    implicit none

    ! Arguments

    real(dp),dimension(:,:),intent(in) :: u   ! The $x$ component of the mean wind field (m/s)
    real(dp),dimension(:,:),intent(in) :: v   ! The $y$ component of the mean wind field (m/s)
    real(dp),dimension(:,:),intent(in) :: h0  ! The topography (m)
    real(dp),               intent(in) :: dx  ! The $x$ grid-spacing (m)
    real(dp),               intent(in) :: dy  ! The $y$ grid-spacing (m)

    ! Returned array

    real(dp),dimension(size(u,1),size(u,2)) :: calc_w0

    ! Internal variables

    integer :: i,j,nx,ny

    ! Beginning of code

    nx=size(u,1) ; ny=size(u,2)

    ! main block loop

    do i=2,nx-1
       do j=2,ny-1
          calc_w0(i,j) = u(i,j)*(h0(i+1,j)-h0(i-1,j))/(2*dx) +v(i,j)*(h0(i,j+1)-h0(i,j-1))/(2.d0*dy)
       enddo
    enddo

    ! top and bottom rows

    do i=2,nx-1
       calc_w0(i,1)  = u(i,1) *(h0(i+1,1) -h0(i-1,1)) /(2*dx) +v(i,1) *(h0(i,2) -h0(i,1))   /dy
       calc_w0(i,ny) = u(i,ny)*(h0(i+1,ny)-h0(i-1,ny))/(2*dx) +v(i,ny)*(h0(i,ny)-h0(i,ny-1))/dy
    enddo

    ! left and right columns

    do j=2,ny-1
       calc_w0(1,j)  = u(1,j) *(h0(2,j) -h0(1,j))   /dx +v(1,j) *(h0(1,j+1) -h0(1,j-1)) /(2*dy)
       calc_w0(nx,j) = u(nx,j)*(h0(nx,j)-h0(nx-1,j))/dx +v(nx,j)*(h0(nx,j+1)-h0(nx,j-1))/(2*dy)
    enddo

    ! corners

    calc_w0(1,1)  = u(1,1)  *(h0(2,1)  -h0(1,1))    /dx +v(1,1)  *(h0(1,2)  -h0(1,1))    /dy
    calc_w0(1,ny) = u(1,ny) *(h0(2,ny) -h0(1,ny))   /dx +v(1,ny) *(h0(1,ny) -h0(1,ny-1)) /dy
    calc_w0(nx,1) = u(nx,1) *(h0(nx,1) -h0(nx-1,1)) /dx +v(nx,1) *(h0(nx,2) -h0(nx,1))   /dy
    calc_w0(nx,ny)= u(nx,ny)*(h0(nx,ny)-h0(nx-1,ny))/dx +v(nx,ny)*(h0(nx,ny)-h0(nx,ny-1))/dy

  end function calc_w0

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  real(dp) function error_func(y)

    ! The error function
    !
    ! The error function may be approximated by:
    ! \[ \mathrm{erf}(y)=1-(\gamma_1 t+\gamma_2 t^2+\gamma_3 t^3)\exp(-y^2)\]
    ! with
    ! \[t=\frac{1}{1+\gamma_0 y}\]  
    ! and
    ! \[\gamma_0 = 0.47047, \]
    ! \[\gamma_1 = 0.3480242, \]
    ! \[\gamma_2 = -0.0958798, \]
    ! \[\gamma_3 = 0.7478556 \]
    ! (from Abramowitz and Stegun 1965)
    ! However, this doesn't seem to be right for $\mathtt{y}<0$, but ok for $\mathtt{y}\geq 0$.
    ! Since the input is always $>0$, this isn't a problem here.
    !*RV The value of the error function at \texttt{y}.

    implicit none

    real(dp),intent(in) :: y ! The independent variable

    real(dp) :: t
    real(dp),parameter :: g0 =  0.47047d0
    real(dp),parameter :: g1 =  0.3480242d0
    real(dp),parameter :: g2 = -0.0958798d0
    real(dp),parameter :: g3 =  0.7478556d0

    t = 1.d0 / (1.d0 + g0*y)

    error_func = 1.d0 - (g1*t + g2*t**2 + g3*t**3)*exp(-y**2)

  end function error_func

end module glint_precip_param
