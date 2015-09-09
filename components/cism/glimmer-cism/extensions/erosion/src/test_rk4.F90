!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   test_rk4.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
! testing rk4 code

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

program testrk4
  use rk4module
  implicit none

  real(kind=dp), dimension(2) :: startx,x
  real(kind=dp) :: h,t0,t1,deltat,t
  integer i, nok, nbad

  external test1

  startx(1) = 2
  startx(2) = 1
  t0 = 0
  t1 = 1.1
  deltat = 0.1

  x = startx
  h = 0.01
  do i=0,int((t1-t0)/deltat)-1
     t = t0+i*deltat
     call odeint(x,t,t+deltat, 0.000001d0, h, 0.d0, nok,nbad, test1)
     write(*,*) t, x(1),x(2), nok,nbad
  end do

end program testrk4

subroutine test1(t,x,dxdt)
  use glimmer_global, only : dp
  implicit none
  real(kind=dp), intent(in) :: t
  real(kind=dp), intent(in), dimension(2) :: x
  real(kind=dp), intent(out), dimension(2) :: dxdt

  dxdt(1) = x(2)
  dxdt(2) = x(1)+x(2)
end subroutine test1
