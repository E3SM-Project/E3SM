!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   test_integrate.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module test_integ

  use glimmer_global, only : dp, sp
  implicit none

contains
  
  ! integrate $\int_0^\pi\a*sin(x+b)dx$

  real(dp) function dfp(x,p)
    implicit none
    real(dp), intent(in) :: x
    real(dp), intent(in), dimension(:) :: p
    
    dfp = p(1)*sin(p(2)*x)
  end function dfp
  
  real(dp) function df(x)
    implicit none
    real(dp), intent(in) :: x

    real(dp), dimension(2) :: p = (/3.d0,0.5d0/)
    
    df = p(1)*sin(p(2)*x)
  end function df

  real(dp) function sfp(x,p)
    implicit none
    real(dp), intent(in) :: x
    real(dp), intent(in), dimension(:) :: p
    
    sfp =  p(1)*sin(p(2)*x)
  end function sfp
  
  real(dp) function sf(x)
    implicit none
    real(dp), intent(in) :: x

    real(dp), dimension(2) :: p = (/3.,0.5/)
    
    sf = p(1)*sin(p(2)*x)
  end function sf

end module test_integ

!TODO - Move this program to another directory so that this directory contains no programs?

program test_integrate

  !*FD test numerical integration schemes
  !WHL - Leaving the sp variables for now since they are used to test
  !      single-precision Romberg integration subroutines

  use test_integ
  use glimmer_physcon, only : pi
  use glimmer_global, only : dp, sp
  use glimmer_integrate
  implicit none

  real(dp), dimension(2) :: dprms = (/3.d0,0.5d0/)
  real(sp), dimension(2) :: sprms = (/3.0,0.50/)

  real(dp) :: dupper = 4.d0*pi
  real(sp) :: supper = 4.0*pi

  write(*,*) 'dp p ',romberg_int(dfp,0.d0,dupper,dprms)
  write(*,*) 'sp p ',romberg_int(sfp,0.,supper,sprms)
  write(*,*) 'dp   ',romberg_int(df,0.d0,dupper)
  write(*,*) 'sp   ',romberg_int(sf,0.,supper)

end program test_integrate



