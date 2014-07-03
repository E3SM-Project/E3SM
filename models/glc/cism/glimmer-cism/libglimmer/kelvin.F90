!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   kelvin.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!> module for calculating zeroth order Kelvin functions and their derivatives.
!! Both single and double precision versions are provided
!!
!! \author Magnus Hagdorn
!! \date June 2000

module kelvin

  implicit none

  integer, private, parameter :: sp=selected_real_kind(6)   !< integer storing kind value for variables with a precision of at least 6 digits (single precision)
  integer, private, parameter :: dp=selected_real_kind(15)  !< integer storing kind value for variables with a precision of at least 15 digits (double precision)
  real(kind=dp), private, parameter :: gamma=0.577215664901532860606512d0 !< Euler's constant
  real(kind=dp), private, parameter :: pi=3.14159265358979323846264338328d0 !< pi
  integer, private :: j_max = 40 !< maximum number of iterations
  real(kind=dp), private :: tolerance=1.d-10 !< the tolerance

  interface ber
     module procedure d_ber, s_ber
  end interface
  interface bei
     module procedure d_bei, s_bei
  end interface
  interface ker
     module procedure d_ker, s_ker
  end interface
  interface kei
     module procedure d_kei, s_kei
  end interface

  interface dber
     module procedure d_dber, s_dber
  end interface
  interface dbei
     module procedure d_dbei, s_dbei
  end interface
  interface dker
     module procedure d_dker, s_dker
  end interface
  interface dkei
     module procedure d_dkei, s_dkei
  end interface

contains
  !> set tolerance and maximum number of iterations
  subroutine set_kelvin(tol, jmax)
    implicit none
    real(kind=dp), intent(in) :: tol
    integer, intent(in) :: jmax
    j_max = jmax
    tolerance = tol
  end subroutine set_kelvin

  function d_ber(x)
    implicit none
    real(kind=dp) :: d_ber
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_d_ber
    real(kind=dp) :: factorial
    real(kind=dp) :: sign
    integer :: j

    p_d_ber = 0.d0
    factorial = 1.d0

    d_ber = 1.d0
    arg = (x/2.d0)**4
    arg_d = arg
    sign = -1.d0

    j=1
    do while (j < j_max)
       p_d_ber = d_ber
       factorial = factorial*2*j*(2*j-1.d0)
       d_ber = d_ber + sign*arg_d/(factorial*factorial)
       if (abs(d_ber-p_d_ber) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function d_ber

  function d_bei(x)
    implicit none
    real(kind=dp) :: d_bei
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_d_bei
    real(kind=dp) :: factorial
    real(kind=dp) :: sign
    integer :: j

    p_d_bei = 1.d12
    factorial = 1.d0

    arg = (x/2.d0)**2
    d_bei = arg
    arg_d = arg*arg*arg
    arg = arg*arg
    sign = -1.d0

    j=1
    do while (j < j_max)
       p_d_bei = d_bei
       factorial = factorial*2*j*(2*j+1.d0)
       d_bei = d_bei + sign*arg_d/(factorial*factorial)
       if (abs(d_bei-p_d_bei) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function d_bei

  function d_ker(x)
    implicit none
    real(kind=dp) :: d_ker
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_d_ker
    real(kind=dp) :: factorial
    real(kind=dp) :: phi
    real(kind=dp) :: sign
    integer :: j

    p_d_ker = 0.d0
    factorial = 1.d0

    arg = (x/2.d0)**4
    arg_d = arg
    sign = -1.d0
    phi = 0.d0
    d_ker = -(log(x/2.d0)+gamma)*d_ber(x)+(pi/4.d0)*d_bei(x)

    j=1
    do while (j < j_max)
       p_d_ker = d_ker
       factorial = factorial*2*j*(2*j-1.d0)
       phi = phi + 1.d0/(2.d0*j-1.d0) + 1.d0/(2.d0*j)
       d_ker = d_ker + sign*phi*arg_d/(factorial*factorial)
       if (abs(d_ker-p_d_ker) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function d_ker

  function d_kei(x)
    implicit none
    real(kind=dp) :: d_kei
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_d_kei
    real(kind=dp) :: factorial
    real(kind=dp) :: phi
    real(kind=dp) :: sign
    integer :: j

    p_d_kei = 0.d0
    factorial = 1.d0

    arg = (x/2.d0)**2
    sign = -1.d0
    phi = 1.d0
    d_kei = -(log(x/2.d0)+gamma)*d_bei(x)-(pi/4.d0)*d_ber(x)+arg
    arg_d = arg
    arg = arg*arg
    arg_d = arg_d*arg

    j=1
    do while (j < j_max)
       p_d_kei = d_kei
       factorial = factorial*2*j*(2*j+1.d0)
       phi = phi + 1.d0/(2.d0*j+1.d0) + 1.d0/(2.d0*j)
       d_kei = d_kei + sign*phi*arg_d/(factorial*factorial)
       if (abs(d_kei-p_d_kei) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function d_kei

  function s_ber(x)
    implicit none
    real(kind=sp) :: s_ber
    real(kind=sp), intent(in) :: x

    s_ber = real(d_ber(real(x,kind=dp)),kind=sp)
  end function s_ber

  function s_bei(x)
    implicit none
    real(kind=sp) :: s_bei
    real(kind=sp), intent(in) :: x

    s_bei = real(d_bei(real(x,kind=dp)),kind=sp)
  end function s_bei

  function s_ker(x)
    implicit none
    real(kind=sp) :: s_ker
    real(kind=sp), intent(in) :: x

    s_ker = real(d_ker(real(x,kind=dp)),kind=sp)
  end function s_ker

  function s_kei(x)
    implicit none
    real(kind=sp) :: s_kei
    real(kind=sp), intent(in) :: x

    s_kei = real(d_kei(real(x,kind=dp)),kind=sp)
  end function s_kei

  function d_dber(x)
    implicit none
    real(kind=dp) :: d_dber
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_d_dber
    real(kind=dp) :: factorial
    real(kind=dp) :: sign
    integer :: j

    p_d_dber = 0.d0
    factorial = 1.d0

    d_dber = 0.d0
    arg = (x/2.d0)**4
    arg_d = (x/2.d0)**3
    sign = -1.d0

    j=1
    do while (j < j_max)
       p_d_dber = d_dber
       factorial = factorial*2*j*(2*j-1.d0)
       d_dber = d_dber + sign*2.d0*j*arg_d/(factorial*factorial)
       if (abs(d_dber-p_d_dber) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function d_dber

  function d_dbei(x)
    implicit none
    real(kind=dp) :: d_dbei
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_d_dbei
    real(kind=dp) :: factorial
    real(kind=dp) :: sign
    integer :: j

    p_d_dbei = 1.d12
    factorial = 1.d0

    arg = (x/2.d0)**4
    arg_d = arg*(x/2.d0)
    d_dbei = (x/2.d0)
    sign = -1.d0

    j=1
    do while (j < j_max)
       p_d_dbei = d_dbei
       factorial = factorial*2*j*(2*j+1.d0)
       d_dbei = d_dbei + sign*(2.d0*j+1.d0)*arg_d/(factorial*factorial)
       if (abs(d_dbei-p_d_dbei) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function d_dbei

  function d_dker(x)
    implicit none
    real(kind=dp) :: d_dker
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_d_dker
    real(kind=dp) :: factorial
    real(kind=dp) :: phi
    real(kind=dp) :: sign
    integer :: j

    p_d_dker = 0.d0
    factorial = 1.d0

    arg = (x/2.d0)**4
    arg_d = (x/2.d0)**3
    sign = -1.d0
    phi = 0.d0
    d_dker = -(log(x/2.d0)+gamma)*d_dber(x)-d_ber(x)/x+(pi/4.d0)*d_dbei(x)

    j=1
    do while (j < j_max)
       p_d_dker = d_dker
       factorial = factorial*2*j*(2*j-1.d0)
       phi = phi + 1.d0/(2.d0*j-1.d0) + 1.d0/(2.d0*j)
       d_dker = d_dker + sign*phi*2.d0*j*arg_d/(factorial*factorial)
       if (abs(d_dker-p_d_dker) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function d_dker

  function d_dkei(x)
    implicit none
    real(kind=dp) :: d_dkei
    real(kind=dp), intent(in) :: x
    
    real(kind=dp) :: arg, arg_d
    real(kind=dp) :: p_d_dkei
    real(kind=dp) :: factorial
    real(kind=dp) :: phi
    real(kind=dp) :: sign
    integer :: j

    p_d_dkei = 0.d0
    factorial = 1.d0

    arg = (x/2.d0)
    sign = -1.d0
    phi = 1.d0
    d_dkei = -(log(x/2.d0)+gamma)*d_dbei(x)-d_bei(x)/x-(pi/4.d0)*d_dber(x)+arg
    arg_d = arg**5
    arg = arg**4

    j=1
    do while (j < j_max)
       p_d_dkei = d_dkei
       factorial = factorial*2*j*(2*j+1.d0)
       phi = phi + 1.d0/(2.d0*j+1.d0) + 1.d0/(2.d0*j)
       d_dkei = d_dkei + sign*phi*(2.d0*j+1.d0)*arg_d/(factorial*factorial)
       if (abs(d_dkei-p_d_dkei) < tolerance) exit
       arg_d = arg_d*arg
       sign = -sign 
       j = j+1
    end do
  end function d_dkei

  function s_dber(x)
    implicit none
    real(kind=sp) :: s_dber
    real(kind=sp), intent(in) :: x

    s_dber = real(d_dber(real(x,kind=dp)),kind=sp)
  end function s_dber

  function s_dbei(x)
    implicit none
    real(kind=sp) :: s_dbei
    real(kind=sp), intent(in) :: x

    s_dbei = real(d_dbei(real(x,kind=dp)),kind=sp)
  end function s_dbei

  function s_dker(x)
    implicit none
    real(kind=sp) :: s_dker
    real(kind=sp), intent(in) :: x

    s_dker = real(d_dker(real(x,kind=dp)),kind=sp)
  end function s_dker

  function s_dkei(x)
    implicit none
    real(kind=sp) :: s_dkei
    real(kind=sp), intent(in) :: x

    s_dkei = real(d_dkei(real(x,kind=dp)),kind=sp)
  end function s_dkei

end module kelvin



