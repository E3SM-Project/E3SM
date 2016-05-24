!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   unittest_glimmer_deriv.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

module unittest_glimmer_deriv
  use glimmer_deriv
  use glimmer_global, only: dp
  use glimmer_physcon, only: pi
  use xls    !TODO - Move xls to libglimmer?
  implicit none
contains

  subroutine setup_sinusoidal_test(f, dfdx, dfdy, dx, dy)
    real(dp), intent(out) :: dx, dy
    real(dp), dimension(10,10,1), intent(out) :: f, dfdx, dfdy

    integer :: i,j
    real(dp) :: x,y

    dx = .1
    dy = .1

    !Set up test functions and their analytical derivatives
    do i = 1,10
       do j = 1,10
          x = dx*(i - 5)
          y = dy*(j - 5)
          f(i,j,1)    = sin(2*pi*x) * cos(3*pi*y)
          dfdx(i,j,1) = 2*pi * cos(2*pi*x) * cos(3*pi*y)
          dfdy(i,j,1) = -3*pi * sin(2*pi*x) * sin(3*pi*y)
       end do
    end do
  end subroutine setup_sinusoidal_test

  subroutine setup_quadratic_test(f, dfdx, dfdy, dx, dy)
    real(dp), intent(out) :: dx, dy
    real(dp), dimension(10,10,1), intent(out) :: f, dfdx, dfdy

    integer :: i,j
    real(dp) :: x,y

    dx = 1
    dy = 1

    do i = 1,10
       do j = 1,10
          x = dx*(i - 5)
          y = dy*(j - 5)
          f(i,j,1)    = x**2 * y**2
          dfdx(i,j,1) = 2 * x * y**2
          dfdy(i,j,1) = 2 * x**2 * y
       end do
    end do
  end subroutine setup_quadratic_test

  subroutine test_glimmer_deriv()
    real(dp), dimension(10,10,1) :: f, dfdx, dfdy
    real(dp) :: dx, dy

    call setup_quadratic_test(f, dfdx, dfdy, dx, dy)

    call write_xls_3d("f.txt",f)
    call write_xls_3d("dfdx.txt",dfdx)
    call write_xls_3d("dfdy.txt",dfdy)
    call check_derivatives(f, dfdx, dfdy, dx, dy)
  end subroutine test_glimmer_deriv

  subroutine check_derivatives(f, dfdx, dfdy, dx, dy)
    real(dp), dimension(:,:,:), intent(in) :: f, dfdx, dfdy
    real(dp), intent(in)::  dx, dy

    integer :: i,j
    real(dp) :: d
    do i = 1, size(f,1)
       do j = 1, size(f,2)
          write(*,*)     i,j,"dx, analytic",dfdx(i,j,1)
          if (i > 1 .and. i < size(f,1)) then
             d = dfdx_3d(f,i,j,1,dx)
             write(*,*) i,j,"dx, centered",d 
          end if

          if (i > 2) then 
             d = dfdx_3d_upwind(f,i,j,1,dx)
             write(*,*) i,j,"dx, upwind  ",d
          end if

          if (i < size(f,1) - 1) then 
             d = dfdx_3d_downwind(f,i,j,1,dx)
             write(*,*) i,j,"dx, downwind",d
          end if

          write(*,*)
          write(*,*)     i,j,"dy, analytic",dfdy(i,j,1)

          if (j > 1 .and. j < size(f,2)) then
             d = dfdy_3d(f,i,j,1,dy)
             write(*,*) i,j,"dy, centered",d 
          end if

          if (j > 2) then 
             d = dfdy_3d_upwind(f,i,j,1,dy)
             write(*,*) i,j,"dy, upwind  ",d
          end if

          if (j < size(f,2) - 1) then 
             d = dfdy_3d_downwind(f,i,j,1,dy)
             write(*,*) i,j,"dy, downwind",d 
          end if

          write(*,*)
       end do
    end do
  end subroutine check_derivatives

end module unittest_glimmer_deriv
