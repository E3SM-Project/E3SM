!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   erosion_advect.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
! for solving dx/dt=v(x,t)

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module erosion_advect

  use glimmer_global, only : dp
  use glimmer_coordinates

  ! private data
  real(kind=dp), private, dimension(:,:), pointer :: ux, uy
  type(coordsystem_type), private :: coords
  real(kind=dp), private :: eps
  
contains
  subroutine er_advect2d_init(velo_coords)
    !initialise advection
    implicit none
    type(coordsystem_type) :: velo_coords ! coordinate system for velo grid
    
    call coordsystem_allocate(velo_coords, ux)
    call coordsystem_allocate(velo_coords, uy)

    coords = velo_coords

    eps = 0.00001*minval(velo_coords%delta%pt)

  end subroutine er_advect2d_init

  subroutine er_advect2d(times, x0d, y0d, x, y)
    use rk4module
    implicit none
    real(kind=dp), intent(in), dimension(:) :: times   ! array of times at which position should be saved
    real(kind=dp), intent(in) :: x0d, y0d              ! initial conditions
    real(kind=dp), intent(out), dimension(:) :: x, y   ! x component of location

    ! local variables
    integer j
    real(kind=dp) t1,t2
    real(kind=dp), dimension(2) :: xtemp
    integer :: nok, nbad

    t1 = 0.
    xtemp(1) = x0d
    xtemp(2) = y0d
    
    do j=1,size(times)
       t2 = times(j)
       call odeint(xtemp,t1,t2,eps,10.d0,0.d0,nok,nbad,interp_velos)
       x(j) = xtemp(1)
       y(j) = xtemp(2)
       
       t1 = times(j)
    end do

  end subroutine er_advect2d

  subroutine interp_velos(t, x, velo)
    use glimmer_interpolate2d
    implicit none
    real(kind=dp), intent(in) :: t
    real(kind=dp), intent(in), dimension(2) :: x
    real(kind=dp), intent(out), dimension(2) :: velo

    ! local vars
    type(coord_point) :: pnt
    type(coord_ipoint), dimension(4) :: nodes
    real(kind=dp), dimension(4) :: weights

    integer k
    
    pnt%pt(:) = x(:)
    call glimmer_bilinear(coords,pnt,nodes,weights)

    velo(:) = 0.d0
    do k=1,4
       velo(1) = velo(1) + weights(k)*ux(nodes(k)%pt(1),nodes(k)%pt(2))
       velo(2) = velo(2) + weights(k)*uy(nodes(k)%pt(1),nodes(k)%pt(2))
    end do
  end subroutine interp_velos
  
  subroutine set_velos(x,y,factor)
    implicit none
    real(kind=dp), dimension(:,:) :: x,y
    real(kind=dp), optional :: factor
    ! local variables
    real(kind=dp) :: f

    if (present(factor)) then
       f = factor
    else
       f = 1.
    end if

    ux = f*x
    uy = f*y
  end subroutine set_velos
end module erosion_advect
