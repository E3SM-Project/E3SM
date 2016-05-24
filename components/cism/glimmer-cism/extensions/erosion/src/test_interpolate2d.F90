!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   test_interpolate2d.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
! testing 2d interpolation stuff

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

program test_interpolate
  use glimmer_coordinates
  use glimmer_sparse
  use glimmer_interpolate2d
  implicit none

  type(coordsystem_type) :: coords,hires_coords
  type(sparse_matrix_type) :: interpolator,f2f_interp
  real(kind=dp) :: delta = 0.1,x,y
  integer i,j
  integer, parameter :: numx = 100
  integer, parameter :: numy = 150
  integer, parameter :: inter_numx = 100
  integer, parameter :: inter_numy = 50

  real(kind=dp), dimension(numx,numy) ::  orig_field
  real(kind=dp), dimension(2*numx,2*numy) ::  highres_field

  real(kind=dp), dimension(inter_numx,inter_numy) :: dispx, dispy, interp_field

  ! setup coordsystem
  coords = coordsystem_new(0.d0,0.d0,delta, delta, numx, numy)
  hires_coords = coordsystem_new(0.d0,0.d0,delta/2., delta/2., 2*numx, 2*numy)

  ! setup data
  do j=1,numy
     y=(j-1)*delta
     do i=1,numx
        x=(i-1)*delta
        orig_field(i,j) = calc_data(x,y)
     end do
  end do

  ! generate random displacement field
  call random_number(dispx)
  dispx = dispx * (numx-1)*delta
  call random_number(dispy)
  dispy = dispy * (numy-1)*delta
  
  open(unit=1,file="interp_disp.data",status="unknown")
  call glimmer_init_bilinear(coords,dispx, dispy,interpolator)
  call glimmer_interpolate(interpolator,orig_field,interp_field)  
  do j=1,inter_numy
     do i=1,inter_numx
        write(1,*) dispx(i,j), dispy(i,j), calc_data(dispx(i,j), dispy(i,j)), interp_field(i,j)
     end do
  end do
  close(1)

  open(unit=1,file="interp_coord.data",status="unknown")
  call glimmer_init_bilinear(coords,hires_coords,f2f_interp)
  call glimmer_interpolate(f2f_interp,orig_field,highres_field)
  do j=1,2*numy-1
     y=0.5*(j-1)*delta
     do i=1,2*numx-1
        x=0.5*(i-1)*delta
        write(1,*) x,y,calc_data(x,y),highres_field(i,j)
     end do
  end do
  close(1)

contains 
  function calc_data(x,y)
    implicit none
    real(kind=dp) :: calc_data
    real(kind=dp), intent(in) :: x,y

    calc_data = sin(x)+cos(y)
  end function calc_data
end program test_interpolate

