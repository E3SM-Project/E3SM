!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_interpolate2d.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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
! module for 2D interpolation

#ifdef HAVE_CONFIG_H
#include <config.inc>
#endif

module glimmer_interpolate2d

  use glimmer_sparse

  interface glimmer_init_bilinear
     module procedure glimmer_init_bilinear_disp, glimmer_init_bilinear_coord
  end interface

contains

  subroutine glimmer_interpolate(interpolator, infield, outfield)
    use glimmer_coordinates
    implicit none
    !*FD interpolate data on infield (defined on inccord system) onto outfield using interpolator
    type(sparse_matrix_type), intent(in) :: interpolator       !*FD sparse matrix containing interpolator
    real(kind=dp), dimension(:,:), intent(in)  :: infield !*FD input data
    real(kind=dp), dimension(:,:), intent(out) :: outfield!*FD output field

    ! local variables
    real(kind=dp), dimension(:), allocatable :: linear_out
    integer i,j, k, numx,numy

    numx = size(outfield,1)
    numy = size(outfield,2)
    
    allocate(linear_out(numx*numy))

    call sparse_matrix_vec_prod(interpolator,pack(infield,.true.), linear_out)

    do j=1,numy
       k = (j-1)*numx
       do i=1,numx
          outfield(i,j) = linear_out(k+i)
       end do
    end do

    deallocate(linear_out)
  end subroutine glimmer_interpolate

  subroutine glimmer_init_bilinear_disp(coord, dispx, dispy, bilinear)
    use glimmer_coordinates
    implicit none
    !*FD initialise sparse matrix defining interpolation given by displacement field
    type(coordsystem_type), intent(in)             :: coord        !*FD coordinate system of the input field
    real(kind=dp), dimension(:,:), intent(in) :: dispx, dispy !*FD displacement field 
    type(sparse_matrix_type)                  :: bilinear     !*FD sparse matrix containing interpolation

    ! local variables
    integer :: i,j, k, lini, numx,numy
    type(coord_point) :: point
    type(coord_ipoint) :: this_node
    type(coord_ipoint), dimension(4) :: nodes
    real(kind=dp), dimension(4) :: weights


    numx = size(dispx,1)
    numy = size(dispx,2)

    call new_sparse_matrix(numx*numy*4,bilinear)

    do j=1,numy
       this_node%pt(2)=j
       do i=1,numx
          this_node%pt(1)=i
          lini = i+(j-1)*numx
          point%pt(1)=dispx(i,j)
          point%pt(2)=dispy(i,j)
          call glimmer_bilinear(coord, point, nodes, weights)
          do k=1,4
             call sparse_insert_val(bilinear,lini,coordsystem_linearise2d(coord,nodes(k)),weights(k))
          end do
       end do
    end do

  end subroutine glimmer_init_bilinear_disp

  subroutine glimmer_init_bilinear_coord(incoord, outcoord, bilinear)
    use glimmer_coordinates
    implicit none
    !*FD initialise sparse matrix defining interpolation between two coordinate systems
    type(coordsystem_type), intent(in) :: incoord      !*FD coordinate system of the input field
    type(coordsystem_type), intent(in) :: outcoord     !*FD coordinate system of the output field
    type(sparse_matrix_type)           :: bilinear     !*FD sparse matrix containing interpolation

    ! local variables
    integer :: i,j, k, lini
    type(coord_point) :: point
    type(coord_ipoint) :: this_node
    type(coord_ipoint), dimension(4) :: nodes
    real(kind=dp), dimension(4) :: weights

    call new_sparse_matrix(outcoord%size%pt(1)*outcoord%size%pt(2)*4,bilinear)    

    do j=1,outcoord%size%pt(2)
       this_node%pt(2)=j
       do i=1,outcoord%size%pt(1)
          this_node%pt(1)=i
          lini = coordsystem_linearise2d(outcoord,this_node)
          point = coordsystem_get_coord(outcoord,this_node)
          call glimmer_bilinear(incoord, point, nodes, weights)
          do k=1,4
             call sparse_insert_val(bilinear,lini,coordsystem_linearise2d(incoord,nodes(k)),weights(k))
          end do
       end do
    end do
  end subroutine glimmer_init_bilinear_coord

  subroutine glimmer_bilinear(coord, point, nodes, weights)
    use glimmer_coordinates
    implicit none
    !*FD bilinear interpolation
    type(coordsystem_type), intent(in)            :: coord     !*FD coordinate system to operate on
    type(coord_point), intent(in)             :: point     !*FD desired point
    type(coord_ipoint), dimension(4), intent(out) :: nodes !*FD array containing indicies into field
    real(kind=dp), dimension(4), intent(out) :: weights   !*FD array of weights

    type(coord_point) :: pnt


    ! get lower left node of cell into which point falls
    nodes(1) = coordsystem_get_llnode(coord,point)

    ! check if point is inside coord-system
    if (.not.coordsystem_point_inside(coord,point)) then
       nodes(1)%pt(:) = 1
       nodes(2)%pt(:) = 1
       nodes(3)%pt(:) = 1
       nodes(4)%pt(:) = 1
       weights(:) = 0.d0
       return
    end if

    nodes(2)%pt(1) = nodes(1)%pt(1) + 1
    nodes(2)%pt(2) = nodes(1)%pt(2)

    nodes(3)%pt(1) = nodes(1)%pt(1) + 1
    nodes(3)%pt(2) = nodes(1)%pt(2) + 1

    nodes(4)%pt(1) = nodes(1)%pt(1)
    nodes(4)%pt(2) = nodes(1)%pt(2) + 1

    pnt%pt(:) = (point%pt(:)-coord%origin%pt(:)-(nodes(1)%pt(:)-1)*coord%delta%pt(:))*coord%delta_r%pt(:)

    if (nodes(1)%pt(1).eq.coord%size%pt(1) .and. nodes(1)%pt(2).eq.coord%size%pt(2)) then
       nodes(2:4)%pt(1) = 1
       nodes(2:4)%pt(2) = 1
       weights(1) = 1.d0
       weights(2:4) = 0.d0
       return
    else if (nodes(1)%pt(1).eq.coord%size%pt(1)) then
       nodes(2)%pt(:) = 1
       nodes(3)%pt(:) = 1
       weights(1) = (1-pnt%pt(2))
       weights(2) = 0.d0
       weights(3) = 0.d0
       weights(4) = pnt%pt(2)
       return
    else if (nodes(1)%pt(2).eq.coord%size%pt(2)) then
       nodes(3)%pt(:) = 1
       nodes(4)%pt(:) = 1
       weights(1) = (1-pnt%pt(1))
       weights(2) = pnt%pt(1)
       weights(3) = 0.d0
       weights(4) = 0.d0
       return
    else
       weights(1) = (1-pnt%pt(1))*(1-pnt%pt(2))
       weights(2) = pnt%pt(1)*(1-pnt%pt(2))
       weights(3) = pnt%pt(1)*pnt%pt(2)
       weights(4) = (1-pnt%pt(1))*pnt%pt(2)
    end if
  end subroutine glimmer_bilinear


end module glimmer_interpolate2d
