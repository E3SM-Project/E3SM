
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glimmer_utils.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!TODO - Is there any other utility code to put in this module?
!> Module containing utility code for GLIMMER.
module glimmer_utils

  use glimmer_global, only: dp

  implicit none

  !TODO - Remove these array_bcs functions?  I can't find where they are used
  interface array_bcs
    module procedure array_bcs1d,array_bcs2d
  end interface

  !TODO - Move check_conformal to glint?  Used by glint_interp only.
  interface check_conformal
    module procedure check_conformal_2d_real
  end interface

contains

  !> Returns the value of a 1D array location,checking first for the boundaries.
  !!
  !! the location is wrapped around the array boundaries until it falls within the array
  !! \author The value of the location in question.
  real(dp) function array_bcs1d(array,i)

    ! Arguments

    real(dp),dimension(:),intent(in) :: array !< The array to be indexed.
    integer,intent(in)               :: i     !< The location to be extracted.

    ! Internal variables

    integer :: n,ii

    n=size(array)
    ii=i

    if ((i<=n).and.(i>=1)) then
      array_bcs1d=array(i)
    endif

    do while (ii>n)
      ii=ii-n
    enddo

    do while (ii<1)
      ii=ii+n
    enddo

    array_bcs1d=array(ii)

  end function array_bcs1d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Returns the value of a 1D array location,checking first for the boundaries.
  !!
  !! the location is wrapped around the array boundaries until it falls within the array
  !! as array_bcs1d but for polar boundary conditions
  !! \author The value of the location in question.
  real(dp) function array_bcs_lats(array,i)


    ! Arguments

    real(dp),dimension(:),intent(in) :: array !< The array to be indexed.
    integer,intent(in) :: i !< The location to be extracted.

    ! Internal variables

    integer :: n,ii

    n=size(array)
    ii=i

    if ((i<=n).and.(i>=1)) then
      array_bcs_lats=array(i)
      return
    endif

    if (ii>n) then
      ii=2*n-ii
      array_bcs_lats=-180.d0+array(ii)
    endif

    if (ii<1) then
      ii=1-ii
      array_bcs_lats=180.d0-array(ii)
    endif

  end function array_bcs_lats

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Returns the value of an array 
  !! location, checking first for the boundaries. 
  !! Over-the-pole boundary conditions are implemented here.
  !! \return The value of the location specified.
  real(dp) function array_bcs2d(array,i,j)

    ! Arguments

    real(dp),dimension(:,:),intent(in) :: array !< Array to be indexed
    integer,intent(in) :: i !< The location to be extracted    
    integer,intent(in) :: j !< The location to be extracted

    ! Internal variables

    integer :: nx,ny,ii,jj

    nx=size(array,1) ; ny=size(array,2)

    if ((i>=1).and.(i<=nx).and.(j>=1).and.(j<=ny)) then
      array_bcs2d=array(i,j)
      return
    endif

    ii=i ; jj=j

    if (jj>ny) then
      jj=2*ny-jj
      ii=ii+nx/2
    endif

    if (jj<1) then
      jj=1-jj
      ii=ii+nx/2
    endif

    do while (ii>nx) 
      ii=ii-nx
    enddo

    do while (ii<1)
      ii=ii+nx
    enddo  

    array_bcs2d=array(ii,jj)

  end function array_bcs2d

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Move to glint?  Used by glint_interp only.

  subroutine fix_bcs2d(i,j,nx,ny)

  !> Adjusts array location indices
  !! so that they fall within the domain.

    integer,intent(inout) :: i !< The location of interest
    integer,intent(inout) :: j !< The location of interest
    integer,intent(in) :: nx  !< The size of the domain (number of points in each direction)
    integer,intent(in) :: ny  !< The size of the domain (number of points in each direction)

    if ((i>=1).and.(i<=nx).and.(j>=1).and.(j<=ny)) return

    if (j>ny) then
      j=2*ny-j
      i=i+nx/2
    endif

    if (j<1) then
      j=1-j
      i=i+nx/2
    endif

    do while (i>nx) 
      i=i-nx
    enddo

    do while (i<1)
      i=i+nx
    enddo  

  end subroutine fix_bcs2d

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine check_conformal_2d_real(array1,array2,label)

  !> Checks that two arrays are of the same size.

    use glimmer_log

    real(dp),dimension(:,:),intent(in) :: array1 !< The array 1 to be checked
    real(dp),dimension(:,:),intent(in) :: array2 !< The array 2 to be checked
    character(*),intent(in),optional :: label    !< Optional label, to facilitate bug tracking if the check fails.

    if ((size(array1,1)/=size(array2,1)).or.(size(array1,2)/=size(array2,2))) then
      if (present(label)) then
        call write_log('Non-conformal arrays. Label: '//label,GM_FATAL,__FILE__,__LINE__)
      else
        call write_log('ERROR: Non-conformal arrays. No label',GM_FATAL,__FILE__,__LINE__)
      endif
    endif

  end subroutine check_conformal_2d_real

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> compute horizontal sum for each vertical level
  !!
  !! Calculates the sum of a given three-dimensional field at each
  !! level. The vertical coordinate of the input is the first index of
  !! the array.
  !! \return
  !! A one-dimensional array of the same size as the first dimension of
  !! inp is returned, containing the sum of inp for 
  !! each level.
  function hsum(inp)


    implicit none

    real(dp),dimension(:,:,:),intent(in) :: inp !< The input array. The first index is the vertical, the othe two horizontal.
    real(dp),dimension(size(inp,dim=1))  :: hsum
  
    integer up

    do up=1,size(inp,dim=1)
       hsum(up) = sum(inp(up,:,:))
    end do

  end function hsum

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !TODO - Move lsum to glint?  Used by glint_interp only.
  !> Calculates the sum of a given two-dimensional field along one axis.
  !! Within GLIMMER, this function calculates the mean vertical profile
  !! in a 2D vertical slice. 
  !! \return
  !! A one-dimensional array of the same size as the first dimension of
  !! inp is returned, containing the sum of inp for 
  !! each row.

  function lsum(inp)


    implicit none

    real(dp),dimension(:,:), intent(in) :: inp !< Input array
    real(dp),dimension(size(inp,dim=1)) :: lsum
    
    lsum = sum(inp(:,:),dim=2)

  end function lsum

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  !> Tridiagonal solver. All input/output arrays should have the 
  !! same number of elements.

  subroutine tridiag(a,b,c,x,y)


    real(dp),dimension(:) :: a !< Lower diagonal; a(1) is ignored.
    real(dp),dimension(:) :: b !< Centre diagonal
    real(dp),dimension(:) :: c !< Upper diagonal; c(n) is ignored.
    real(dp),dimension(:) :: x !< Unknown vector
    real(dp),dimension(:) :: y !< Right-hand side

    real(dp),dimension(size(a)) :: aa
    real(dp),dimension(size(a)) :: bb

    integer :: n,i

    n=size(a)

    aa(1) = c(1)/b(1)
    bb(1) = y(1)/b(1)

    do i=2,n
       aa(i) = c(i)/(b(i)-a(i)*aa(i-1))
       bb(i) = (y(i)-a(i)*bb(i-1))/(b(i)-a(i)*aa(i-1))
    end do
    
    x(n) = bb(n)

    do i=n-1,1,-1
       x(i) = bb(i)-aa(i)*x(i+1)
    end do

  end subroutine tridiag

end module glimmer_utils
