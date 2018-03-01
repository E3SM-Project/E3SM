!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glide_lithot1d.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

#include "glide_mask.inc"

! module for 1D temperature calculations in the upper lithosphere

!TODO - This is a local calculation and should be parallel-friendly,
!       but has not yet been tested in parallel.

module glide_lithot1d

  implicit none

contains

  subroutine init_lithot1d(model)

    use glide_types
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    ! allocate memory for 1D code
    allocate(model%lithot%rhs(model%lithot%nlayer))
    allocate(model%lithot%subd(model%lithot%nlayer))
    allocate(model%lithot%diag(model%lithot%nlayer))
    allocate(model%lithot%supd(model%lithot%nlayer))
    
    ! setup coefficient matrix
    model%lithot%subd(:) =    - model%lithot%zfactors(1,:)
    model%lithot%diag(:) = 1. + model%lithot%zfactors(2,:)
    model%lithot%supd(:) =    - model%lithot%zfactors(3,:)
    ! and the boundary conditions
    ! top face
    ! simply match air temperature where no ice and basal temperature where ice
    model%lithot%subd(1) = 0.
    model%lithot%diag(1) = 1.
    model%lithot%supd(1) = 0.
    ! bottom face
    ! keep constant
    model%lithot%subd(model%lithot%nlayer) = 0.
    model%lithot%diag(model%lithot%nlayer) = 1.
    model%lithot%supd(model%lithot%nlayer) = 0.
  end subroutine init_lithot1d

  subroutine calc_lithot1d(model)
    use glide_types
    use glimmer_utils, only: tridiag
    !use glide_mask
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    integer i,j,k

!TODO - I think these loops can be left as is for parallel code.
!       Local calculation, so no issues with computing in halo cells.

    ! loop over grid
    do j=1,model%general%nsn
       do i=1,model%general%ewn
          ! calculate RHS for upper BC
          if (GLIDE_IS_GROUND(model%geometry%thkmask(i,j)) .and. .not. GLIDE_IS_THIN(model%geometry%thkmask(i,j)) ) then
             model%lithot%rhs(1) = model%temper%temp(model%general%upn,i,j) ! ice basal temperature
             model%lithot%mask(i,j) = .true.
          else
             if (model%lithot%mask(i,j)) then
                if (GLIDE_IS_OCEAN(model%geometry%thkmask(i,j))) then
                   model%lithot%rhs(1) = model%lithot%mart
                else if (GLIDE_IS_LAND(model%geometry%thkmask(i,j))) then
                   model%lithot%rhs(1) = model%climate%artm(i,j) ! air temperature outside ice sheet
                end if
             end if
          end if

          if (model%lithot%mask(i,j)) then
             ! calculate RHS for rest
             do k=2,model%lithot%nlayer-1
                model%lithot%rhs(k) = - model%lithot%subd(k)*model%lithot%temp(i,j,k-1) &
                     + (2.-model%lithot%diag(k))*model%lithot%temp(i,j,k) &
                     - model%lithot%supd(k)*model%lithot%temp(i,j,k+1)
             end do
             model%lithot%rhs(model%lithot%nlayer) = model%lithot%temp(i,j,model%lithot%nlayer)

             ! solve tri-diagonal matrix eqn
             call tridiag(model%lithot%subd(1:), &
                  model%lithot%diag(:), &
                  model%lithot%supd(:model%lithot%nlayer), &
                  model%lithot%temp(i,j,:) ,                 &
                  model%lithot%rhs(:))
          end if
       end do
    end do
  end subroutine calc_lithot1d

  subroutine finalise_lithot1d(model)
    use glide_types
    implicit none
    type(glide_global_type),intent(inout) :: model       !*FD model instance

    deallocate(model%lithot%rhs)
    deallocate(model%lithot%subd)
    deallocate(model%lithot%diag)
    deallocate(model%lithot%supd)
  end subroutine finalise_lithot1d

end module glide_lithot1d
