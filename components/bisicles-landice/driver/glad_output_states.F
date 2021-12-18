!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glad_output_states.F90 - part of the Community Ice Sheet Model (CISM)  
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

module glad_output_states

  ! This module defines routines for computing the output state variables that CISM sends
  ! to a climate model.

  use glimmer_global, only : dp
  use glimmer_paramets, only : thk0
  use glide_types, only : glide_global_type, glide_geometry
  
  implicit none
  private

  public :: set_output_states  ! set state fields output to a climate model

contains

  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine set_output_states(instance, &
       ice_covered, topo, ice_sheet_grid_mask)

    use glad_type, only : glad_instance
    
    ! Arguments ----------------------------------------------------------------------------

    type(glad_instance), intent(in) :: instance
    real(dp),dimension(:,:),intent(out) :: ice_covered  ! whether each grid cell is ice-covered [0,1]
    real(dp),dimension(:,:),intent(out) :: topo         ! output surface elevation (m)
    real(dp),dimension(:,:),intent(out) :: ice_sheet_grid_mask !mask of ice sheet grid coverage

    ! Internal variables ----------------------------------------------------------------------

    integer :: nxl, nyl  ! local grid dimensions
    integer :: i, j      ! indices

    ! Begin subroutine code -------------------------------------------------------------------
    
    ! Initialize arrays. This shouldn't be necessary (because, if the below code is
    ! correct, all points should be explicitly assigned some value), but adds some safety
    ! in case any bugs creep into the below code.
    ice_covered(:,:) = 0.d0
    topo(:,:) = 0.d0
    ice_sheet_grid_mask(:,:) = 0.d0
    
    nxl = instance%lgrid%size%pt(1)
    nyl = instance%lgrid%size%pt(2)

    do j = 1, nyl
       do i = 1, nxl
          if (is_in_active_grid(instance%model%geometry, i, j)) then
             ice_sheet_grid_mask(i,j) = 1.d0

             if (is_ice_covered(instance%model%geometry, i, j)) then
                ice_covered(i,j) = 1.d0
             else
                ice_covered(i,j) = 0.d0
             end if

             ! Note that we use the same method for computing topo whether this point is
             ! ice-covered or ice-free. This is in contrast to the method for computing
             ! ice-free topo in glint_upscaling_gcm.
             topo(i,j) = thk0 * instance%model%geometry%usrf(i,j)
             
          else
             ! Note that this logic implies that if (in theory) we had an ice-covered
             ! point outside the "active grid", it will get classified as ice-free for
             ! these purposes. This mimics the logic currently in glint_upscaling_gcm.
             ice_sheet_grid_mask(i,j) = 0.d0
             ice_covered(i,j) = 0.d0
             topo(i,j) = 0.d0
          end if

       end do
    end do

  end subroutine set_output_states

  
  !===================================================================

  logical function is_in_active_grid(geometry, i, j)
    ! Return true if the given point is inside the "active grid". The active grid includes
    ! any point that can receive a positive surface mass balance, which includes any
    ! point classified as land or ice sheet.
    type(glide_geometry), intent(in) :: geometry
    integer, intent(in) :: i, j  ! point of interest

    real(dp) :: usrf     ! surface elevation (m)

    ! TODO(wjs, 2015-03-18) Could the logic here be replaced by the use of some existing
    ! mask? For now I am simply re-implementing the logic that was in glint.

    usrf = thk0 * geometry%usrf(i,j)

    if (usrf > 0.d0) then
       ! points not at sea level are assumed to be land or ice sheet
       is_in_active_grid = .true.
    else
       is_in_active_grid = .false.
    end if

  end function is_in_active_grid
  
  !===================================================================

  logical function is_ice_covered(geometry, i, j)
    ! Return true if the given point is ice-covered

    use glad_constants, only : min_thck

    type(glide_geometry), intent(in) :: geometry
    integer, intent(in) :: i, j  ! point of interest

    real(dp) :: thck     ! ice thickness (m)

    ! TODO(wjs, 2015-03-18) The logic here should probably be replaced by the use of some
    ! existing mask. For now I am simply re-implementing the logic that was in glint.

    thck = thk0 * geometry%thck(i,j)

    if (thck > min_thck) then
       is_ice_covered = .true.
    else
       is_ice_covered = .false.
    end if

  end function is_ice_covered

end module glad_output_states
