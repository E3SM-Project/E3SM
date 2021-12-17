!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glissade_grid_operators.F90 - part of the Community Ice Sheet Model (CISM)  
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
!
! This module contains various grid operators for the Glissade dycore, including routines 
! for computing gradients and interpolating between staggered and unstaggered grids.
!
! Author: William Lipscomb
!         Los Alamos National Laboratory
!         Group T-3, MS B216
!         Los Alamos, NM 87545
!         USA
!         <lipscomb@lanl.gov>
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module glissade_grid_operators

    use glimmer_global, only: dp
    use glimmer_log
    use glide_types  ! HO_GRADIENT_MARGIN_*
    use parallel

    implicit none

    private
    public :: glissade_stagger, glissade_unstagger,    &
              glissade_centered_gradient, glissade_upstream_gradient,    &
              glissade_gradient_at_edges, glissade_vertical_average

    logical, parameter :: verbose_gradient = .false.

contains

!----------------------------------------------------------------------------

  subroutine glissade_stagger(nx,           ny,        &
                              var,          stagvar,   &
                              ice_mask,     stagger_margin_in)

    !----------------------------------------------------------------
    ! Given a variable on the unstaggered grid (dimension nx, ny), interpolate
    ! to find values on the staggered grid (dimension nx-1, ny-1).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx,ny), intent(in) ::    &
       var                      ! unstaggered field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    integer, dimension(nx,ny), intent(in), optional ::        &
       ice_mask                 ! = 1 where values are included in the average, else = 0
                                ! Typically ice_mask = 1 where ice is present (or thck > thklim), else = 0
                                ! Note: ice_mask is not needed if stagger_margin = 0

    integer, intent(in), optional ::   &
       stagger_margin_in        ! 0 = use all values when interpolating
                                !   may be appropriate when computing stagusrf and stagthck on land
                                ! 1 = use only values where ice_mask = 1
                                !   preferable for tracers (e.g., temperature, flwa) and ocean margins

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp) :: sumvar, summask
    integer :: stagger_margin

    if (present(stagger_margin_in)) then
       stagger_margin = stagger_margin_in
    else
       stagger_margin = 0  ! default is to average over all cells, including those where ice is absent
    endif

    if (stagger_margin == 1 .and. .not.present(ice_mask)) then
       call write_log('Must pass in ice_mask to compute staggered field with stagger_margin = 1', GM_FATAL)
    endif

    stagvar(:,:) = 0.d0

    if (stagger_margin == 0) then

       ! Average over all four neighboring cells

       do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          stagvar(i,j) = (var(i,j+1) + var(i+1,j+1) + var(i,j) + var(i+1,j)) / 4.d0
       enddo
       enddo  

    elseif (stagger_margin == 1) then

       ! Average over cells where ice_mask = 1

       do j = 1, ny-1     ! all vertices
       do i = 1, nx-1
          sumvar = ice_mask(i,j+1)*var(i,j+1) + ice_mask(i+1,j+1)*var(i+1,j+1)  &
                 + ice_mask(i,j)  *var(i,j)   + ice_mask(i+1,j)  *var(i+1,j)
          summask = real(ice_mask(i,j+1) + ice_mask(i+1,j+1) + ice_mask(i,j) + ice_mask(i+1,j), dp)
          if (summask > 0.d0) stagvar(i,j) = sumvar / summask
       enddo
       enddo  

    endif

  end subroutine glissade_stagger

!----------------------------------------------------------------------------

  subroutine glissade_unstagger(nx,           ny,          &
                                stagvar,      unstagvar,   &
                                vmask,        stagger_margin_in)

    !----------------------------------------------------------------
    ! Given a variable on the staggered grid (dimension nx-1, ny-1), interpolate
    ! to find values on the staggered grid (dimension nx, ny).
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny                   ! horizontal grid dimensions

    real(dp), dimension(nx-1,ny-1), intent(in) ::    &
       stagvar                  ! staggered field, defined at cell vertices

    real(dp), dimension(nx,ny), intent(out) ::    &
       unstagvar                ! unstaggered field, defined at cell centers

    integer, dimension(nx-1,ny-1), intent(in), optional  ::        &
       vmask                    ! = 1 for vertices where the value is used in the average, else = 0
                                ! Note: The user needs to compute this mask in the calling subroutine.
                                !       It will likely be based on the scalar ice mask, but the details are left open.

    integer, intent(in), optional ::   &
       stagger_margin_in        ! 0 = use all values when interpolating
                                ! 1 = use only values where vmask = 1

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j
    real(dp) :: sumvar, summask
    integer :: stagger_margin

    if (present(stagger_margin_in)) then
       stagger_margin = stagger_margin_in
    else
       stagger_margin = 0  ! default is to average over all cells, including those where ice is absent
    endif

    if (stagger_margin == 1 .and. .not.present(vmask)) then
       call write_log('Must pass in vmask to compute unstaggered field with stagger_margin = 1', GM_FATAL)
    endif

    unstagvar(:,:) = 0.d0

    if (stagger_margin == 0) then

       ! Average over all four neighboring cells

       do j = 2, ny-1   ! loop does not include outer row of cells
       do i = 2, nx-1
          unstagvar(i,j) = (stagvar(i,j) + stagvar(i-1,j) + stagvar(i,j-1) + stagvar(i-1,j-1)) / 4.d0
       enddo
       enddo  

    elseif (stagger_margin == 1) then

       ! Average over vertices with vmask = 1

       do j = 2, ny-1   ! loop does not include outer row of cells
       do i = 2, nx-1
          sumvar = vmask(i-1,j)  *stagvar(i-1,j)   + vmask(i,j)  *stagvar(i,j)  &
                 + vmask(i-1,j-1)*stagvar(i-1,j-1) + vmask(i,j-1)*stagvar(i,j-1)  
          summask = real(vmask(i-1,j) + vmask(i,j) + vmask(i-1,j-1) + vmask(i,j-1), dp)
          if (summask > 0.d0) unstagvar(i,j) = sumvar / summask
       enddo
       enddo  

    endif

    ! Fill in halo values
    call parallel_halo(unstagvar)

  end subroutine glissade_unstagger

!****************************************************************************

  subroutine glissade_centered_gradient(nx,           ny,        &
                                        dx,           dy,        &
                                        field,                   &
                                        df_dx,        df_dy,     &
                                        ice_mask,                &
                                        gradient_margin_in,      &
                                        usrf,                    &
                                        land_mask,               &
                                        max_slope)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient is evaluated at the four neighboring points and is second-order accurate.
    !
    ! There are several choices for computing gradients at the ice margin:
    !
    ! HO_MARGIN_GRADIENT_ALL = 0: All neighbor values are used to compute the gradient, including 
    !  values in ice-free cells.  This convention is used by Glide, but performs poorly for 
    !  ice shelves with a sudden drop in ice thickness and surface elevation at the margin.
    !
    ! HO_MARGIN_GRADIENT_ICE_LAND = 1: Values in ice-covered and/or land cells are used to compute 
    !  the gradient, but values in ice-free ocean cells are ignored.  Where required values are 
    !  missing, the gradient is set to zero. For land-based problems this reduces to option (0) 
    !  (except where ice-free land rises above the ice surface), and for ocean-based problems 
    !  this reduces to option (2).
    !
    !  NOTE: Ice-free land cells contribute to the gradient only where their elevation lies below
    !        the elevation of the adjacent ice-covered cell. This constraint prevents nunataks from
    !        causing large, unrealistic velocities.
    !
    ! HO_MARGIN_GRADIENT_ICE_ONLY = 2: Only values in ice-covered cells (i.e., cells with thck > thklim) 
    !  are used to compute gradients.  Where required values are missing, the gradient is set to zero.
    !  This option works well at shelf margins but less well for land margins (e.g., the Halfar test case).
    !
    ! Since option (1) generally works well at both land and shelf boundaries, it is the default.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       field                    ! input scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components of input field, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells)
                                !    if one or more values is masked out, construct df_fx and df_dy from the others
                                ! 2: use values in ice-covered cells only
                                !    if one or more values is masked out, construct df_fx and df_dy from the others

    real(dp), dimension(nx,ny), intent(in), optional ::       &
       usrf                     ! ice surface elevation (required for gradient_margin = HO_GRADIENT_ICE_LAND)

    integer, dimension(nx,ny), intent(in), optional ::        &
       land_mask                ! = 1 for land cells, else = 0 (required for gradient_margin = HO_GRADIENT_ICE_LAND)

    real(dp), intent(in), optional :: &
       max_slope                ! maximum slope allowed for surface gradient computations (unitless)

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: gradient_margin

    integer :: i, j

    logical, dimension(nx-1,ny) ::    &
       edge_mask_x               ! edge mask for computing df/dx

    logical, dimension(nx,ny-1) ::    &
       edge_mask_y               ! edge mask for computing df/dy

    real(dp) :: df_dx_north, df_dx_south  ! df_dx at neighboring edges
    real(dp) :: df_dy_east, df_dy_west    ! df_dx at neighboring edges

    !WHL - debug
    real(dp) :: dfdx, dfdy
    integer :: edge_count

    !--------------------------------------------------------
    !   Gradient at vertex(i,j) is based on f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)
    !--------------------------------------------------------

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_ICE_LAND
    endif

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       edge_mask_x(:,:) = .true.       ! true for all edges
       edge_mask_y(:,:) = .true.

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_LAND) then

       if (present(land_mask) .and. present(usrf)) then

          call glissade_edgemask_gradient_margin_ice_land(nx,          ny,         &
                                                          ice_mask,    land_mask,  &
                                                          usrf,                    &
                                                          edge_mask_x, edge_mask_y)
       else
          call write_log('Must pass in land mask and usrf to compute centered gradient with gradient_margin = 1', GM_FATAL)
       endif   ! present(land_mask)

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       ! mask for east and west cell edges
       do j = 1, ny
          do i = 1, nx-1
             if (ice_mask(i,j)==1  .and. ice_mask(i+1,j)==1) then
                edge_mask_x(i,j) = .true.
             else
                edge_mask_x(i,j) = .false.
             endif
          enddo
       enddo
       
       ! mask for north and south edges
       do j = 1, ny-1
          do i = 1, nx
             if (ice_mask(i,j)==1  .and. ice_mask(i,j+1)==1) then
                edge_mask_y(i,j) = .true.
             else
                edge_mask_y(i,j) = .false.
             endif
          enddo
       enddo
       
    endif  ! gradient_margin

    !WHL - debug - Count number of edges that require slope-limiting
    if (present(max_slope)) then
       edge_count = 0
       do j = nhalo+1, ny-nhalo
          do i = nhalo+1, nx-nhalo
             dfdx = (field(i+1,j) - field(i,j)) / dx
             if (abs(dfdx) > max_slope .and. edge_mask_x(i,j)) then
                edge_count = edge_count + 1
             endif
             dfdy = (field(i,j+1) - field(i,j)) / dy
             if (abs(dfdy) > max_slope .and. edge_mask_y(i,j)) then
                edge_count = edge_count + 1
             endif
          enddo
       enddo
       edge_count = parallel_reduce_sum(edge_count)
       if (main_task) then
          print*, 'Number of edges:', (nx-2*nhalo)*(ny-2*nhalo)*2
          print*, 'Limit slope: edge_count =', edge_count
       endif
    endif

    ! compute gradient at vertices by averaging gradient at adjacent edges
    ! ignore edges with edge_mask = 0

    do j = 1, ny-1
       do i = 1, nx-1

          ! df/dx
          df_dx_north = (field(i+1,j+1) - field(i,j+1))/ dx
          df_dx_south = (field(i+1,j)   - field(i,j))  / dx

          if (present(max_slope)) then

             if (df_dx_north > 0.0d0) then
                df_dx_north = min(df_dx_north, max_slope)
             else
                df_dx_north = max(df_dx_north, -max_slope)
             endif

             if (df_dx_south > 0.0d0) then
                df_dx_south = min(df_dx_south, max_slope)
             else
                df_dx_south = max(df_dx_south, -max_slope)
             endif

          endif

          if (edge_mask_x(i,j) .and. edge_mask_x(i,j+1)) then
             df_dx(i,j) = (df_dx_north + df_dx_south) / 2.d0
          elseif (edge_mask_x(i,j)) then
             df_dx(i,j) = df_dx_south
          elseif (edge_mask_x(i,j+1)) then
             df_dx(i,j) = df_dx_north
          else
             df_dx(i,j) = 0.d0
          endif

          ! df/dy
          df_dy_east = (field(i+1,j+1) - field(i+1,j))/ dy
          df_dy_west = (field(i,j+1)   - field(i,j))  / dy

          if (present(max_slope)) then

             if (df_dy_east > 0.0d0) then
                df_dy_east = min(df_dy_east, max_slope)
             else
                df_dy_east = max(df_dy_east, -max_slope)
             endif

             if (df_dy_west > 0.0d0) then
                df_dy_west = min(df_dy_west, max_slope)
             else
                df_dy_west = max(df_dy_west, -max_slope)
             endif

          endif

          if (edge_mask_y(i,j) .and. edge_mask_y(i+1,j)) then
             df_dy(i,j) = (df_dy_east + df_dy_west) / 2.d0
          elseif (edge_mask_y(i,j)) then
             df_dy(i,j) = df_dy_west
          elseif (edge_mask_y(i+1,j)) then
             df_dy(i,j) = df_dy_east
          else
             df_dy(i,j) = 0.d0
          endif

       enddo  ! i
    enddo     ! j
   
    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Centered gradient:'
       print*, ' '
       print*, 'df_dx:'
       do j = ny-1, 1, -1
!!          do i = 1, nx-1
          do i = 1, nx/2
             write(6,'(f9.6)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo
       
       print*, ' '
       print*, 'df_dy:'
       do j = ny-1, 1, -1
!!          do i = 1, nx-1
          do i = 1, nx/2
             write(6,'(f9.6)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_centered_gradient

!****************************************************************************

  subroutine glissade_upstream_gradient(nx,           ny,        &
                                        dx,           dy,        &
                                        field,                   &
                                        df_dx,        df_dy,     &
                                        ice_mask,     usrf,      &
                                        gradient_margin_in,      &
                                        accuracy_flag_in,        &
                                        land_mask,               &
                                        max_slope)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    !  compute its gradient (df_dx, df_dy) on the staggered grid (dimension nx-1, ny-1).
    ! The gradient can be evaluated at one upstream edge (for first-order accuracy) 
    !  or at two upstream edges (for second-order accuracy).
    ! The reason to take a one-sided gradient is to damp checkerboard noise
    !  that often arises with a centered gradient.
    !
    ! Note: Upstream is defined by the direction of higher surface elevation.
    !  For df_dx, the edge gradients are upstream in the y direction,
    !  and for df_dy, the edge gradients are upstream in the x direction.
    !
    ! See comments in subroutine glissade_centered_gradient about the 
    !  various values of gradient_margin. 
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       field                    ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny-1), intent(out) ::    &
       df_dx, df_dy             ! gradient components, defined at cell vertices

    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    real(dp), dimension(nx,ny), intent(in) ::       &
       usrf                     ! ice surface elevation (required to determine upstream direction)

    integer, intent(in), optional ::    &
       accuracy_flag_in         ! = 1 for 1st order, 2 for 2nd order

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells)
                                !    if one or more values is masked out, construct df_fx and df_dy from the others
                                ! 2: use values in ice-covered cells only
                                !    if one or more values is masked out, construct df_fx and df_dy from the others

    integer, dimension(nx,ny), intent(in), optional ::        &
       land_mask                ! = 1 for land cells, else = 0 (required for gradient_margin = HO_GRADIENT_ICE_LAND)

    real(dp), intent(in), optional :: &
       max_slope               ! maximum slope allowed for surface gradient computations (unitless)

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: gradient_margin, accuracy_flag
    integer :: i, j
    integer :: summask
    real(dp) :: sum1, sum2

    logical, dimension(nx-1,ny) ::    &
       edge_mask_x               ! edge mask for computing df/dx

    logical, dimension(nx,ny-1) ::    &
       edge_mask_y               ! edge mask for computing df/dy

    real(dp) :: df_dx_north, df_dx_north2
    real(dp) :: df_dx_south, df_dx_south2
    real(dp) :: df_dy_east, df_dy_east2
    real(dp) :: df_dy_west, df_dy_west2

    !--------------------------------------------------------
    !   First-order upstream gradient at vertex(i,j) is based on two points out of f(i:i+1,j:j+1)
    ! 
    !   (i,j+1)  |  (i+1,j+1)
    !   -------(i,j)----------
    !   (i,j)    |  (i+1,j)
    !
    !   Second-order gradient is based on four points in the upstream direction
    !--------------------------------------------------------

    if (present(accuracy_flag_in)) then
       accuracy_flag = accuracy_flag_in
    else
       accuracy_flag = 2   ! default to second-order
    endif

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_ICE_LAND
    endif

    ! Set integer edge mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       edge_mask_x(:,:) = .true.       ! true for all edges
       edge_mask_y(:,:) = .true.

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_LAND) then

       if (present(land_mask)) then

          call glissade_edgemask_gradient_margin_ice_land(nx,          ny,         &
                                                          ice_mask,    land_mask,  &
                                                          usrf,                    &
                                                          edge_mask_x, edge_mask_y)
       else
          call write_log('Must pass in land mask to compute upstream gradient with gradient_margin = 1', GM_FATAL)
       endif   ! present(land_mask)

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       ! mask for east and west cell edges
       do j = 1, ny
          do i = 1, nx-1
             if (ice_mask(i,j)==1  .and. ice_mask(i+1,j)==1) then
                edge_mask_x(i,j) = .true.
             else
                edge_mask_x(i,j) = .false.
             endif
          enddo
       enddo
       
       ! mask for north and south cell edges
       do j = 1, ny-1
          do i = 1, nx
             if (ice_mask(i,j)==1  .and. ice_mask(i,j+1)==1) then
                edge_mask_y(i,j) = .true.
             else
                edge_mask_y(i,j) = .false.
             endif
          enddo
       enddo
       
    endif  ! gradient_margin

    if (accuracy_flag == 1) then   ! first-order accurate

       do j = 1, ny-1
          do i = 1, nx-1

             if (edge_mask_x(i,j) .or. edge_mask_x(i,j+1)) then
                
                ! Compute df_dx by taking upstream gradient

                df_dx_north = (field(i+1,j+1) - field(i,j+1)) / dx
                df_dx_south = (field(i+1,j) - field(i,j)) / dx

                if (present(max_slope)) then

                   if (df_dx_north > 0.0d0) then
                      df_dx_north = min(df_dx_north, max_slope)
                   else
                      df_dx_north = max(df_dx_north, -max_slope)
                   endif
                   
                   if (df_dx_south > 0.0d0) then
                      df_dx_south = min(df_dx_south, max_slope)
                   else
                      df_dx_south = max(df_dx_south, -max_slope)
                   endif

                endif

                sum1 = usrf(i+1,j+1) + usrf(i,j+1)
                sum2 = usrf(i+1,j) + usrf(i,j)

                if (sum1 > sum2) then   ! north is upstream; use north edge gradient if possible

                   if (edge_mask_x(i,j+1)) then
                      df_dx(i,j) = df_dx_north
                   else
                      df_dx(i,j) = df_dx_south
                   endif

                else  ! south is upstream; use south edge gradient if possible

                   if (edge_mask_x(i,j)) then
                      df_dx(i,j) = df_dx_south
                   else
                      df_dx(i,j) = df_dx_north
                   endif
                   
                endif   ! sum1 > sum2

             else    ! both adjacent edge masks = F; punt

                df_dx(i,j) = 0.d0

             endif   ! adjacent edge_mask = T

             if (edge_mask_y(i,j) .or. edge_mask_y(i+1,j)) then
 
                ! Compute df_dy by taking upstream gradient
             
                df_dy_east = (field(i+1,j+1) - field(i+1,j)) / dy
                df_dy_west = (field(i,j+1) - field(i,j)) / dy

                if (present(max_slope)) then

                   if (df_dy_east > 0.0d0) then
                      df_dy_east = min(df_dy_east, max_slope)
                   else
                      df_dy_east = max(df_dy_east, -max_slope)
                   endif

                   if (df_dy_west > 0.0d0) then
                      df_dy_west = min(df_dy_west, max_slope)
                   else
                      df_dy_west = max(df_dy_west, -max_slope)
                   endif
                   
                endif

                sum1 = usrf(i+1,j+1) + usrf(i+1,j)
                sum2 = usrf(i,j+1) + usrf(i,j)
             
                if (sum1 > sum2) then   ! east is upstream; use east edge gradient if possible

                   if (edge_mask_y(i+1,j)) then
                      df_dy(i,j) = df_dy_east
                   else
                      df_dy(i,j) = df_dy_west
                   endif

                else  ! west is upstream; use west edge gradient if possible

                   if (edge_mask_y(i,j)) then
                      df_dy(i,j) = df_dy_west
                   else
                      df_dy(i,j) = df_dy_east
                   endif

                endif   ! sum1 > sum2

             else    ! both adjacent edge masks = F; punt
                
                df_dy(i,j) = 0.d0

             endif   ! adjacent edge mask = T

          enddo
       enddo

    else    ! second-order accurate

       do j = 2, ny-2   ! loop does not include all of halo
          do i = 2, nx-2

             if (edge_mask_x(i,j) .or. edge_mask_x(i,j+1)) then
                
                ! Compute df_dx by taking upstream gradient

                df_dx_north2 = (field(i+1,j+2) - field(i,j+2)) / dx
                df_dx_north  = (field(i+1,j+1) - field(i,j+1)) / dx
                df_dx_south  = (field(i+1,j)   - field(i,j))   / dx
                df_dx_south2 = (field(i+1,j-1) - field(i,j-1)) / dx

                if (present(max_slope)) then

                   if (df_dx_north > 0.0d0) then
                      df_dx_north = min(df_dx_north, max_slope)
                   else
                      df_dx_north = max(df_dx_north, -max_slope)
                   endif
                   
                   if (df_dx_north2 > 0.0d0) then
                      df_dx_north2 = min(df_dx_north2, max_slope)
                   else
                      df_dx_north2 = max(df_dx_north2, -max_slope)
                   endif

                   if (df_dx_south > 0.0d0) then
                      df_dx_south = min(df_dx_south, max_slope)
                   else
                      df_dx_south = max(df_dx_south, -max_slope)
                   endif

                   if (df_dx_south2 > 0.0d0) then
                      df_dx_south2 = min(df_dx_south2, max_slope)
                   else
                      df_dx_south2 = max(df_dx_south2, -max_slope)
                   endif

                endif

                sum1 = usrf(i+1,j+1) + usrf(i,j+1) + usrf(i+1,j+2) + usrf(i,j+2)
                sum2 = usrf(i+1,j) + usrf(i,j) + usrf(i+1,j-1) + usrf(i,j-1)

                if (sum1 > sum2) then   ! north is upstream; use north edge gradients if possible

                   if (edge_mask_x(i,j+1) .and. edge_mask_x(i,j+2)) then
                      df_dx(i,j) = 1.5d0 * df_dx_north - 0.5d0 * df_dx_north2
                   elseif (edge_mask_x(i,j+1)) then   ! revert to first order
                      df_dx(i,j) = df_dx_north
                   else      ! first-order downstream
                      df_dx(i,j) = df_dx_south
                   endif
                                            
                else    ! south is upstream; use south edge gradients if possible

                   if (edge_mask_x(i,j) .and. edge_mask_x(i,j-1)) then
                      df_dx(i,j) = 1.5d0 * df_dx_south - 0.5d0 * df_dx_south2
                   elseif (edge_mask_x(i,j)) then   ! revert to first order
                      df_dx(i,j) = df_dx_south
                   else      ! first-order downstream
                      df_dx(i,j) = df_dx_north
                   endif
                   
                endif   ! sum1 > sum2

             else   ! both adjacent edge masks = F; punt

                df_dx(i,j) = 0.d0

             endif  ! adjacent edge mask = T

             if (edge_mask_y(i,j) .or. edge_mask_y(i+1,j)) then

                ! Compute df_dy by taking upstream gradient

                df_dy_east2 = (field(i+2,j+1) - field(i+2,j)) / dy
                df_dy_east  = (field(i+1,j+1) - field(i+1,j)) / dy
                df_dy_west  = (field(i,j+1)   - field(i,j))   / dy
                df_dy_west2 = (field(i-1,j+1) - field(i-1,j)) / dy

                if (present(max_slope)) then
                   
                   if (df_dy_east > 0.0d0) then
                      df_dy_east = min(df_dy_east, max_slope)
                   else
                      df_dy_east = max(df_dy_east, -max_slope)
                   endif

                   if (df_dy_east2 > 0.0d0) then
                      df_dy_east2 = min(df_dy_east2, max_slope)
                   else
                      df_dy_east2 = max(df_dy_east2, -max_slope)
                   endif
                   
                   if (df_dy_west > 0.0d0) then
                      df_dy_west = min(df_dy_west, max_slope)
                   else
                      df_dy_west = max(df_dy_west, -max_slope)
                   endif

                   if (df_dy_west2 > 0.0d0) then
                      df_dy_west2 = min(df_dy_west2, max_slope)
                   else
                      df_dy_west2 = max(df_dy_west2, -max_slope)
                   endif

                endif
                
                ! determine upstream direction

                sum1 = usrf(i+1,j+1) + usrf(i+1,j) + usrf(i+2,j+1) + usrf(i+2,j)
                sum2 = usrf(i,j+1) + usrf(i,j) + usrf(i-1,j+1) + usrf(i-1,j)
             
                if (sum1 > sum2) then  ! east is upstream; use east edge gradients if possible

                   if (edge_mask_y(i+1,j) .and. edge_mask_y(i+2,j)) then
                      df_dy(i,j) = 1.5d0 * df_dy_east - 0.5d0 * df_dy_east2
                   elseif (edge_mask_y(i+1,j)) then   ! revert to first order
                      df_dy(i,j) = df_dy_east
                   else    ! first-order downstream
                      df_dy(i,j) = df_dy_west
                   endif

                else   ! west is upstream; use west edge gradients if possible

                   if (edge_mask_y(i,j) .and. edge_mask_y(i-1,j)) then
                      df_dy(i,j) = 1.5d0 * df_dy_west - 0.5d0 * df_dy_west2
                   elseif (edge_mask_y(i,j)) then   ! revert to first order
                      df_dy(i,j) = df_dy_west
                   else    ! first_order downstream
                      df_dy(i,j) = df_dy_east
                   endif

                endif   ! sum1 > sum2

             else       ! both adjacent edge masks = F; punt

                df_dy(i,j) = 0.d0

             endif      ! adjacent edge mask = T

          enddo     ! i
       enddo     ! j

       ! fill in halo values
       call staggered_parallel_halo(df_dx)
       call staggered_parallel_halo(df_dy)

    endif   ! first or second order accurate

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'upstream df_dx:'
       do j = ny-2, 2, -1
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo

       print*, ' '
       print*, 'upstream df_dy:'
       do j = ny-2, 2, -1
          do i = 1, nx-1
             write(6,'(f7.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo

    endif

  end subroutine glissade_upstream_gradient

!****************************************************************************

  subroutine glissade_gradient_at_edges(nx,           ny,        &
                                        dx,           dy,        &
                                        field,                   &
                                        df_dx,        df_dy,     &
                                        gradient_margin_in,      &
                                        ice_mask,     land_mask, &
                                        usrf,                    &
                                        max_slope)

    !----------------------------------------------------------------
    ! Given a scalar variable f on the unstaggered grid (dimension nx, ny),
    ! compute its gradient (df_dx, df_dy) at cell edges (i.e., the C grid):
    ! df_dx at the midpoint of the east edge and df_dy at the midpoint of
    ! the north edge.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions

    real(dp), intent(in) ::     &
       dx, dy                   ! grid cell length and width                                           
                                ! assumed to have the same value for each grid cell  

    real(dp), dimension(nx,ny), intent(in) ::       &
       field                    ! scalar field, defined at cell centers

    real(dp), dimension(nx-1,ny), intent(out) ::    &
       df_dx                    ! x gradient component, defined on east cell edge

    real(dp), dimension(nx,ny-1), intent(out) ::    &
       df_dy                    ! y gradient component, defined on north cell edge

    integer, intent(in), optional ::    &
       gradient_margin_in       ! 0: use all values when computing gradient (including zeroes where ice is absent)
                                ! 1: use values in ice-covered and/or land cells (but not ocean cells)
                                !    if one or more values is masked out, set gradient to zero
                                ! 2: use values in ice-covered cells only
                                !    if one or more values is masked out, set gradient to zero

    integer, dimension(nx,ny), intent(in), optional ::        &
       ice_mask,     &          ! = 1 where ice is present, else = 0
       land_mask                ! = 1 for land cells, else = 0 (required for gradient_margin = HO_GRADIENT_ICE_LAND)

    real(dp), dimension(nx,ny), intent(in), optional ::       &
       usrf                     ! ice surface elevation (required for gradient_margin = HO_GRADIENT_ICE_LAND)

    real(dp), intent(in), optional :: &
       max_slope               ! maximum slope allowed for surface gradient computations (unitless)

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: gradient_margin
    integer :: i, j

    logical, dimension(nx-1,ny) ::    &
       edge_mask_x               ! edge mask for computing df/dx

    logical, dimension(nx,ny-1) ::    &
       edge_mask_y               ! edge mask for computing df/dy

    !--------------------------------------------------------
    !   Gradient at east edge(i,j) is based on f(i:i+1,j)
    !   Gradient at north edge(i,j) is based on f(i,j:j+1)
    !
    !   |             |
    !   |   (i,j+1)   |
    !   |             |
    !   |             |
    !   ----df_dy------------------
    !   |             |  
    !   |             |
    !   |   (i,j)   df_dx   (i+1,j)
    !   |             |
    !   |             |
    !   |--------------
    !
    !--------------------------------------------------------

    if (present(gradient_margin_in)) then
       gradient_margin = gradient_margin_in
    else
       gradient_margin = HO_GRADIENT_MARGIN_ICE_LAND
    endif

    ! Set integer edge mask based on gradient_margin.

    if (gradient_margin == HO_GRADIENT_MARGIN_ALL) then

       edge_mask_x(:,:) = .true.       ! true for all edges
       edge_mask_y(:,:) = .true.

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_LAND) then

       if (present(land_mask) .and. present(usrf)) then

          call glissade_edgemask_gradient_margin_ice_land(nx,          ny,         &
                                                          ice_mask,    land_mask,  &
                                                          usrf,                    &
                                                          edge_mask_x, edge_mask_y)
       else
          call write_log('Must pass in land mask and usrf to compute edge gradient with gradient_margin = 1', GM_FATAL)
       endif   ! present(land_mask)

    elseif (gradient_margin == HO_GRADIENT_MARGIN_ICE_ONLY) then

       ! mask for east and west cell edges
       do j = 1, ny
          do i = 1, nx-1
             if (ice_mask(i,j)==1  .and. ice_mask(i+1,j)==1) then
                edge_mask_x(i,j) = .true.
             else
                edge_mask_x(i,j) = .false.
             endif
          enddo
       enddo
       
       ! mask for north and south cell edges
       do j = 1, ny-1
          do i = 1, nx
             if (ice_mask(i,j)==1  .and. ice_mask(i,j+1)==1) then
                edge_mask_y(i,j) = .true.
             else
                edge_mask_y(i,j) = .false.
             endif
          enddo
       enddo
       
    endif  ! gradient_margin

    ! Compute the gradients where edge_ mask = 1

    ! df_dx
    do j = 1, ny
       do i = 1, nx-1
          if (edge_mask_x(i,j)) then
             df_dx(i,j) = (field(i+1,j) - field(i,j)) / dx
          else
             df_dx(i,j) = 0.d0
          endif
       enddo    ! i
    enddo       ! j

    ! df_dy
    do j = 1, ny-1
       do i = 1, nx
          if (edge_mask_y(i,j)) then
             df_dy(i,j) = (field(i,j+1) - field(i,j)) / dy
          else
             df_dy(i,j) = 0.d0
          endif
       enddo    ! i
    enddo       ! j

    if (present(max_slope)) then

       ! limit df_dx
       do j = 1, ny
          do i = 1, nx-1
             if (df_dx(i,j) > 0.0d0) then
                df_dx(i,j) = min(df_dx(i,j), max_slope)
             else
                df_dx(i,j) = max(df_dx(i,j), -max_slope)
             endif
          enddo
       enddo
       
       ! limit df_dy
       do j = 1, ny-1
          do i = 1, nx
             if (df_dy(i,j) > 0.0d0) then
                df_dy(i,j) = min(df_dy(i,j), max_slope)
             else
                df_dy(i,j) = max(df_dy(i,j), -max_slope)
             endif
          enddo
       enddo

    endif  ! present(max_slope)

    if (verbose_gradient .and. main_task) then
       print*, ' '
       print*, 'Edge gradient:'
       print*, ' '
       print*, 'df_dx:'
       do j = ny-1, 1, -1
          do i = 1, nx
             write(6,'(f8.4)',advance='no') df_dx(i,j)
          enddo
          print*, ' '
       enddo
       print*, ' '
       print*, 'df_dy:'
       do j = ny, 1, -1
          do i = 1, nx-1
             write(6,'(f8.4)',advance='no') df_dy(i,j)
          enddo
          print*, ' '
       enddo
    endif

  end subroutine glissade_gradient_at_edges

!****************************************************************************

  subroutine glissade_edgemask_gradient_margin_ice_land(nx,          ny,         &
                                                        ice_mask,    land_mask,  &
                                                        usrf,                    &
                                                        edge_mask_x, edge_mask_y)
    
    !----------------------------------------------------------------
    ! Compute edge masks required for option gradient_margin = HO_MARGIN_GRADIENT_ICE_LAND.
    ! Values in ice-covered and/or land cells are used to compute the gradient, but values 
    !  in ice-free ocean cells are ignored.  Where required values are missing, the gradient 
    !  is set to zero.  
    ! NOTE: Ice-free land cells contribute to the gradient only where their elevation lies below
    !        the elevation of the adjacent ice-covered cell. This constraint prevents nunataks from
    !        causing large, unrealistic velocities.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::      &
       nx, ny                   ! horizontal grid dimensions
   
    integer, dimension(nx,ny), intent(in) ::        &
       ice_mask                 ! = 1 where ice is present, else = 0

    integer, dimension(nx,ny), intent(in)  ::        &
       land_mask                ! = 1 for land cells, else = 0

    real(dp), dimension(nx,ny), intent(in)  ::       &
       usrf                     ! ice surface elevation

    logical, dimension(nx-1,ny), intent(out) ::    &
       edge_mask_x               ! edge mask for computing df/dx

    logical, dimension(nx,ny-1), intent(out) ::    &
       edge_mask_y               ! edge mask for computing df/dy

    !----------------------------------------------------------------
    ! Local variables
    !----------------------------------------------------------------

    integer :: i, j

    !WHL - debug
    ! Set to true to use new edge mask, in which ice-free land cells contribute to the gradient
    !  only if adjacent to an ice-covered cell.
    ! Set to false to use the old method where all ice-free land cells contribute to the gradient.
    !TODO - Remove the old edge mask option?

    logical, parameter :: new_edgemask = .true.

 if (new_edgemask) then

    ! compute mask for east and west cell edges
    do j = 1, ny
       do i = 1, nx-1
          if (( ice_mask(i,j)==1 .and.  ice_mask(i+1,j)==1)  .or.  &
              ( ice_mask(i,j)==1 .and. land_mask(i+1,j)==1 .and. usrf(i,j)   >= usrf(i+1,j)) .or.  &
              (land_mask(i,j)==1 .and.  ice_mask(i+1,j)==1 .and. usrf(i+1,j) >= usrf(i,j))) then
             edge_mask_x(i,j) = .true.
          else
             edge_mask_x(i,j) = .false.
          endif
       enddo
    enddo
    
    ! compute mask for north and south cell edges
    do j = 1, ny-1
       do i = 1, nx
          if (( ice_mask(i,j)==1 .and.  ice_mask(i,j+1)==1)  .or.  &
              ( ice_mask(i,j)==1 .and. land_mask(i,j+1)==1 .and. usrf(i,j)   >= usrf(i,j+1)) .or.  &
              (land_mask(i,j)==1 .and.  ice_mask(i,j+1)==1 .and. usrf(i,j+1) >= usrf(i,j))) then
             edge_mask_y(i,j) = .true.
          else
             edge_mask_y(i,j) = .false.
          endif
       enddo
    enddo

 else   ! old edge mask

    ! compute mask for east and west cell edges
    do j = 1, ny
       do i = 1, nx-1
           if ((ice_mask(i,j)==1  .and. ice_mask(i+1,j)==1)  .or.  &
               (ice_mask(i,j)==1  .and. land_mask(i+1,j)==1) .or.  &
               (land_mask(i,j)==1 .and. ice_mask(i+1,j)==1)  .or.  &
               (land_mask(i,j)==1 .and. land_mask(i+1,j)==1)) then
             edge_mask_x(i,j) = .true.
          else
             edge_mask_x(i,j) = .false.
          endif
       enddo
    enddo
          
    ! compute mask for north and south cell edges
    do j = 1, ny-1
       do i = 1, nx
          if ((ice_mask(i,j)==1  .and. ice_mask(i,j+1)==1)  .or.  &
              (ice_mask(i,j)==1  .and. land_mask(i,j+1)==1) .or.  &
              (land_mask(i,j)==1 .and. ice_mask(i,j+1)==1)  .or.  &
              (land_mask(i,j)==1 .and. land_mask(i,j+1)==1)) then
             edge_mask_y(i,j) = .true.
          else
             edge_mask_y(i,j) = .false.
          endif
       enddo
    enddo

 endif  ! new_edgemask

  end subroutine glissade_edgemask_gradient_margin_ice_land

!****************************************************************************

  subroutine glissade_vertical_average(nx,         ny,        &
                                       nz,         sigma,     &
                                       mask,                  &
                                       var,        var_2d)

    !----------------------------------------------------------------
    ! Compute the vertical average of a given variable.
    ! Note: It is assumed that the variable is defined at layer midpoints,
    !       and hence has vertical dimension (nz-1).
    ! Note: This subroutine will work for variables on the staggered
    !       horizontal grid if stagthck is passed in place of thck.
    !----------------------------------------------------------------

    !----------------------------------------------------------------
    ! Input-output arguments
    !----------------------------------------------------------------

    integer, intent(in) ::   &
       nx, ny,            &     ! horizontal grid dimensions
       nz                       ! number of vertical levels

    real(dp), dimension(nz), intent(in) ::    &
       sigma                    ! sigma vertical coordinate

    logical, dimension(nx, ny), intent(in) ::    &
       mask                     ! compute var_2d where mask = .true.

    real(dp), dimension(nz-1,nx, ny), intent(in) ::    &
       var                      ! 3D field to be averaged vertically

    real(dp), dimension(nx, ny), intent(out) ::    &
       var_2d                   ! 2D vertically averaged field

    !--------------------------------------------------------
    ! Local variables
    !--------------------------------------------------------

    integer :: i, j, k

    do j = 1, ny
       do i = 1, nx

          var_2d(i,j) = 0.d0

          if (mask(i,j)) then
             do k = 1, nz-1
                var_2d(i,j) = var_2d(i,j) + var(k,i,j) * (sigma(k+1) - sigma(k))
             enddo
          endif

       enddo
    enddo

  end subroutine glissade_vertical_average

!****************************************************************************

  end module glissade_grid_operators

!****************************************************************************
