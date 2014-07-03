!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!                                                             
!   glam_grid_operators.F90 - part of the Glimmer Community Ice Sheet Model (Glimmer-CISM)  
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

!TODO - Make sure all these subroutines are currently used, or might be in the future.

! Various grid operators for glide dycore, including routines for computing gradients
! and switching between staggered and unstaggered grids

#ifdef HAVE_CONFIG_H
#include "config.inc"
#endif

#include "glide_mask.inc"

!*FD This module contains functions for computing derivatives numerically, both
!*FD for a single value and for an entire field.
!*FD Note that this module is written with the first index in a matrix corresponding
!*FD to the x (east-west) coordinate.  If this is not the case (i.e. if the first
!*FD index corresponds to the y (north-south) coordinate), then transposition
!*FD will be necessary.  Simply ask for the y-derivative when you mean to ask for
!*FD the x-derivative, and vice versa.

module glam_grid_operators

    use glimmer_global, only: sp, dp

    implicit none

contains

!----------------------------------------------------------------------------

  subroutine glam_geometry_derivs(model)

    ! Compute derivatives of the ice and bed geometry, as well as averaging
    ! them onto the staggered grid

    use glide_types, only: glide_global_type
    use glide_grid_operators, only: stagvarb  ! can we remove this?
    implicit none

    type(glide_global_type), intent(inout) :: model

    call stagthickness(model%geometry% thck, &
                       model%geomderv%stagthck,&
                       model%general%ewn, &
                       model%general%nsn, &
                       model%geometry%usrf, &
                       model%numerics%thklim, &
                       model%geometry%thkmask)

!TODO: Should these calls to stagvarb be replaced by calls to df_field_2d_staggered?    
    call stagvarb(model%geometry%lsrf, &
                  model%geomderv%staglsrf,&
                  model%general%ewn, &
                  model%general%nsn)

    call stagvarb(model%geometry%topg, &
                  model%geomderv%stagtopg,&
                  model%general%ewn, &
                  model%general%nsn)

    model%geomderv%stagusrf = model%geomderv%staglsrf + model%geomderv%stagthck

    call df_field_2d_staggered(model%geometry%usrf, &
                               model%numerics%dew, model%numerics%dns, &
                               model%geomderv%dusrfdew, &
                               model%geomderv%dusrfdns, &
                               model%geometry%thck,     &
                               model%numerics%thklim )

    call df_field_2d_staggered(model%geometry%thck, &
                              model%numerics%dew, model%numerics%dns, &
                              model%geomderv%dthckdew, &
                              model%geomderv%dthckdns, &
                              model%geometry%thck,     &
                              model%numerics%thklim )
 
    !Make sure that the derivatives are 0 where staggered thickness is 0

    where (model%geomderv%stagthck == 0.d0)
           model%geomderv%dusrfdew = 0.d0
           model%geomderv%dusrfdns = 0.d0
           model%geomderv%dthckdew = 0.d0
           model%geomderv%dthckdns = 0.d0
    endwhere

    model%geomderv%dlsrfdew = model%geomderv%dusrfdew - model%geomderv%dthckdew
    model%geomderv%dlsrfdns = model%geomderv%dusrfdns - model%geomderv%dthckdns
      
    !Compute second derivatives.
    
    call d2f_field_stag(model%geometry%usrf, model%numerics%dew, model%numerics%dns, &
                        model%geomderv%d2usrfdew2, model%geomderv%d2usrfdns2, &
                        .false., .false.)

    call d2f_field_stag(model%geometry%thck, model%numerics%dew, model%numerics%dns, &
                        model%geomderv%d2thckdew2, model%geomderv%d2thckdns2, &
                        .false., .false.)
 
  end subroutine glam_geometry_derivs

!----------------------------------------------------------------------------

  subroutine stagthickness(ipvr,opvr,ewn,nsn,usrf,thklim,mask)

  !! A special staggering algorithm that is meant to conserve mass when operating on thickness fields.
  !! This incorporates Ann LeBroque's nunatak fix and the calving front fix.

!NOTE: This subroutine, used by the glam HO dycore, is different from stagvarb, 
!      which is used by the glide SIA dycore.  Here, zero-thickness values are
!      ignored when thickness is averaged over four adjacent grid cells.
!      In stagvarb, zero-thickness values are included in the average.
!      The glam approach works better for calving. 

    implicit none 

    real(dp), intent(out), dimension(:,:) :: opvr 
    real(dp), intent(in), dimension(:,:) :: ipvr
    
    real(dp), intent(in), dimension(:,:) :: usrf
    real(dp), intent(in) :: thklim
    integer, intent(in), dimension(:,:) :: mask
    
    integer :: ewn,nsn,ew,ns,n
    real(dp) :: tot

        do ns = 1,nsn-1
            do ew = 1,ewn-1

                !If any of our staggering points are shelf front, ignore zeros when staggering
                !if (any(GLIDE_IS_CALVING(mask(ew:ew+1, ns:ns+1)))) then  ! in contact with the ocean
                !Use the "only nonzero thickness" staggering criterion for ALL marginal ice. For
                ! reasons that are not entirely clear, this corrects an error whereby the land ice 
                ! margin is defined incorrectly as existing one grid cell too far inland from where 
                ! it should be.  

                if (any(GLIDE_HAS_ICE(mask(ew:ew+1,ns:ns+1)))) then
                    n = 0
                    tot = 0
                    if (abs(ipvr(ew,ns)) > 0.0d0  )then
                        tot = tot + ipvr(ew,ns)
                        n   = n   + 1
                    end if
                    if (abs(ipvr(ew+1,ns)) > 0.0d0  )then
                        tot = tot + ipvr(ew+1,ns)
                        n   = n   + 1
                    end if
                    if (abs(ipvr(ew,ns+1)) > 0.0d0  )then
                        tot = tot + ipvr(ew,ns+1)
                        n   = n   + 1
                    end if
                    if (abs(ipvr(ew+1,ns+1)) > 0.0d0  )then
                        tot = tot + ipvr(ew+1,ns+1)
                        n   = n   + 1
                    end if
                    if (n > 0) then
                        opvr(ew,ns) = tot/n
                    else
                        opvr(ew,ns) = 0
                    end if

                !The following cases relate to Anne LeBroque's fix for nunataks
                !ew,ns cell is ice free:
                else if (ipvr(ew,ns) <= thklim .and. &
                   ((usrf(ew,ns) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >= thklim) &
                    .or. (usrf(ew,ns) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >= thklim))) then
                        opvr(ew,ns) = 0.0

                !ew+1,ns cell is ice free:
                else if (ipvr(ew+1,ns) <= thklim .and. &
                    ((usrf(ew+1,ns) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim) &
                    .or. (usrf(ew+1,ns) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim))) then
                        opvr(ew,ns) = 0.0
    
                !ew,ns+1 cell is ice free:
                else if (ipvr(ew,ns+1) <= thklim .and. &
                    ((usrf(ew,ns+1) >= usrf(ew,ns) .and. ipvr(ew,ns) >= thklim) &
                    .or. (usrf(ew,ns+1) >= usrf(ew+1,ns+1) .and. ipvr(ew+1,ns+1) >= thklim))) then
                        opvr(ew,ns) = 0.0
    
                !ew+1,ns+1 cell is ice free:
                else if (ipvr(ew+1,ns+1) <= thklim .and. &
                    ((usrf(ew+1,ns+1) >= usrf(ew+1,ns) .and. ipvr(ew+1,ns) >=thklim) &
                    .or. (usrf(ew+1,ns+1) >= usrf(ew,ns+1) .and. ipvr(ew,ns+1) >=thklim))) then
                        opvr(ew,ns) = 0.0
               
!                !Standard Staggering   !! Not needed if only-nonzero-thickness staggering scheme is used
!                else
!                        opvr(ew,ns) = (ipvr(ew+1,ns) + ipvr(ew,ns+1) + &
!                                       ipvr(ew+1,ns+1) + ipvr(ew,ns)) / 4.0d0

                end if
  
        end do
    end do

  end subroutine stagthickness

!----------------------------------------------------------------------------

    !------------------------------------------------------------------
    !First Derivative Estimates, Second Order, 2D
    !------------------------------------------------------------------

    !*FD Computes derivative fields of the given function.
    subroutine df_field_2d(f,  &
                           deltax,      deltay,      &
                           out_dfdx,    out_dfdy,    &
                           direction_x, direction_y)

      use parallel
        implicit none
        real(dp), dimension(:, :), intent(in) :: f
        real(dp), intent(in) :: deltax, deltay
        real(dp), dimension(:, :), intent(out) :: out_dfdx, out_dfdy
        real(dp), dimension(:, :), intent(in), optional  :: direction_x, direction_y
        
        logical :: upwind !Whether or not directions for upwinding were provided

        integer :: grad_x, grad_y !Whether to upwind or downwind at the current point

        integer :: nx, ny, x, y

        !Get the size of the field we're working with
        nx = size(f, 1)
        ny = size(f, 2)
        
        upwind = present(direction_x) .and. present(direction_y)

        !For now, we'll use the function calls defined above.
        !Later on we might want to refactor?

!LOOP: all scalar points (uses upwinding and downwinding to avoid stepping out of bounds)
        do x = 1, nx
            do y = 1, ny
                grad_x = 0
                grad_y = 0
                if (upwind) then
                    if (direction_x(x,y) < 0 .and. x > 2) then !Upstream case
                        grad_x = -1
                    else if(direction_x(x,y) > 0 .and. x < nx - 1) then !Downstream case
                        grad_x = 1
                    end if

                    if (direction_y(x,y) < 0 .and. y > 2) then !Upstream case
                        grad_y = -1
                    else if(direction_y(x,y) > 0 .and. y < ny - 1) then !Downstream case
                        grad_y = 1
                    end if
                end if
  
                !For each of the variables in x, y, check whether or not
                !we need to use an upwinding or downwinding differentiation
                !scheme.
                if (x == 1 .or. grad_x > 0) then
                    out_dfdx(x, y) = dfdx_2d_downwind(f, x, y, deltax)
                else if (x == nx .or. grad_x < 0) then
                    out_dfdx(x, y) = dfdx_2d_upwind(f, x, y, deltax)
                else
                    out_dfdx(x, y) = dfdx_2d(f, x, y, deltax)
                end if
                        
                if (y == 1 .or. grad_y > 0) then
                    out_dfdy(x, y) = dfdy_2d_downwind(f, x, y, deltay)
                elseif (y == ny .or. grad_y < 0) then
                    out_dfdy(x, y) = dfdy_2d_upwind(f, x, y, deltay)
                else
                    out_dfdy(x, y) = dfdy_2d(f, x, y, deltay)
                end if
                        
            end do  
        end do

!NOTE:  If halo updates are needed, they should be done at a higher level.

!!        call parallel_halo(out_dfdx)
!!        call parallel_halo(out_dfdy)
        
    end subroutine df_field_2d

!----------------------------------------------------------------------------

    !*FD Computes derivative fields of the given function.  Places the result
    !*FD on a staggered grid.  If periodic in one dimension is set, that 
    !*FD dimension for derivatives must be the same size as the value's dimension.
    !*FD Otherwise, it should be one less

!TODO - This is the subroutine used to compute dusrfdew/ns, dthkdew/ns.
!       Not sure if mods are needed for parallel code.
!TODO - thck and thklim are no longer used in this subroutine.  Remove from argument list.

    subroutine df_field_2d_staggered(f,                  &
                                     deltax,   deltay,   &
                                     out_dfdx, out_dfdy, &
                                     thck,     thklim )

        implicit none
        real(dp), dimension(:, :), intent(in) :: f, thck    ! unstaggered grid
        real(dp), intent(in) :: deltax, deltay, thklim
        real(dp), dimension(:, :), intent(out) :: out_dfdx, out_dfdy  ! staggered grid
        
        integer :: nx, ny, x, y
        
        !Get the size of the field we're working with
        nx = size(f, 1)
        ny = size(f, 2)

        ! intialize to zeros
        out_dfdx = 0.0d0
        out_dfdy = 0.0d0
        
        ! *SFP* old subroutine calls, commented out below but still available, 
        ! use centered diffs on normal thck / surf grids but do nothing special at lateral
        ! boundaries where centered diffs might give unreasonable values (e.g., due to jumping
        ! from a region of non-zero to zero thickness / elevation). New calls access new 
        ! subroutines which attempt to correct for this if/when possible using approx., first-order
        ! accurate one-sided diffs.

!TODO - Use old calls or new calls?

        do x = 1, nx - 1  ! We go to nx - 1 because we're using a staggered grid
            do y = 1, ny - 1
                out_dfdx(x,y) = dfdx_2d_stag(f, x, y, deltax) !*SFP* old call
                out_dfdy(x,y) = dfdy_2d_stag(f, x, y, deltay) !*SFP* old call
!                out_dfdx(x,y) = dfdx_2d_stag_os(f, x, y, deltax, thck, thklim )
!                out_dfdy(x,y) = dfdy_2d_stag_os(f, x, y, deltay, thck, thklim )
            end do
        end do

!TODO - Remove this chunk of code
!               !Deal with periodic boundary conditions.  We will do so by
!               !providing another set of values at the end of each dimension
!               !that contains the derivative of the value off the edge of the
!               !grid.  Because this set of values is at the end, when
!               !x = nx, x+1 = 1.  This identity has been hard-coded below.
!               if (periodic_x) then
!                       do y = 1, ny - 1
!                               out_dfdx(nx,y) = -(f(1, y) + f(1, y+1) - f(nx, y) - f(nx, y+1))/(2*deltax)
!                               out_dfdy(nx,y) = dfdy_2d_stag(f, nx, y, deltay)
!                       end do
!               end if
!               
!               if (periodic_y) then
!                       do x = 1, nx - 1
!                           out_dfdx(x,ny) = dfdx_2d_stag(f, x, ny, deltax)
!                               out_dfdy(x,ny) = -(f(x, 1) + f(x+1, 1) - f(x,ny) - f(x+1, ny))/(2*deltay)
!                       end do
!               end if
!               
!               !Do the corner that hasn't been done if both dimensions are periodic
!               if (periodic_x .and. periodic_y) then
!                       out_dfdx(nx,ny) = (f(1, ny) + f(1, 1) - f(nx, ny) - f(nx, 1))/(2*deltax)
!                       out_dfdy(nx,ny) = (f(nx, 1) + f(1, 1) - f(nx,ny)  - f(1, ny))/(2*deltay)
!               end if
!               
        end subroutine df_field_2d_staggered

!----------------------------------------------------------------------------

!TODO - This 3D subroutine is never called.  Remove it?
 
    subroutine df_field_3d(f,                                  &
                           deltax,      deltay,      deltaz,   &
                           out_dfdx,    out_dfdy,    out_dfdz, &
                           direction_x, direction_y)

    !*FD Computes derivative fields of the given function.
    !*FD The z axis is computed on an irregular grid.

        implicit none
        real(dp), dimension(:, :, :), intent(in) :: f
        real(dp), intent(in) :: deltax, deltay
        real(dp), dimension(:), intent(in) :: deltaz
        real(dp), dimension(:, :, :), intent(out) :: out_dfdx, out_dfdy, out_dfdz

        !Field containing the direction that derivatives should be upwinded in.
        !If 0, centered differences are used.  If negative, then upwinded
        !derivatives (approaching from the negative side) are used.  If
        !positive, then downwinded derivatives (approaching from the positive
        !side) are used.
        real(dp), dimension(:,:), optional :: direction_x, direction_y

 
        integer :: grad_x, grad_y !Sign of the gradient, used for determining upwinding
        integer :: nx, ny, nz, x, y, z
        logical :: upwind

        !Get the size of the field we're working with
        nx = size(f, 2)
        ny = size(f, 3)
        nz = size(f, 1)
        
        upwind = present(direction_x) .and. present(direction_y)

        !For now, we'll use the function calls defined above.
        !Later on we might want to refactor?

!LOOP: all scalar points
!      uses upwinding and downwinding to avoid going out of bounds
        do x = 1, nx
                do y = 1, ny
                        grad_x = 0
                        grad_y = 0
                        if (upwind) then
                            if (direction_x(x,y) < 0 .and. x > 2) then !Upstream case
                                grad_x = -1
                            else if(direction_x(x,y) > 0 .and. x < nx - 1) then !Downstream case
                                grad_x = 1
                            end if

                            if (direction_y(x,y) < 0 .and. y > 2) then !Upstream case
                                grad_y = -1
                            else if(direction_y(x,y) > 0 .and. y < ny - 1) then !Downstream case
                                grad_y = 1
                            end if
                        end if
                        
                        do z = 1, nz
                                !For each of the variables in x, y, check whether or not
                                !we need to use an upwinding or downwinding differentiation
                                !scheme.
                                if (x == 1 .or. grad_x > 0) then
                                        out_dfdx(z, x, y) = dfdx_3d_downwind(f, x, y, z, deltax)
                                        !out_dfdx(x, y, z) = (f(x+1,y,z) - f(x,y,z))/deltax
                                else if (x == nx .or. grad_x < 0) then
                                        out_dfdx(z, x, y) = dfdx_3d_upwind(f, x, y, z, deltax)
                                        !out_dfdx(x, y, z) = (f(x,y,z) - f(x-1,y,z))/deltax
                                else
                                        out_dfdx(z, x, y) = dfdx_3d(f, x, y, z, deltax)
                                end if
                                if (y == 1 .or. grad_y > 0) then
                                        out_dfdy(z, x, y) = dfdy_3d_downwind(f, x, y, z, deltay)
                                        !out_dfdy(x, y, z) = (f(x,y+1,z) - f(x,y,z))/deltay
                                else if (y == ny .or. grad_y < 0) then
                                        out_dfdy(z, x, y) = dfdy_3d_upwind(f, x, y, z, deltay)
                                        !out_dfdy(x, y, z) = (f(x,y,z) - f(x,y-1,z))/deltay
                                else
                                        out_dfdy(z, x, y) = dfdy_3d(f, x, y, z, deltay)
                                end if
                                if (z == 1) then
                                        out_dfdz(z, x, y) = dfdz_3d_downwind_irregular(f, x, y, z, deltaz)
                                else if (z == nz) then
                                        out_dfdz(z, x, y) = dfdz_3d_upwind_irregular(f, x, y, z, deltaz)
                                else
                                        out_dfdz(z, x, y) = dfdz_3d_irregular(f, x, y, z, deltaz)
                                end if
                        end do
                end do  
        end do
        
    end subroutine df_field_3d

!----------------------------------------------------------------------------

!TODO - This 3D subroutine is never called.  Remove it?

    subroutine df_field_3d_stag(f,                                  &
                                deltax,      deltay,      deltaz,   &
                                out_dfdx,    out_dfdy,    out_dfdz)

        !*FD Computes the derivative fields of the given function.  The X and Y
        !*FD derivatives are computed on a staggered grid.  The Z derivative
        !*FD is computed on a nonstaggered but irregular grid.  This means that,
        !*FD if an array of dimensions (n1, n2, n3), the output arrays should
        !*FD be of size (n1 - 1, n2 - 1, n3)

        implicit none
        real(dp), dimension(:, :, :), intent(in) :: f
        real(dp), intent(in) :: deltax, deltay
        real(dp), dimension(:), intent(in) :: deltaz
        real(dp), dimension(:, :, :), intent(out) :: out_dfdx, out_dfdy, out_dfdz
        
        real(dp), dimension(4) :: zDerivs !Temporarily holds derivatives in Z to average
        integer :: nx, ny, nz, x, y, z

        !Get the size of the field we're working with
        nx = size(f, 1)
        ny = size(f, 2)
        nz = size(f, 3)
        
!LOOP: all scalar points
!      uses upwinding and downwinding to avoid going out of bounds

        do x = 1, nx - 1
                do y = 1, ny - 1
                        do z = 1, nz
                                !We will never have to compute upstream and downstream
                                !derivatives in the horizontal (avoided by the staggered scheme),
                                !but we will in the vertical.
                                        out_dfdx(x,y,z) = dfdx_3d_stag(f, x, y, z, deltax)
                                        out_dfdy(x,y,z) = dfdy_3d_stag(f, x, y, z, deltay)
                                        
                                        !Even though we are not staggering in the vertical, the points
                                        !we compute the derivatives at are still staggered in the
                                        !horizontal.  We'll solve this by computing four
                                        !derivatives horizontally around the point requested
                                        !and averaging the results
                                if (z == 1) then
                                    zDerivs(1) = dfdz_3d_downwind_irregular(f, x, y, z, deltaz)
                                    zDerivs(2) = dfdz_3d_downwind_irregular(f, x+1, y, z, deltaz)
                                    zDerivs(3) = dfdz_3d_downwind_irregular(f, x, y+1, z, deltaz)
                                    zDerivs(4) = dfdz_3d_downwind_irregular(f, x+1, y+1, z, deltaz)
                                else if (z == nz) then
                                    zDerivs(1) = dfdz_3d_upwind_irregular(f, x, y, z, deltaz)
                                    zDerivs(2) = dfdz_3d_upwind_irregular(f, x+1, y, z, deltaz)
                                    zDerivs(3) = dfdz_3d_upwind_irregular(f, x, y+1, z, deltaz)
                                    zDerivs(4) = dfdz_3d_upwind_irregular(f, x+1, y+1, z, deltaz)
                                else
                                    zDerivs(1) = dfdz_3d_irregular(f, x, y, z, deltaz)
                                    zDerivs(2) = dfdz_3d_irregular(f, x+1, y, z, deltaz)
                                    zDerivs(3) = dfdz_3d_irregular(f, x, y+1, z, deltaz)
                                    zDerivs(4) = dfdz_3d_irregular(f, x+1, y+1, z, deltaz)
                                end if
                                out_dfdz(x, y, z) = (zDerivs(1) + zDerivs(2) + zDerivs(3) + zDerivs(4)) / 4
                        end do
                end do
        end do          
        
    end subroutine df_field_3d_stag

!----------------------------------------------------------------------------

!TODO - Check the rest of this module for unused functions we might want to remove

    !*FD Computes derivative with respect to x at a given point.
    !*FD Applies periodic boundary conditions if needed.

    function dfdx_2d(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_2d
        
        dfdx_2d = (-.5/delta)*f(i-1, j) + (.5/delta)*f(i+1, j)
        !write(*,*), i, j, f(i,j), ip1, im1, delta, dfdx_2d
    end function dfdx_2d

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to y at a given point

    function dfdy_2d(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_2d
        
        integer :: jp1, jm1
        jp1 = j + 1
        jm1 = j - 1
        if (jp1 == size(f, 2)+1) jp1 = 2
        if (jm1 == 0) jm1 = size(f, 2)-1
        
        dfdy_2d = (-.5/delta)*f(i, j-1) + (.5/delta)*f(i, j+1)
    end function dfdy_2d

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to x at the equivalent
    !*FD point on a staggered grid.

    function dfdx_2d_stag(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_2d_stag
        dfdx_2d_stag = (f(i+1, j) + f(i+1, j+1) - f(i, j) - f(i, j+1))/(2*delta) 
    end function dfdx_2d_stag

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to y at the equivalent
    !*FD point on a staggered grid.

    function dfdy_2d_stag(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_2d_stag
        dfdy_2d_stag = (f(i, j+1) + f(i+1, j+1) - f(i,j) - f(i+1, j))/(2*delta)
    end function dfdy_2d_stag

!----------------------------------------------------------------------------

    function dfdx_2d_stag_os(f_in, i, j, delta, thck, thklim )

    !*SFP* altered/expanded version of above function that uses approx. one-sided
    ! diffs at physical domain edges so as not to overesimate grads there.
        implicit none
        real(dp), dimension(:,:), intent(in) :: f_in, thck
        integer, intent(in) :: i, j
        real(dp), intent(in) :: delta, thklim
        real(dp) :: dfdx_2d_stag_os

        real(dp), dimension(2,2) :: f_array, thck_array
        real(dp), dimension(1:size(f_in,1),1:size(f_in,2)) :: f
        real(dp) :: f_min

        ! initialize vars/arrays to zeros
        dfdx_2d_stag_os = 0.0d0; f_array = 0.0d0; f_min = 0.0d0

        f = f_in

        where( thck <= thklim )
            f = 0
        end where

        f_array(1,1) = f(i,j); f_array(2,1) = f(i+1,j); f_array(1,2) = f(i,j+1); f_array(2,2) = f(i+1,j+1);

        if( sum( f_array/ f_array, MASK = f_array /= 0.0d0 ) == 4.0 )then

            ! normal differencing for interior points
            dfdx_2d_stag_os = (f(i+1,j) + f(i+1,j+1) - f(i,j) - f(i,j+1))/(2*delta)

        elseif( sum( f_array/ f_array, MASK = f_array /= 0.0d0 ) == 3.0 )then

            ! corner; use 2x next closest value
            if( f(i,j) == f_min )then     ! southwest corner point missing: apply value from s.e. point
                dfdx_2d_stag_os = ( f(i+1,j+1) + f(i+1,j) - 2.0*f(i,j+1) )/(2*delta)
            elseif( f(i+1,j) == f_min )then ! southeast corner point missing: apply value from s.w. point
                dfdx_2d_stag_os = ( 2.0*f(i+1,j+1) - f(i,j) - f(i,j+1) )/(2*delta)
            elseif( f(i,j+1) == f_min )then ! northwest corner point missing: apply value from n.e. point
                dfdx_2d_stag_os = ( f(i+1,j+1) + f(i+1,j) - 2.0*f(i,j))/(2*delta)
            elseif( f(i+1,j+1) == f_min )then ! northeast corner point missing: apply value from n.w. point
                dfdx_2d_stag_os = ( 2.0*f(i+1,j) - f(i,j) - f(i,j+1) )/(2*delta)
            endif

        elseif( sum( f_array/ f_array, MASK = f_array /= 0.0d0 ) == 2.0 )then

            ! side; back up and take gradient from points one set of cells in OR use only the single set of
            ! cells available along the differencing direction 
            if( f(i,j) == f_min .and. f(i,j+1) == f_min )then   ! west cells empty
                dfdx_2d_stag_os = (f(i+2,j) + f(i+2,j+1) - f(i+1,j+1) - f(i+1,j))/(2*delta)
            elseif( f(i+1,j) == f_min .and. f(i+1,j+1) == f_min )then   ! east cells empty
                dfdx_2d_stag_os = (f(i,j) + f(i,j+1) - f(i-1,j) - f(i-1,j+1))/(2*delta)
            elseif( f(i,j+1) == f_min .and. f(i+1,j+1) == f_min )then   ! north cells empty
                dfdx_2d_stag_os = (f(i+1,j) - f(i,j) )/(delta)
            elseif( f(i,j) == f_min .and. f(i+1,j) == f_min )then   ! south cells empty
                dfdx_2d_stag_os = (f(i+1,j+1) - f(i,j+1) )/(delta)
            endif

        elseif( sum( f_array/ f_array, MASK = f_array /= 0.0d0 ) == 1.0 )then

            ! isolated; treat by assuming it is part of a 3 block for which the rest of the values are not contained in
            ! the local 2x2 block with indices i:i+1, j:j+1 
            if( f(i,j) /= f_min .and.  f(i+1,j) == f_min .and. f(i+1,j+1) == f_min .and. f(i,j+1) == f_min)then
            ! a northeast corner
                dfdx_2d_stag_os = ( f(i,j) - f(i-1,j) ) / (delta)
            elseif( f(i,j) == f_min .and.  f(i+1,j) /= f_min .and. f(i+1,j+1) == f_min .and. f(i,j+1) == f_min)then
            ! a northwest corner
                dfdx_2d_stag_os = ( f(i+2,j) - f(i+1,j) ) / (delta)
            elseif( f(i,j) == f_min .and.  f(i+1,j) == f_min .and. f(i+1,j+1) /= f_min .and. f(i,j+1) == f_min)then
            ! a southwest corner
                dfdx_2d_stag_os = ( f(i+2,j+1) - f(i+1,j+1) ) / (delta)
            elseif( f(i,j) == f_min .and.  f(i+1,j) == f_min .and. f(i+1,j+1) == f_min .and. f(i,j+1) /= f_min)then
            ! a southeast corner
                dfdx_2d_stag_os = ( f(i,j+1) - f(i-1,j+1) ) / (delta)
            endif

        endif

    end function dfdx_2d_stag_os

!----------------------------------------------------------------------------

    function dfdy_2d_stag_os(f_in, i, j, delta, thck, thklim )

    !*SFP* altered/expanded version of above function that uses approx. one-sided
    ! diffs at physical domain edges so as not to overesimate grads there.

        implicit none
        real(dp), dimension(:,:), intent(in) :: f_in, thck
        integer, intent(in) :: i, j
        real(dp), intent(in) :: delta, thklim
        real(dp) :: dfdy_2d_stag_os

        real(dp), dimension(2,2) :: f_array, thck_array
        real(dp), dimension(1:size(f_in,1),1:size(f_in,2)) :: f
        real(dp) :: f_min

        ! initialize to zeros
        dfdy_2d_stag_os = 0.0d0; f_array = 0.0d0; f_min = 0.0d0

        f = f_in

        where( thck <= thklim )
            f = 0
        end where

        f_array(1,1) = f(i,j); f_array(2,1) = f(i+1,j); f_array(1,2) = f(i,j+1); f_array(2,2) = f(i+1,j+1);

        !TODO - Change reals to dp
        if( sum( f_array/ f_array, MASK = f_array /= 0.0d0 ) == 4.0 )then

            ! normal differencing for interior points
            dfdy_2d_stag_os = (f(i,j+1) + f(i+1,j+1) - f(i,j) - f(i+1,j))/(2*delta)

        elseif( sum( f_array/ f_array, MASK = f_array /= 0.0d0 ) == 3.0 )then

            ! corner; use 2x next closest value
            if( f(i,j) == f_min )then     ! southwest corner point missing: apply value from s.e. point
                dfdy_2d_stag_os = (f(i,j+1) + f(i+1,j+1) - 2.0*f(i+1,j))/(2*delta)
            elseif( f(i+1,j) == f_min )then ! southeast corner point missing: apply value from s.w. point
                dfdy_2d_stag_os = (f(i,j+1) + f(i+1,j+1) - 2.0*f(i,j))/(2*delta)
            elseif( f(i,j+1) == f_min )then ! northwest corner point missing: apply value from n.e. point
                dfdy_2d_stag_os = ( 2.0*f(i+1,j+1) - f(i,j) - f(i+1,j))/(2*delta)
            elseif( f(i+1,j+1) == f_min )then ! northeast corner point missing: apply value from n.w. point
                dfdy_2d_stag_os = ( 2.0*f(i,j+1) - f(i,j) - f(i+1,j))/(2*delta)
            endif

        elseif( sum( f_array/ f_array, MASK = f_array /= 0.0d0 ) == 2.0 )then

            ! side; back up and take gradient from points one set of cells in OR use only the single set of
            ! cells available along the differencing direction 
            if( f(i,j) == f_min .and. f(i,j+1) == f_min )then   ! west cells empty
                dfdy_2d_stag_os = (f(i+1,j+1) - f(i+1, j))/(delta)
            elseif( f(i+1,j) == f_min .and. f(i+1,j+1) == f_min )then   ! east cells empty
                dfdy_2d_stag_os = (f(i,j+1) - f(i,j) )/(delta)
            elseif( f(i,j+1) == f_min .and. f(i+1,j+1) == f_min )then   ! north cells empty
                dfdy_2d_stag_os = (f(i,j) + f(i+1,j) - f(i,j-1) - f(i+1,j-1))/(2*delta)
            elseif( f(i,j) == f_min .and. f(i+1,j) == f_min )then   ! south cells empty
                dfdy_2d_stag_os = (f(i,j+2) + f(i+1,j+2) - f(i,j+1) - f(i+1,j+1))/(2*delta)
            endif

        elseif( sum( f_array/ f_array, MASK = f_array /= 0.0d0 ) == 1.0 )then

            ! isolated; treat by assuming it is part of a 3 block for which the rest of the values are not contained within
            ! the local 2x2 block with indices i:i+1, j:j+1 
            if( f(i,j) /= f_min .and.  f(i+1,j) == f_min .and. f(i+1,j+1) == f_min .and. f(i,j+1) == f_min )then
            ! a northeast corner
                dfdy_2d_stag_os = ( f(i,j) - f(i,j-1) ) / (delta)
            elseif( f(i,j) == f_min .and.  f(i+1,j) /= f_min .and. f(i+1,j+1) == f_min .and. f(i,j+1) == f_min )then
            ! a northwest corner
                dfdy_2d_stag_os = ( f(i+1,j) - f(i+1,j-1) ) / (delta)
            elseif( f(i,j) == f_min .and.  f(i+1,j) == f_min .and. f(i+1,j+1) /= f_min .and. f(i,j+1) == f_min )then
            ! a southwest corner
                dfdy_2d_stag_os = ( f(i+1,j+2) - f(i+1,j+1) ) / (delta)
            elseif( f(i,j) == f_min .and.  f(i+1,j) == f_min .and. f(i+1,j+1) == f_min .and. f(i,j+1) /= f_min )then
            ! a southeast corner
                dfdy_2d_stag_os = ( f(i,j+2) - f(i,j+1) ) / (delta)
            endif

        endif

    end function dfdy_2d_stag_os

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to x at the given point
    !*FD using an upwind method (suitable for maximum boundaries)

    function dfdx_2d_upwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_2d_upwind
        dfdx_2d_upwind = (.5 * f(i-2,j) - 2 * f(i-1, j) + 1.5 * f(i, j))/delta
    end function dfdx_2d_upwind

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to y at the given point
    !*FD using an upwind method (suitable for maximum boundaries)

    function dfdy_2d_upwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_2d_upwind
        dfdy_2d_upwind = (.5 * f(i,j-2) - 2 * f(i, j-1) + 1.5 * f(i, j))/delta
    end function dfdy_2d_upwind

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to x at the given point
    !*FD using a downwind method (suitable for minimum boundaries)

    function dfdx_2d_downwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_2d_downwind
        dfdx_2d_downwind = (-1.5 * f(i, j) + 2 * f(i+1, j) - .5 * f(i+2, j))/delta
    end function dfdx_2d_downwind 

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to y at the given point
    !*FD using a downwind method (suitable for minimum boundaries)

    function dfdy_2d_downwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_2d_downwind
        dfdy_2d_downwind = (-1.5 * f(i, j) + 2 * f(i, j+1) - .5 * f(i, j+2))/delta
    end function dfdy_2d_downwind

!----------------------------------------------------------------------------

    !------------------------------------------------------------------
    !First Derivative Estimates, Second Order, 3D
    !------------------------------------------------------------------

    !*FD Computes derivative with respect to x at a given point

    function dfdx_3d(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_3d
        dfdx_3d = (-.5/delta)*f(k, i-1, j)  + (.5/delta)*f(k, i+1, j)
    end function dfdx_3d

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to y at a given point

    function dfdy_3d(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_3d
        dfdy_3d = (-.5/delta)*f(k, i, j-1) + (.5/delta)*f(k, i, j+1)
    end function dfdy_3d
    
!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to z at a given point
    !*FD where the Z axis uses an irregular grid defined by \ittext{deltas}.
    !*FD This derivative is given by the formula:

    function dfdz_3d_irregular(f, i, j, k, dz)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), dimension(:), intent(in) :: dz
        real(dp) :: dfdz_3d_irregular

        dfdz_3d_irregular = f(k-1,i,j)*(dz(k) - dz(k+1))/((dz(k) - dz(k-1))*(dz(k+1)-dz(k-1))) + &
                            f(k,  i,j)*(dz(k+1)-2*dz(k)+dz(k-1))/((dz(k)-dz(k-1))*(dz(k+1)-dz(k))) + &
                            f(k+1,i,j)*(dz(k)-dz(k-1))/((dz(k+1)-dz(k))*(dz(K+1)-dz(k-1)))
    end function
    
!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to z at a given point using an upwinding
    !*FD scheme.  The Z axis uses an irregular grid defined by \iittext{deltas}.

    function dfdz_3d_upwind_irregular(f, i, j, k, deltas)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), dimension(:), intent(in) :: deltas
        real(dp) :: dfdz_3d_upwind_irregular
        real(dp) :: zkMinusZkm1, zkMinusZkm2, zkm1MinusZkm2
        zkMinusZkm1 = deltas(k) - deltas(k-1)
        zkMinusZkm2 = deltas(k) - deltas(k-2)
        zkm1MinusZkm2 = deltas(k-1) - deltas(k-2)
        
        dfdz_3d_upwind_irregular = f(k-2, i, j) * zkMinusZkm1 / (zkm1MinusZkm2 * zkMinusZkm2) - &
                                   f(k-1, i, j) * zkMinusZkm2 / (zkMinusZkm1 * zkm1MinusZkm2) + &
                                   f(k,   i, j) * (2*deltas(k) - deltas(k-1) - deltas(k-2)) / (zkMinusZkm1 * zkMinusZkm2)
    end function
    
!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to z at a given point using a downwinding
    !*FD scheme.  The Z axis uses an irregular grid defined by \iittext{deltas}.

    function dfdz_3d_downwind_irregular(f, i, j, k, deltas)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), dimension(:), intent(in) :: deltas
        real(dp) :: dfdz_3d_downwind_irregular
        real(dp) :: zkp1MinusZk, zkp2MinusZk, zkp2MinusZkp1
        zkp1MinusZk = deltas(k+1) - deltas(k)
        zkp2MinusZk = deltas(k+2) - deltas(k)
        zkp2MinusZkp1 = deltas(k+2) - deltas(k+1)
        
        dfdz_3d_downwind_irregular =f(k,   i, j) * (-zkp1MinusZk - zkp2MinusZk)/(zkp1MinusZk * zkp2MinusZk) + &
                                    f(k+1, i, j) * zkp2MinusZk / (zkp2MinusZkp1 * zkp1MinusZk) - &
                                    f(k+2, i, j) * zkp1MinusZk / (zkp2MinusZkp1 * zkp2MinusZk)
    end function
    
!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to x at the equivalent
    !*FD point on a staggered grid.

    function dfdx_3d_stag(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_3d_stag
        dfdx_3d_stag = (f(k, i+1, j) + f(k, i+1, j+1) - f(k, i, j) - f(k, i, j+1))/(2*delta) 
    end function dfdx_3d_stag

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to y at the equivalent
    !*FD point on a staggered grid.

    function dfdy_3d_stag(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_3d_stag
        dfdy_3d_stag = (f(k, i, j+1) + f(k, i+1, j+1) - f(k, i, j) - f(k, i+1, j))/(2*delta)
    end function dfdy_3d_stag

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to x at the given point
    !*FD using an upwind method (suitable for maximum boundaries)

    function dfdx_3d_upwind(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_3d_upwind
        dfdx_3d_upwind = (.5 * f(k, i-2, j) - 2 * f(k, i-1, j) + 1.5 * f(k, i, j))/delta
    end function dfdx_3d_upwind

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to y at the given point
    !*FD using an upwind method (suitable for maximum boundaries)

    function dfdy_3d_upwind(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_3d_upwind
        dfdy_3d_upwind = (.5 * f(k, i, j-2) - 2 * f(k, i, j-1) + 1.5 * f(k, i, j))/delta
    end function dfdy_3d_upwind

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to x at the given point
    !*FD using a downwind method (suitable for minimum boundaries)

    function dfdx_3d_downwind(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j, k
        real(dp), intent(in) :: delta
        real(dp) :: dfdx_3d_downwind
        dfdx_3d_downwind = (-1.5 * f(k, i, j) + 2 * f(k, i+1, j) - .5 * f(k, i+2, j))/delta
    end function dfdx_3d_downwind 

!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to y at the given point
    !*FD using a downwind method (suitable for minimum boundaries)

    function dfdy_3d_downwind(f, i, j, k, delta)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta
        real(dp) :: dfdy_3d_downwind
        dfdy_3d_downwind = (-1.5 * f(k, i, j) + 2 * f(k, i, j+1) - .5 * f(k, i, j+2))/delta
    end function dfdy_3d_downwind
    
!----------------------------------------------------------------------------

    !------------------------------------------------------------------
    !Second Derivative Estimates, Second Order
    !------------------------------------------------------------------
    
    !*FD Computes 2nd derivative with respect to x at the given point

    function d2fdx2_2d(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdx2_2d     
        d2fdx2_2d = (f(i+1,j) + f(i-1,j) - 2 * f(i, j))/(delta*delta)
    end function d2fdx2_2d

!----------------------------------------------------------------------------
    
    function d2fdx2_2d_downwind(f,i,j,delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdx2_2d_downwind   

        d2fdx2_2d_downwind = (3*f(i, j) - 7*f(i+1, j) + 5*f(i+2, j) - f(i+3, j)) / (2*delta**2)

    end function d2fdx2_2d_downwind

!----------------------------------------------------------------------------

    function d2fdx2_2d_upwind(f,i,j,delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdx2_2d_upwind 

        d2fdx2_2d_upwind = (3*f(i, j) - 7*f(i-1, j) + 5*f(i-2, j) - f(i-3, j)) / (2*delta**2)

    end function d2fdx2_2d_upwind

!----------------------------------------------------------------------------

    function d2fdy2_2d_downwind(f,i,j,delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdy2_2d_downwind   

        d2fdy2_2d_downwind = (3*f(i, j) - 7*f(i, j+1) + 5*f(i, j+2) - f(i, j+3)) / (2*delta**2)

    end function d2fdy2_2d_downwind

!----------------------------------------------------------------------------

    function d2fdy2_2d_upwind(f,i,j,delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdy2_2d_upwind 

        d2fdy2_2d_upwind = (3*f(i, j) - 7*f(i, j-1) + 5*f(i, j-2) - f(i, j-3)) / (2*delta**2)

    end function d2fdy2_2d_upwind

!----------------------------------------------------------------------------

    function d2fdx2_2d_stag(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdx2_2d_stag
        
        !This formula can be derived using two central differences
        !(i to i+2, and i-1 to i+1) to get the derivative at
        !i and i+1, then applying a central difference to that
        !in order to get the 2nd derivative at a staggered point
        d2fdx2_2d_stag = sum(f(i+2, j:j+1) + f(i-1, j:j+1) - f(i+1, j:j+1) - f(i, j:j+1))/(4*delta**2)
    end function d2fdx2_2d_stag

!----------------------------------------------------------------------------

    function d2fdx2_2d_stag_downwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdx2_2d_stag_downwind
        
        d2fdx2_2d_stag_downwind = sum(3*f(i, j:j+1) - 7*f(i+1, j:j+1) + 5*f(i+2, j:j+1) - f(i+3, j:j+1)) / (4*delta**2)
        end function d2fdx2_2d_stag_downwind
        
        function d2fdx2_2d_stag_upwind(f, i, j, delta)
            implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdx2_2d_stag_upwind
        
        d2fdx2_2d_stag_upwind = sum(-3*f(i+1, j:j+1) + 7*f(i, j:j+1) - 5*f(i-1, j:j+1) + f(i-2, j:j+1)) / (4*delta**2)
    end function d2fdx2_2d_stag_upwind

!----------------------------------------------------------------------------

    !*FD Computes 2nd derivative with respect to y at the given point

    function d2fdy2_2d(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdy2_2d
        d2fdy2_2d = (f(i, j+1) + f(i, j-1) - 2 * f(i, j))/(delta*delta)
    end function d2fdy2_2d
    
!----------------------------------------------------------------------------

    function d2fdy2_2d_stag(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdy2_2d_stag
        
        !This formula can be derived using two central differences
        !(i to i+2, and i-1 to i+1) to get the derivative at
        !i and i+1, then applying a central difference to that
        !in order to get the 2nd derivative at a staggered point
        d2fdy2_2d_stag = sum(f(i:i+1, j+2) + f(i:i+1, j-1) - f(i:i+1, j+1) - f(i:i+1, j))/(4*delta**2)
    end function d2fdy2_2d_stag
    
!----------------------------------------------------------------------------

    function d2fdy2_2d_stag_downwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdy2_2d_stag_downwind
        
        d2fdy2_2d_stag_downwind = sum(3*f(i:i+1, j) - 7*f(i:i+1, j+1) + 5*f(i:i+1, j+2) - f(i:i+1, j+3)) / (4*delta**2)
    end function d2fdy2_2d_stag_downwind
        
!----------------------------------------------------------------------------

    function d2fdy2_2d_stag_upwind(f, i, j, delta)
        implicit none
        real(dp), dimension(:,:), intent(in) :: f
        integer, intent(in) :: i,j
        real(dp), intent(in) :: delta
        real(dp) :: d2fdy2_2d_stag_upwind
        
        d2fdy2_2d_stag_upwind = sum(-3*f(i:i+1, j+1) + 7*f(i:i+1, j) - 5*f(i:i+1, j-1) + f(i:i+1, j-2)) / (4*delta**2)
    end function d2fdy2_2d_stag_upwind
   
!----------------------------------------------------------------------------

    subroutine d2f_field(f, deltax, deltay, d2fdx2, d2fdy2, direction_x, direction_y)

        use parallel
!!        use glimmer_horiz_bcs, only: horiz_bcs_unstag_scalar
        implicit none 

        real(dp), intent(out), dimension(:,:) :: d2fdx2, d2fdy2
        real(dp), intent(in), dimension(:,:) :: f    ! unstaggered grid
        real(dp), intent(in) :: deltax, deltay
        real(dp), intent(in), dimension(:,:), optional :: direction_x, direction_y
        integer :: i,j

!LOOP: all scalar points
!      uses upwinding and downwinding to avoid going out of bounds

        do i = 1,size(f,1)
            do j = 1,size(f,2)

                !unstaggered grid
                if (i == 1) then
                    d2fdx2(i,j) = d2fdx2_2d_downwind(f,i,j,deltax)
                else if (i == size(f,1)) then
                    d2fdx2(i,j) = d2fdx2_2d_upwind(f,i,j,deltax)
                else
                    if (present(direction_x)) then
                        if (direction_x(i,j) > 0) then
                            d2fdx2(i,j) = d2fdx2_2d_downwind(f,i,j,deltax)
                        else if (direction_x(i,j) < 0) then
                            d2fdx2(i,j) = d2fdx2_2d_upwind(f,i,j,deltax)
                        else
                            d2fdx2(i,j) = d2fdx2_2d(f,i,j,deltax)
                        end if
                    else
                        d2fdx2(i,j) = d2fdx2_2d(f,i,j,deltax)
                    end if
                end if
                
                if (j == 1) then
                    d2fdy2(i,j) = d2fdy2_2d_downwind(f,i,j,deltax)
                else if (j == size(f,2)) then
                    d2fdy2(i,j) = d2fdy2_2d_upwind(f,i,j,deltax)
                else
                    if (present(direction_y)) then
                        if (direction_y(i,j) > 0) then
                            d2fdy2(i,j) = d2fdy2_2d_downwind(f,i,j,deltax)
                        else if (direction_y(i,j) < 0) then
                            d2fdy2(i,j) = d2fdy2_2d_upwind(f,i,j,deltax)
                        else
                            d2fdy2(i,j) = d2fdy2_2d(f,i,j,deltax)
                        end if
                    else
                        d2fdy2(i,j) = d2fdy2_2d(f,i,j,deltax)
                    end if
                end if
            end do
        end do

!TODO:  If halo updates are needed, they should be done at a higher level.

        call parallel_halo(d2fdx2)
!        call horiz_bcs_unstag_scalar(d2fdx2)
        call parallel_halo(d2fdy2)
!        call horiz_bcs_unstag_scalar(d2fdy2)

    end subroutine d2f_field

!----------------------------------------------------------------------------

!TODO: Rewrite this using the existing derivative machinery?

    subroutine d2f_field_stag(f, deltax, deltay, d2fdx2, d2fdy2, periodic_x, periodic_y)

    implicit none 

    real(dp), intent(out), dimension(:,:) :: d2fdx2, d2fdy2
    real(dp), intent(in), dimension(:,:) :: f
    real(dp), intent(in) :: deltax, deltay
    logical :: periodic_x, periodic_y

    real(dp) :: dewsq4, dnssq4
    integer :: ew,ns

    integer :: pt(2)
    integer :: nsn
    integer :: ewn
    
    nsn = size(f,2)
    ewn = size(f,1)

    dewsq4 = 4.0d0 * deltax * deltax
    dnssq4 = 4.0d0 * deltay * deltay

    d2fdx2 = 0.d0
    d2fdy2 = 0.d0
 
!LOOP - not sure what bounds should be in this subroutine

    do ns = 2, nsn-2
      do ew = 2, ewn-2
        d2fdx2(ew,ns) = centerew(ew,ns)
        d2fdy2(ew,ns) = centerns(ew,ns)
      end do
    end do

! *** 2nd order boundaries using upwinding

    do ew = 1, ewn-1, ewn-2

      pt = whichway(ew)

      do ns = 2, nsn-2 
          d2fdx2(ew,ns) = boundyew(pt,ns)
          d2fdy2(ew,ns) = centerns(ew,ns)
      end do

    end do

    do ns = 1, nsn-1, nsn-2

      pt = whichway(ns)

      do ew = 2, ewn-2
          d2fdx2(ew,ns) = centerew(ew,ns)
          d2fdy2(ew,ns) = boundyns(pt,ew)
      end do

    end do

    do ns = 1, nsn-1, nsn-2
      do ew = 1, ewn-1, ewn-2
          pt = whichway(ew)
          d2fdx2(ew,ns) = boundyew(pt,ns)
          pt = whichway(ns)
          d2fdy2(ew,ns) = boundyns(pt,ew)
      end do
    end do

  contains

!----------------------------------------------------------------------------

    function centerew(ew,ns)

      implicit none

      real(dp) :: centerew
      integer ns,ew

      centerew = (sum(f(ew+2,ns:ns+1)) + sum(f(ew-1,ns:ns+1)) - &
                  sum(f(ew+1,ns:ns+1)) - sum(f(ew,ns:ns+1))) / dewsq4
    
    end function centerew 

!----------------------------------------------------------------------------

    function centerns(ew,ns)

      implicit none

      real(dp) :: centerns
      integer ns,ew

      centerns = (sum(f(ew:ew+1,ns+2)) + sum(f(ew:ew+1,ns-1)) - &
                  sum(f(ew:ew+1,ns+1)) - sum(f(ew:ew+1,ns))) / dnssq4
  
    end function centerns 

!----------------------------------------------------------------------------

    function boundyew(pt,ns)

      implicit none

      integer, intent(in) :: pt(2)
      real(dp) :: boundyew
      integer ns

      boundyew = pt(1) * (3.0d0 * sum(f(pt(2),ns:ns+1)) - 7.0d0 * sum(f(pt(2)+pt(1),ns:ns+1)) + &
                 5.0d0 * sum(f(pt(2)+2*pt(1),ns:ns+1)) - sum(f(pt(2)+3*pt(1),ns:ns+1))) / dewsq4

    end function boundyew

!----------------------------------------------------------------------------

    function boundyns(pt,ew)

      implicit none

      integer, intent(in) :: pt(2)
      real(dp) :: boundyns
      integer ew

      boundyns = pt(1) * (3.0d0 * sum(f(ew:ew+1,pt(2))) - 7.0d0 * sum(f(ew:ew+1,pt(2)+pt(1))) + &
                 5.0d0 * sum(f(ew:ew+1,pt(2)+2*pt(1))) - sum(f(ew:ew+1,pt(2)+3*pt(1)))) / dnssq4

    end function boundyns

!----------------------------------------------------------------------------

    function whichway(i)

      implicit none

      integer, intent(in) :: i
      integer :: whichway(2) 

      if (i == 1) then 
        whichway = (/1,1/)
      else
        whichway = (/-1,i+1/)
      end if

    end function whichway

!----------------------------------------------------------------------------
!TODO: Not sure what this code is doing here.

!       real(dp), dimension(:,:), intent(in) :: f
!       real(dp), dimension(:,:), intent(out) :: d2fdx2, d2fdy2
!       real(dp), intent(in) :: deltax, deltay
!       logical :: periodic_x, periodic_y
!       
!       integer :: nx, x, ny, y
!       
!       nx = size(f, 1)
!       ny = size(f, 2)
!       
!       !NOTE: See the field 1st derivative staggered function for
!       !a discussion of periodic boundary conditions
!       
!       !First compute the values that do not fall on any boundaries
!       !This is the same regardless of whether periodic boundary
!       !conditions are used
!       do x = 1, nx-1
!             do y = 1, ny-1
!                 if (x == 1) then
!                     d2fdx2(1,y) = d2fdx2_2d_stag_downwind(f, 1, y, deltax)
!                 else if (x == nx - 1) then
!                     d2fdx2(nx-1, y) = d2fdx2_2d_stag_upwind(f, nx-1, y, deltax)
!                 else
!                     d2fdx2(x,y) = d2fdx2_2d_stag(f, x, y, deltax)
!                 end if
!                 
!                 if (y == 1) then
!                     d2fdy2(x,1) = d2fdy2_2d_stag_downwind(f, x, 1, deltay)
!                 else if (y == ny - 1) then
!                     d2fdy2(x, ny-1) = d2fdy2_2d_stag_upwind(f, x, ny-1, deltay)
!                 else
!                     d2fdy2(x,y) = d2fdy2_2d_stag(f, x, y, deltay)
!                 end if
!               end do
!       end do
!       
!       !If we are not using periodic boundary conditions, then we need
!       !to use an upwinding scheme to get the values when x = 1, y = 1,
!       !x = nx - 1, or y = ny - 1
!               !If we are using periodic boundary conditions, then compute the
!       !boundaries with input from the periodic conditions.  We do not
!       !upwind or downwind.  Also, because an extra set of values around
!       !the edges is necessary to correctly maintain periodicity,
!       !we fill in values where x = nx and where y = ny (whereas we
!       !do not with nonperiodic boundaries, as the staggered grid
!       !points fall strictly in the interior of the nonstaggered
!       !grid)
!       do y = 1, ny - 2
!               if (.not.(periodic_x)) then
!                       d2fdx2(1,y) = d2fdx2_2d_stag_downwind(f, 1, y, deltax)
! 
!                       d2fdx2(nx-1, y) = d2fdx2_2d_stag_upwind(f, nx-1, y, deltax)
!                       
!                   else
!                       !Because of the periodicity, I will simply copy the appropriate values
!                       !(e.g. u(1) = u(n-2), u(n-1) = u(2)
!                       d2fdx2(1,y) = d2fdx2(nx-2,y)
!                       d2fdx2(nx-1,y) = d2fdx2(2,y)
!               end if
!               d2fdy2(1,y) = d2fdy2_2d_stag(f, 1, y, deltay)    
!               d2fdy2(nx-1, y) = d2fdy2_2d_stag(f, nx-1, y, deltay)    
!       end do
!       
!       !See comments for the periodic x boundary case above; the same
!       !principles apply here.
!       do x=1, nx-2
!               if (.not.(periodic_y)) then
!                       d2fdy2(x,1) = d2fdy2_2d_stag_downwind(f, x, 1, deltay)
!                       d2fdy2(x, ny-1) = d2fdy2_2d_stag_upwind(f, x, ny-1, deltay)
!               else
!                               d2fdy2(x,1) = d2fdy2(x,ny-2)
!                               d2fdy2(x,nx-1) = d2fdy2(x,2)
!               end if
!               d2fdx2(x,1) = d2fdx2_2d_stag(f, x, 1, deltax)
!               d2fdx2(x, ny-1) = d2fdx2_2d_stag(f, x, ny-1, deltax)
!       end do
!       
!       
!       !TODO: Change this to use the scheme above
!       !We have neglected so far to take care of the four points that occur at the
!       !maximum two indices in x and y.  If no periodic boundaries are being used,
!       !we compute the value zt (nx-1, ny-1) using upwinding schemes.
!       if (.not. periodic_x .and. .not. periodic_y) then
!               d2fdx2(nx-1, ny-1) = d2fdx2_2d_stag_upwind(f, nx-1, ny-1, deltax)
!               d2fdy2(nx-1, ny-1) = d2fdy2_2d_stag_upwind(f, nx-1, ny-1, deltay)
!       else if (.not. periodic_x) then
!               !y is periodic - this means we need to compute the derivative
!               !for x=nx-1 and y=ny, ny-1.  We will copy and paste
!               !y derivatives (for now), as above, and upwind
!               !the x derivatives
!               d2fdx2(nx-1, ny-1) = d2fdx2_2d_stag_upwind(f, nx-1, ny-1, deltax)
!               d2fdy2(nx-1, ny-1) = sum(f(nx-1:nx, 1) + f(nx-1:nx, ny-2) - f(nx-1:nx, ny) - f(nx-1:nx, ny-1))/(4*deltay**2)
!               
!               
!               d2fdx2(nx-1, ny)  = d2fdx2_2d_stag_upwind(f, nx-1, ny, deltax)
!               d2fdy2(nx-1, ny)  = sum(f(nx-1:nx, 2) + f(nx-1:nx, ny-1) - f(nx-1:nx, 1) - f(nx-1:nx, ny))/(4*deltay**2)
!               
!       else if (.not. periodic_y) then
!               !See comments for the periodic y case above - we are basically using the same
!               !logic with x and y swapped
!               d2fdx2(nx-1, ny-1) = sum(f(1, ny-1:ny) + f(nx-2, ny-1:ny) - f(nx, ny-1:ny) - f(nx-1, ny-1:ny))/(4*deltax**2)
!               d2fdy2(nx-1, ny-1) = d2fdy2_2d_stag_upwind(f, nx-1, ny-1, deltay)
!               
!               d2fdx2(nx, ny-1) = sum(f(2, ny-1:ny) + f(nx-1, ny-1:ny) - f(1, ny-1:ny) - f(nx, ny-1:ny))/(4*deltax**2)
!               d2fdy2(nx, ny-1) = d2fdy2_2d_stag_upwind(f, nx-1, ny-1, deltay)
!       else
!               !X and Y are periodic; we will use the periodic forms of the above differences
!               !Some of these will get very funky because 
!                       d2fdx2(nx-1, ny-1) = sum(f(1, ny-1:ny) + f(nx-1, ny-1:ny) - f(nx, ny-1:ny) - f(nx-1, ny-1:ny))/(4*deltax**2)
!                       d2fdy2(nx-1, ny-1) = sum(f(nx-1:nx, 1) + f(nx-1:nx, ny-2) - f(nx-1:nx, ny) - f(nx-1:nx, ny-1))/(4*deltay**2)
!               
!               d2fdx2(nx, ny-1) = sum(f(2, ny-1:ny) + f(nx-1, ny-1:ny) - f(1, ny-1:ny) - f(nx, ny-1:ny))/(4*deltax**2)
!               d2fdy2(nx, ny-1) = ((f(nx, 1) + f(nx, ny-2) - f(nx, ny) - f(nx, ny-1)) + &
!                                    (f(1, 1) + f(1, ny-2) - f(1, ny) - f(1, ny-1)))/(4*deltay**2)
!               
!               d2fdy2(nx-1, ny)  = ((f(1, ny) + f(nx-1, ny) - f(nx, ny) - f(nx-1, ny)) + &
!                                    (f(1, 1) + f(nx-1, 1) - f(nx, 1) - f(nx-1, 1)))/(4*deltax**2)
!                       d2fdy2(nx-1, ny)  = sum(f(nx-1:nx, 2) + f(nx-1:nx, ny-1) - f(nx-1:nx, 1) - f(nx-1:nx, ny))/(4*deltay**2)
!                       
!                       d2fdx2(nx, ny) = ((f(2, ny) + f(nx-1, ny) - f(1, ny) - f(nx, ny)) + (f(2, 1) + f(nx-1, 1) - f(1, 1) - f(nx, 1))) / (4*deltax**2)
!                       d2fdy2(nx, ny) = ((f(nx, 2) + f(nx, ny-1) - f(nx, 1) - f(nx, ny)) + (f(1, 2) + f(1, ny-1) - f(1, 1) - f(1, ny)))/(4*deltay**2)
!               
!       end if
        
    end subroutine d2f_field_stag
    
!----------------------------------------------------------------------------

    !*FD Computes derivative taken first w.r.t x, then to y at the given point.

    function d2fdxy_3d(f, i, j, k, delta_x, delta_y)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta_x, delta_y
        real(dp) :: d2fdxy_3d
        
        d2fdxy_3d = (f(k, i-1, j-1) - f(k, i-1, j+1) - f(k, i+1, j-1) + f(k, i+1, j+1))/(4*delta_x*delta_y) 
    end function d2fdxy_3d
    
!----------------------------------------------------------------------------

    function d2fdxz_3d(f, i, j, k, delta_x, dz)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta_x
        real(dp), dimension(:), intent(in) :: dz
        real(dp) :: d2fdxz_3d
        
        d2fdxz_3d = (.5/delta_x) * ( &
          (f(k-1, i+1, j) - f(k-1, i-1, j)) * (dz(k) - dz(k+1)) / ( (dz(k) - dz(k-1)) * (dz(k+1) - dz(k-1)) ) + &
          (f(k,   i+1, j) - f(k,   i-1, j)) * (dz(k+1) + dz(k-1) - 2*dz(k)) / ( (dz(k) - dz(k-1)) * (dz(k+1) - dz(k)) ) + &
          (f(k+1, i+1, j) - f(k+1, i-1, j)) * (dz(k) - dz(k-1)) / ( (dz(k+1) - dz(k)) * (dz(k+1) - dz(k-1)) ) )
        end function d2fdxz_3d
        
        function d2fdyz_3d(f, i, j, k, delta_x, dz)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), intent(in) :: delta_x
        real(dp), dimension(:), intent(in) :: dz
        real(dp) :: d2fdyz_3d
        
        d2fdyz_3d = (.5/delta_x) * ( &
          (f(k-1, i, j+1) - f(k-1, i, j-1)) * (dz(k) - dz(k+1)) / ( (dz(k) - dz(k-1)) * (dz(k+1) - dz(k-1)) ) + &
          (f(k,   i, j+1) - f(k,   i, j-1)) * (dz(k+1) + dz(k-1) - 2*dz(k)) / ( (dz(k) - dz(k-1)) * (dz(k+1) - dz(k)) ) + &
          (f(k+1, i, j+1) - f(k+1, i, j-1)) * (dz(k) - dz(k-1)) / ( (dz(k+1) - dz(k)) * (dz(k+1) - dz(k-1)) ) )
        end function d2fdyz_3d
        
!----------------------------------------------------------------------------

    !*FD Computes derivative with respect to z at a given point
    !*FD where the Z axis uses an irregular grid defined by \ittext{deltas}.

    function d2fdz2_3d_irregular(f, i, j, k, deltas)
        implicit none
        real(dp), dimension(:,:,:), intent(in) :: f
        integer, intent(in) :: i,j,k
        real(dp), dimension(:), intent(in) :: deltas
        real(dp) :: d2fdz2_3d_irregular
        real(dp) :: zkMinusZkp1, zkMinusZkm1, zkp1MinusZkm1, zkp1MinusZk
        
        zkMinusZkp1 = deltas(k) - deltas(k+1)
        zkMinusZkm1 = deltas(k) - deltas(k-1)
        zkp1MinusZkm1 = deltas(k+1) - deltas(k-1)
        zkp1MinusZk = -1 * zkMinusZkp1
        
        
        d2fdz2_3d_irregular = 2 * f(k-1, i, j) / (zkMinusZkm1 * zkp1MinusZkm1) - &
                              2 * f(k,   i, j) / (zkp1MinusZk * zkMinusZkm1) + &
                              2 * f(k+1, i, j) / (zkp1Minuszk * zkp1MinusZkm1)    
    end function d2fdz2_3d_irregular

!---------------------------------------------------------------------------------

end module glam_grid_operators

!----------------------------------------------------------------------------
